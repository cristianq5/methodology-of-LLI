# -*- coding: utf-8 -*-
# Test sintético de las métricas de copa (lli_core.R: metricas_copa y perfil_anchura) calculadas A MANO
# (réplica de python/pipeline/tests/test_metricas_copa.py: mismo caso y mismos valores esperados).
#
# Franjas de 0,1 m desde 0,3 m (franja j = [0,3 + 0,1 j, 0,4 + 0,1 j)); puntos a la altura central de cada franja y repartidos
# en t de -w/2 a +w/2 (anchura estricta = w).
#   Rebanada 0 (s = 0,05): franjas 2..7 con anchuras 0,6 0,7 0,8 0,9 1,0 0,5 (10 puntos por franja, 60 en total).
#       área = 0,1 x 4,5 = 0,45 m² ; volumen = 0,045 m³ ; base de copa = 0,5 m ; MH = 1,05 ; profundidad = 0,55
#       anchura máxima = 1,0 en la franja 6 (centro 0,95 m) ; MW = 4,5 / 6 = 0,75 ; TRV (4 m) = 1,05 x 0,75 x 10 000 / 4 = 1968,75
#   Rebanada 1 (s = 0,15): franjas 0, 1, 3 y 4 con anchura 1,0 (15 puntos por franja); la franja 2 está vacía (hueco).
#       área = 0,4 ; base = 0,3 ; MH = 0,75 ; profundidad = 0,45 ; anchura máxima 1,0 en la franja 0 (centro 0,35 m)
#   Perfil de la hilera (2 rebanadas): franja 0 -> ocupada en 1 de 2, anchura 1,0, contando huecos 0,5;
#       franja 2 -> solo la rebanada 0 (0,6), contando huecos 0,3; franja 3 -> (0,7 + 1,0) / 2 = 0,85.
# Uso:  Rscript R/pipeline/tests/test_metricas_copa.R
a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
DIR_SCRIPTS <- dirname(dirname(normalizePath(sub("^--file=", "", a[1]))))
source(file.path(DIR_SCRIPTS, "utils.R")); source(file.path(DIR_SCRIPTS, "lli_core.R"))

cerca <- function(x, y, tol = 1e-8) length(x) == length(y) && all(is.finite(x)) && all(abs(x - y) <= tol + 1e-5 * abs(y))
comprobar <- function(cond, msg) if (!isTRUE(cond)) stop("FALLA: ", msg, call. = FALSE)

construir <- function() {
  s <- c(); t <- c(); h <- c()
  for (k in 1:6) {
    j <- (2:7)[k]; w <- c(0.6, 0.7, 0.8, 0.9, 1.0, 0.5)[k]
    tt <- seq(-w / 2, w / 2, length.out = 10)
    s <- c(s, rep(0.05, 10)); t <- c(t, tt); h <- c(h, rep(0.3 + 0.1 * j + 0.05, 10))
  }
  for (j in c(0, 1, 3, 4)) {
    tt <- seq(-0.5, 0.5, length.out = 15)
    s <- c(s, rep(0.15, 15)); t <- c(t, tt); h <- c(h, rep(0.3 + 0.1 * j + 0.05, 15))
  }
  list(s = s, t = t, h = h)
}

d <- construir()
P <- perfiles_rebanada(d$s, d$t, d$h, s_ini = 0.0, n_reb = 2, min_pts_franja = 5, min_pts_reb = 50, n_franjas = 10)
m <- metricas_copa(P, 4.0)
comprobar(cerca(m$area_seccion_m2, c(0.45, 0.4)), paste("area", toString(m$area_seccion_m2)))
comprobar(cerca(m$vol_copa_m3, c(0.045, 0.04)), paste("volumen", toString(m$vol_copa_m3)))
comprobar(cerca(m$h_base_copa_m, c(0.5, 0.3)), paste("base", toString(m$h_base_copa_m)))
comprobar(cerca(P$MH_max, c(1.05, 0.75)), paste("MH", toString(P$MH_max)))
comprobar(cerca(m$prof_copa_m, c(0.55, 0.45)), paste("profundidad", toString(m$prof_copa_m)))
comprobar(cerca(m$anchura_max_m, c(1.0, 1.0)), paste("anchura máxima", toString(m$anchura_max_m)))
comprobar(cerca(m$h_anchura_max_m, c(0.95, 0.35)), paste("altura de la anchura máxima", toString(m$h_anchura_max_m)))
comprobar(cerca(m$TRV_m3_ha[1], 1968.75), paste("TRV", toString(m$TRV_m3_ha)))
comprobar(all(is.na(metricas_copa(P)$TRV_m3_ha)), "sin distancia entre hileras no debe haber TRV")
perf <- perfil_anchura(P)
comprobar(identical(as.integer(sort(perf$franja)), 0:7), paste("franjas", toString(perf$franja)))
f_ <- function(j) perf[perf$franja == j]
comprobar(f_(0)$n_ocupadas == 1 && cerca(f_(0)$anchura_media_m, 1.0), "franja 0")
comprobar(cerca(f_(0)$anchura_media_con_huecos_m, 0.5), "franja 0 con huecos")
comprobar(cerca(f_(2)$anchura_media_m, 0.6) && cerca(f_(2)$anchura_media_con_huecos_m, 0.3), "franja 2")
comprobar(cerca(f_(3)$anchura_media_m, 0.85) && f_(3)$frac_ocupada == 1.0, "franja 3")
comprobar(cerca(f_(0)$h_ini_m, 0.3) && cerca(f_(7)$h_fin_m, 1.1), "límites de franja")
por_tramo <- perfil_anchura(P, c(0, 1))
comprobar(setequal(unique(por_tramo$grupo), c(0, 1)), "grupos del perfil por tramo")
comprobar(max(por_tramo$franja[por_tramo$grupo == 1]) == 4, "el tramo 1 llega hasta la franja 4")
cat(sprintf("OK: área [%s], base [%s], TRV %.2f, perfil franja 3 = %.3f m\n", toString(m$area_seccion_m2), toString(m$h_base_copa_m),
            m$TRV_m3_ha[1], f_(3)$anchura_media_m))
