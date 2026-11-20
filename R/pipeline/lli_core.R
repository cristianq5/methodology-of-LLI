# lli_core.R: núcleo numérico del LLI sin E/S (réplica de lli_core.py)
# MH, MW, envolvente de copa, clasificación de rayos (Bt / Bi), P, L y LLI.
# Marco de hilera: s (a lo largo), t (perpendicular, + izquierda), h = HAG.
#   LLI = MH * MW * (1 - P),  P = Bi / Bt  con  Bi = rayos que ATRAVIESAN la copa  (P = porosidad, L = 1 - P).
# Convención de matrices: celdas [rebanada, franja] = matriz n_reb x nf (índice lineal i + n_reb*j + 1, i,j desde 0).

# --------------------------------------------------------------------------------------
# MH y MW por rebanada, y envolvente por celda (rebanada x franja)
# --------------------------------------------------------------------------------------
perfiles_rebanada <- function(s, t, h, s_ini, n_reb, rebanada = 0.1, franja = 0.1, umbral = 0.3,
                              min_pts_franja = 5, min_pts_reb = 50, mh_percentil = 99.0, mw_percentil_t = c(2.5, 97.5),
                              n_franjas = NULL, env_percentil_t = c(1.0, 99.0)) {
  ok <- h >= umbral
  s <- s[ok]; t <- t[ok]; h <- h[ok]
  i <- floor((s - s_ini) / rebanada)
  j <- floor((h - umbral) / franja + 1e-9)
  v <- i >= 0 & i < n_reb & j >= 0
  i <- i[v]; j <- j[v]; t <- t[v]; h <- h[v]
  nf <- if (is.null(n_franjas)) (if (length(j)) max(j) + 1L else 1L) else n_franjas
  j <- pmin(j, nf - 1L)
  P <- list(n_reb = n_reb, nf = nf, rebanada = rebanada, franja = franja, umbral = umbral, s_ini = s_ini)
  P$n_pts <- tabulate(i + 1L, n_reb)
  tam <- n_reb * nf
  mk <- function() rep(NA_real_, tam)
  tmin <- mk(); tmax <- mk(); tlo <- mk(); thi <- mk(); elo <- mk(); ehi <- mk()
  n_cel <- tabulate(i + j * n_reb + 1L, tam)
  if (length(i)) {
    dt <- data.table(cel = i + j * n_reb + 1L, t = t)
    pr <- c(mw_percentil_t, env_percentil_t) / 100
    g <- dt[, {
      q <- if (.N >= min_pts_franja) quantile(t, pr, type = 7, names = FALSE) else rep(NA_real_, 4)
      .(tmn = min(t), tmx = max(t), q1 = q[1], q2 = q[2], q3 = q[3], q4 = q[4])
    }, by = cel]
    tmin[g$cel] <- g$tmn; tmax[g$cel] <- g$tmx; tlo[g$cel] <- g$q1; thi[g$cel] <- g$q2; elo[g$cel] <- g$q3; ehi[g$cel] <- g$q4
  }
  mm <- function(x) matrix(x, n_reb, nf)
  P$n_cel <- mm(n_cel); P$tmin <- mm(tmin); P$tmax <- mm(tmax); P$tlo <- mm(tlo); P$thi <- mm(thi); P$elo <- mm(elo); P$ehi <- mm(ehi)
  P$ocupada <- P$n_cel >= min_pts_franja
  ancho_e <- ifelse(P$ocupada, P$tmax - P$tmin, NA_real_)
  ancho_p <- ifelse(P$ocupada, P$thi - P$tlo, NA_real_)
  fm <- function(a) { r <- rowMeans(a, na.rm = TRUE); r[is.nan(r)] <- NA_real_; r }
  P$MW_max <- fm(ancho_e); P$MW_p <- fm(ancho_p)
  mh_max <- rep(NA_real_, n_reb); mh_p <- rep(NA_real_, n_reb)
  if (length(i)) {
    gh <- data.table(i = i, h = h)[, .(mx = max(h), mp = quantile(h, mh_percentil / 100, type = 7, names = FALSE)), by = i]
    mh_max[gh$i + 1L] <- gh$mx; mh_p[gh$i + 1L] <- gh$mp
  }
  P$MH_max <- mh_max; P$MH_p <- mh_p
  malas <- P$n_pts < min_pts_reb
  P$MH_max[malas] <- NA; P$MH_p[malas] <- NA; P$MW_max[malas] <- NA; P$MW_p[malas] <- NA
  P
}

# --------------------------------------------------------------------------------------
# Métricas de copa (réplica de anchuras_franja, metricas_copa y perfil_anchura de lli_core.py)
# --------------------------------------------------------------------------------------
# Anchura de cada celda (matriz n_reb x nf): estricta (tmax - tmin) y robusta (percentiles de t).
# NA si la franja no está ocupada (< min_pts_franja puntos) o si la rebanada no es válida (< min_pts_reb puntos).
anchuras_franja <- function(P) {
  valida <- matrix(is.finite(P$MH_max), P$n_reb, P$nf)          # misma validez para todas las franjas de la rebanada
  m <- P$ocupada & valida
  list(ae = ifelse(m, P$tmax - P$tmin, NA_real_), ap = ifelse(m, P$thi - P$tlo, NA_real_))
}

# Métricas de copa por rebanada, a partir de la misma malla rebanada x franja que MH y MW (sin leer otra vez la nube).
#   area_seccion_m2 : área de la sección = suma de (anchura de franja x franja) de las franjas ocupadas (estricta; _p robusta).
#   vol_copa_m3     : area_seccion_m2 x rebanada.
#   h_base_copa_m   : borde inferior de la franja ocupada más baja (si la copa baja de 'umbral', vale 'umbral').
#   prof_copa_m     : MH_max - h_base_copa_m.
#   anchura_max_m   : anchura de la franja más ancha; h_anchura_max_m = altura (centro) de esa franja.
#   TRV_m3_ha       : Tree Row Volume = MH x MW x 10 000 / distancia entre hileras (estricta; TRV_p con MH_p, MW_p).
# Devuelve una lista de vectores (n_reb).
metricas_copa <- function(P, dist_hileras = NA_real_) {
  a <- anchuras_franja(P); ae <- a$ae; ap <- a$ap
  valida <- is.finite(P$MH_max)
  occ <- matrix(is.finite(ae), P$n_reb, P$nf)
  area_e <- ifelse(valida, rowSums(ifelse(occ, ae, 0)) * P$franja, NA_real_)
  area_p <- ifelse(valida, rowSums(ifelse(is.finite(ap), ap, 0)) * P$franja, NA_real_)
  hay <- rowSums(occ) > 0 & valida
  # índice (desde 0) de la primera franja ocupada y de la franja más ancha (primera en caso de empate, como numpy.argmax)
  j_min <- ifelse(hay, apply(occ, 1, function(r) if (any(r)) which(r)[1] else 1L) - 1L, 0L)
  ae_m <- matrix(ifelse(occ, ae, -Inf), P$n_reb, P$nf)
  j_anch <- ifelse(hay, apply(ae_m, 1, which.max) - 1L, 0L)
  h_base <- ifelse(hay, P$umbral + j_min * P$franja, NA_real_)
  anch_max <- ifelse(hay, ae[cbind(seq_len(P$n_reb), j_anch + 1L)], NA_real_)
  h_anch <- ifelse(hay, P$umbral + (j_anch + 0.5) * P$franja, NA_real_)
  dist_hileras <- as.numeric(dist_hileras)
  trv <- P$MH_max * P$MW_max * 1e4 / dist_hileras
  trv_p <- P$MH_p * P$MW_p * 1e4 / dist_hileras
  trv[!is.finite(trv)] <- NA_real_; trv_p[!is.finite(trv_p)] <- NA_real_
  list(area_seccion_m2 = area_e, area_seccion_p_m2 = area_p,
       vol_copa_m3 = area_e * P$rebanada, vol_copa_p_m3 = area_p * P$rebanada,
       h_base_copa_m = h_base, prof_copa_m = P$MH_max - h_base,
       anchura_max_m = anch_max, h_anchura_max_m = h_anch,
       TRV_m3_ha = trv, TRV_p_m3_ha = trv_p)
}

# Perfil de anchura por altura: para cada franja y cada grupo de rebanadas (NULL = toda la hilera; si no, un vector (n_reb) con
# el número de grupo de cada rebanada, p. ej. el tramo de 1,5 m), sobre las rebanadas válidas:
#   n_rebanadas, n_ocupadas, frac_ocupada,
#   anchura_media_m            : media de la anchura de la franja donde está ocupada (estricta),
#   anchura_media_p_m          : idem con la anchura robusta,
#   anchura_media_con_huecos_m : media contando 0 donde la franja está vacía (= anchura_media_m x frac_ocupada).
# Solo hasta la franja ocupada más alta del grupo. Devuelve un data.table (una fila por grupo y franja; columna 'grupo').
perfil_anchura <- function(P, grupos = NULL) {
  a <- anchuras_franja(P); ae <- a$ae; ap <- a$ap
  valida <- is.finite(P$MH_max)
  g <- if (is.null(grupos)) rep(0L, P$n_reb) else as.vector(grupos)
  filas <- list()
  for (k in sort(unique(g[valida]))) {
    m <- valida & g == k
    e <- ae[m, , drop = FALSE]; p_ <- ap[m, , drop = FALSE]
    occ <- is.finite(e)
    if (!any(occ)) next
    j_top <- max(which(colSums(occ) > 0)) - 1L
    n_m <- sum(m)
    for (j in 0:j_top) {
      n_o <- sum(occ[, j + 1L])
      pj <- p_[, j + 1L]
      filas[[length(filas) + 1L]] <- list(grupo = k, franja = j,
        h_ini_m = round(P$umbral + j * P$franja, 3), h_fin_m = round(P$umbral + (j + 1) * P$franja, 3),
        n_rebanadas = n_m, n_ocupadas = n_o, frac_ocupada = n_o / n_m,
        anchura_media_m = if (n_o) mean(e[, j + 1L], na.rm = TRUE) else NA_real_,
        anchura_media_p_m = if (any(is.finite(pj))) mean(pj, na.rm = TRUE) else NA_real_,
        anchura_media_con_huecos_m = sum(e[, j + 1L], na.rm = TRUE) / n_m)
    }
  }
  rbindlist(filas)
}

# (tmin, tmax, ocupada) por celda según config (estricta o percentil)
envolvente <- function(P, estadistico = "max") {
  if (estadistico == "max") return(list(tmin = P$tmin, tmax = P$tmax, occ = P$ocupada))
  list(tmin = P$elo, tmax = P$ehi, occ = P$ocupada & is.finite(P$elo))
}

# --------------------------------------------------------------------------------------
# Clasificación de rayos
# --------------------------------------------------------------------------------------
# O (N x 3): origen (s, t, z); D (N x 3): dirección unitaria; R (N): distancia al eco (Inf = sin eco)
# dtm_fn(s, t) -> z del suelo. Devuelve lista de vectores por rayo: valido, entra, atraviesa, lado, i, j.
clasificar_rayos <- function(O, D, R, env, P, dtm_fn, paso = 0.025, az_max_deg = 25.0, el_rango_deg = c(-60, 15),
                             origen_t_abs = c(0.8, 3.0), alcance_max = 30.0, tol_extremo = NULL, sub_lote = 15000L) {
  tmin <- env$tmin; tmax <- env$tmax; occ <- env$occ
  n_reb <- nrow(occ); nf <- ncol(occ)
  N <- length(R)
  out <- list(valido = logical(N), entra = logical(N), atraviesa = logical(N), lado = integer(N), i = rep(-1L, N), j = rep(-1L, N))
  if (N == 0) return(out)
  side <- as.integer(sign(O[, 2]))
  ta <- abs(O[, 2])
  el <- asin(pmin(pmax(D[, 3], -1), 1)) * 180 / pi
  hacia <- (D[, 2] * side) < -1e-6
  az <- atan2(abs(D[, 1]), -D[, 2] * side) * 180 / pi
  v <- ta >= origen_t_abs[1] & ta <= origen_t_abs[2] & hacia & az <= az_max_deg & el >= el_rango_deg[1] & el <= el_rango_deg[2]
  v[is.na(v)] <- FALSE
  out$valido <- v; out$lado <- side
  if (!any(v)) return(out)
  tmin_f <- ifelse(occ, tmin, Inf); tmax_f <- ifelse(occ, tmax, -Inf)
  tm <- if (any(occ)) max(abs(c(tmin[occ], tmax[occ]))) else 0
  if (!is.finite(tm) || tm == 0) return(out)
  K <- as.integer(ceiling(2 * tm / paso)) + 1L
  tol <- if (is.null(tol_extremo)) paso else tol_extremo
  s_ini <- P$s_ini; reb <- P$rebanada; fr <- P$franja; um <- P$umbral
  dentro <- function(sk_, tk_, hk_, infl) {
    ii <- floor((sk_ - s_ini) / reb)
    jj <- floor((hk_ - um) / fr + 1e-9)
    okc <- !is.na(ii) & !is.na(jj) & ii >= 0 & ii < n_reb & jj >= 0 & jj < nf
    ii2 <- ifelse(okc, ii, 0); jj2 <- ifelse(okc, jj, 0)
    li <- ii2 + n_reb * jj2 + 1
    ins <- okc & occ[li] & tk_ >= tmin_f[li] - infl & tk_ <= tmax_f[li] + infl
    list(ins = ins, ii = ii2, jj = jj2)
  }
  iv_all <- which(v)
  kk <- 0:(K - 1L)
  for (a in seq(1L, length(iv_all), by = sub_lote)) {
    iv <- iv_all[a:min(a + sub_lote - 1L, length(iv_all))]
    n <- length(iv)
    Rv <- pmin(R[iv], alcance_max)
    Ov <- O[iv, , drop = FALSE]; Dv <- D[iv, , drop = FALSE]; sv <- side[iv]
    tk <- outer(sv, tm - paso * kk)                    # de lado propio a lado opuesto (n x K)
    lam <- (tk - Ov[, 2]) / Dv[, 2]
    sk <- Ov[, 1] + lam * Dv[, 1]
    zk <- Ov[, 3] + lam * Dv[, 3]
    hk <- zk - matrix(dtm_fn(as.vector(sk), as.vector(tk)), n, K)
    d0 <- dentro(sk, tk, hk, 0.0)
    ins <- d0$ins & lam > 0 & lam <= Rv
    ins[is.na(ins)] <- FALSE
    fin_t <- Ov[, 2] + Rv * Dv[, 2]; fin_s <- Ov[, 1] + Rv * Dv[, 1]
    fin_h <- Ov[, 3] + Rv * Dv[, 3] - dtm_fn(fin_s, fin_t)
    hay_eco <- is.finite(R[iv]) & R[iv] <= alcance_max
    df <- dentro(fin_s, fin_t, fin_h, tol)
    fins <- df$ins & hay_eco
    alguno <- rowSums(ins) > 0
    primero <- max.col(ins * 1, ties.method = "first")
    entra <- alguno | fins
    li <- cbind(seq_len(n), primero)
    i_ent <- ifelse(alguno, d0$ii[li], df$ii)
    j_ent <- ifelse(alguno, d0$jj[li], df$jj)
    atr <- entra & !fins
    out$entra[iv] <- entra; out$atraviesa[iv] <- atr
    out$i[iv] <- ifelse(entra, i_ent, -1L); out$j[iv] <- ifelse(entra, j_ent, -1L)
  }
  out
}

# Bt y Bi por (lado, rebanada, franja) -> arrays dim c(2, n_reb, nf); lado 1 = izquierdo (+), 2 = derecho (-)
acumular_conteos <- function(cl, n_reb, nf, conteos = NULL) {
  if (is.null(conteos)) conteos <- list(Bt = array(0, c(2, n_reb, nf)), Bi = array(0, c(2, n_reb, nf)))
  m <- cl$entra
  if (any(m)) {
    l <- as.integer(cl$lado[m] < 0)
    lin <- l + 2L * (cl$i[m] + n_reb * cl$j[m]) + 1L
    n <- 2L * n_reb * nf
    conteos$Bt <- conteos$Bt + array(tabulate(lin, n), c(2, n_reb, nf))
    a <- cl$atraviesa[m]
    if (any(a)) conteos$Bi <- conteos$Bi + array(tabulate(lin[a], n), c(2, n_reb, nf))
  }
  conteos
}

# P por lado y rebanada: solo franjas con Bt >= min_rayos_celda. modo "suma_rayos" | "media_franjas"
# Devuelve Pl (2 x n_reb), nfi (2 x n_reb), Pm (n_reb; NA si falta un lado)
p_por_rebanada <- function(Bt, Bi, min_rayos_celda = 30, min_franjas = 3, modo = "suma_rayos", j_min = 0) {
  fi <- Bt >= min_rayos_celda
  if (j_min > 0) fi[, , seq_len(j_min)] <- FALSE
  sm <- function(a) apply(a, c(1, 2), sum)
  nfi <- sm(fi * 1)
  if (modo == "suma_rayos") {
    bt <- sm(ifelse(fi, Bt, 0)); bi <- sm(ifelse(fi, Bi, 0))
    Pl <- ifelse(bt > 0, bi / pmax(bt, 1), NA_real_)
  } else {
    r <- ifelse(fi, Bi / pmax(Bt, 1), 0)
    Pl <- ifelse(nfi > 0, sm(r) / pmax(nfi, 1), NA_real_)
  }
  Pl <- ifelse(nfi >= min_franjas, Pl, NA_real_)
  Pm <- ifelse(is.finite(Pl[1, ]) & is.finite(Pl[2, ]), (Pl[1, ] + Pl[2, ]) / 2, NA_real_)
  list(Pl = Pl, nfi = nfi, Pm = Pm)
}

# Cobertura vertical de P (QA-05): h_max (2 x n_reb), frac (2 x n_reb), frac_min (n_reb)
cobertura_P <- function(Bt, MH, umbral, franja, min_rayos_celda = 30, j_min = 0) {
  fi <- Bt >= min_rayos_celda
  if (j_min > 0) fi[, , seq_len(j_min)] <- FALSE
  nf <- dim(Bt)[3]
  jj <- array(rep(seq_len(nf) - 1L, each = 2L * dim(Bt)[2]), dim(Bt))
  jmax <- apply(ifelse(fi, jj, -1L), c(1, 2), max)
  h_max <- ifelse(jmax >= 0, umbral + (jmax + 1) * franja, NA_real_)
  nfi <- apply(fi * 1, c(1, 2), sum)
  frac <- pmin(pmax(sweep(nfi * franja, 2, (MH - umbral), "/"), 0), 1)
  frac[!is.finite(frac)] <- NA_real_
  list(h_max = h_max, frac = frac, frac_min = pmin(frac[1, ], frac[2, ]))
}

lli <- function(MH, MW, P) MH * MW * (1.0 - P)

# Máscara de valores dentro de media +- n_de desviaciones típicas (ignora NA; DE poblacional como numpy.nanstd)
filtro_atipicos <- function(x, n_de = 2.5) {
  x <- as.numeric(x)
  ok <- is.finite(x)
  m <- mean(x[ok]); d <- sqrt(mean((x[ok] - m)^2))
  if (!is.finite(d) || d == 0) return(ok)
  ok & abs(x - m) <= n_de * d
}
