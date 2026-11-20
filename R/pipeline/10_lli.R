# 10_lli.R: LLI = MH * MW * (1 - P) por rebanada; filtro +-2,5 DE por hilera; media por hilera y por tramos de 1,5 m
# (réplica de 10_lli.py). Se calcula con MH/MW estrictos y percentil, y con P de cada vía. QA de signo de L.
if (!exists("DIR_SCRIPTS")) { a <- grep("^--file=", commandArgs(FALSE), value = TRUE); DIR_SCRIPTS <- dirname(normalizePath(sub("^--file=", "", a[1]))) }
source(file.path(DIR_SCRIPTS, "utils.R")); source(file.path(DIR_SCRIPTS, "lli_core.R"))

# Métricas de copa de 09 (lli_core.R::metricas_copa): se agregan sin filtro de atípicos (las rebanadas sin copa ya son NA).
COPA_MEDIA <- c("area_seccion_m2", "area_seccion_p_m2", "h_base_copa_m", "prof_copa_m", "anchura_max_m", "h_anchura_max_m",
                "TRV_m3_ha", "TRV_p_m3_ha")
COPA_SUMA <- c("vol_copa_m3", "vol_copa_p_m3")

# Medias (_media) y volúmenes totales (_total, m³, y _por_m, m³ por metro de hilera) de un grupo de rebanadas (réplica de agregar_copa).
agregar_copa <- function(g, cfg) {
  f <- list()
  for (c_ in intersect(COPA_MEDIA, names(g))) f[[paste0(c_, "_media")]] <- nanmean(g[[c_]])
  largo <- nrow(g) * cfg$metricas$rebanada_m
  for (c_ in intersect(COPA_SUMA, names(g))) {
    tot <- if (all(is.na(g[[c_]]))) NA_real_ else sum(g[[c_]], na.rm = TRUE)      # como pandas sum(min_count = 1)
    f[[paste0(c_, "_total")]] <- tot
    f[[paste0(c_, "_por_m")]] <- if (largo > 0) tot / largo else NA_real_
  }
  if ("dist_hileras_m" %in% names(g)) f$dist_hileras_m <- g$dist_hileras_m[1]
  f
}

main_10 <- function(cfg = cargar_config()) {
  out <- dir_salida(cfg, "10_lli")
  c <- cfg$lli
  ent <- file.path(dir_salida(cfg, "09_metricas"), "metricas_rebanadas.csv")
  con_paso(cfg, "10_lli", entrada = ent, parametros = c, fn = function(p) {
    df <- fread(ent, colClasses = list(character = "id_fila"))
    df[, id_fila := as.character(id_fila)]
    atr <- intersect(atributos_croquis(cfg), names(df))     # atributos del archivo vectorial (variedad, bloque...)
    p$n_antes <- nrow(df)
    variantes <- list(estricta = c("MH_max", "MW_max"), percentil = c("MH_p", "MW_p"))
    fuentes <- c("P", intersect(c("P_h06", "P_mejor"), names(df)), intersect(c("P_pcap", "P_rayos", "P_pcap_todas"), names(df)))
    cols <- character(0)
    for (vn in names(variantes)) for (pf in fuentes) {
      nom <- if (pf == "P") paste0("LLI_", vn) else paste0("LLI_", vn, "_", substring(pf, 3))
      df[[nom]] <- lli(df[[variantes[[vn]][1]]], df[[variantes[[vn]][2]]], df[[pf]])
      cols <- c(cols, nom)
    }
    df[, L := 1.0 - P]
    # filtro de atípicos por hilera
    for (col in cols) df[[paste0(col, "_ok")]] <- unlist(lapply(split(df[[col]], factor(df$id_fila, levels = unique(df$id_fila))), filtro_atipicos, n_de = c$filtro_atipicos_de), use.names = FALSE)
    df[, tramo := as.integer(floor((s_ini - min(s_ini) + 1e-9) / cfg$hilera$tramo_m)), by = id_fila]
    agg_h <- list(); agg_t <- list()
    for (k in unique(df$id_fila)) {
      g <- df[id_fila == k]
      n_p <- sum(!is.na(g$P))
      mot <- g[is.na(P), P_motivo_nan]; mot[is.na(mot)] <- ""
      vc <- if (length(mot)) { tb <- table(factor(mot, levels = unique(mot))); tb[order(-tb)] } else integer(0)
      estado <- if (n_p == nrow(g)) "completa" else if (n_p == 0) paste0("SIN P: ", if (length(vc)) names(vc)[1] else "motivo no registrado") else
        sprintf("parcial: %d/%d rebanadas con P; sin P por: %s", n_p, nrow(g), paste(sprintf("%s (%d)", names(vc), as.integer(vc)), collapse = "; "))
      vcm <- table(g$via_P_mejor[nzchar(g$via_P_mejor)]); vcm <- vcm[order(-vcm)]
      fila <- c(list(id_fila = k), lapply(stats::setNames(atr, atr), function(a_) g[[a_]][1]), list(n_rebanadas = nrow(g), n_rebanadas_con_P = n_p, estado_P = estado,
                   n_ventanas_validas = g$n_ventanas_validas[1], n_rebanadas_P_mejor = sum(!is.na(g$P_mejor)),
                   origen_P_mejor = paste(sprintf("%s=%d", names(vcm), as.integer(vcm)), collapse = "; "),
                   P_media = nanmean(g$P), P_h06_media = nanmean(g$P_h06),
                   h_max_P_izq_media = nanmean(g$h_max_P_izq), h_max_P_der_media = nanmean(g$h_max_P_der),
                   h_max_P_izq_mediana = nanmedian(g$h_max_P_izq), h_max_P_der_mediana = nanmedian(g$h_max_P_der),
                   cob_MH_media = nanmean(g$cob_MH), MH_max_media = nanmean(g$MH_max), MW_max_media = nanmean(g$MW_max),
                   MH_p_media = nanmean(g$MH_p), MW_p_media = nanmean(g$MW_p)), agregar_copa(g, cfg))
      for (col in cols) {
        v <- g[[col]][g[[paste0(col, "_ok")]]]
        fila[[paste0(col, "_media")]] <- nanmean(v)
        fila[[paste0(col, "_n")]] <- length(v)
        fila[[paste0(col, "_descartadas")]] <- sum(!is.na(g[[col]])) - length(v)
      }
      agg_h[[length(agg_h) + 1]] <- fila
      for (tr in sort(unique(g$tramo))) {
        gt <- g[tramo == tr]
        if (nrow(gt) < c$min_rebanadas_tramo) next
        f2 <- c(list(id_fila = k), lapply(stats::setNames(atr, atr), function(a_) g[[a_]][1]), list(tramo = tr, s_ini = min(gt$s_ini), s_fin = max(gt$s_ini) + cfg$metricas$rebanada_m, n_rebanadas = nrow(gt),
                   P_media = nanmean(gt$P), P_h06_media = nanmean(gt$P_h06), h_max_P_izq_media = nanmean(gt$h_max_P_izq),
                   h_max_P_der_media = nanmean(gt$h_max_P_der), cob_MH_media = nanmean(gt$cob_MH)), agregar_copa(gt, cfg))
        for (col in cols) f2[[col]] <- nanmean(gt[[col]][gt[[paste0(col, "_ok")]]])
        agg_t[[length(agg_t) + 1]] <- f2
      }
    }
    H <- rbindlist(agg_h, fill = TRUE); T <- rbindlist(agg_t, fill = TRUE)
    fwrite(df, file.path(out, "lli_rebanadas.csv")); fwrite(H, file.path(out, "lli_hileras.csv")); fwrite(T, file.path(out, "lli_tramos.csv"))
    res <- list(hileras = H)
    rng <- list()
    for (k in unique(T$id_fila)) { g <- T[id_fila == k]; rng[[k]] <- lapply(intersect(c("LLI_estricta", "LLI_percentil", "LLI_estricta_h06"), names(g)), function(cc) range(g[[cc]], na.rm = TRUE)) }
    res$rango_tramos_1_5m <- rng
    lo <- c$rango_esperado_m2[1]; hi <- c$rango_esperado_m2[2]
    res$avisos_rango <- Filter(Negate(is.null), lapply(seq_len(nrow(H)), function(i) {
      v <- H$LLI_estricta_media[i]
      if (is.na(v) || !(v >= lo && v <= hi)) sprintf("%s: LLI_estricta %s fuera de [%g,%g]", H$id_fila[i], format(v, digits = 3), lo, hi) else NULL }))
    if (isTRUE(c$qa_signo_L)) {
      dens <- df$n_pts / (df$MH_max * df$MW_max)
      cor_ <- list()
      for (nm in c("densidad_puntos", "MH", "MW")) {
        x <- switch(nm, densidad_puntos = dens, MH = df$MH_max, MW = df$MW_max)
        ok <- is.finite(x) & is.finite(df$L)
        cor_[[nm]] <- if (sum(ok) > 5) cor(x[ok], df$L[ok]) else NA
      }
      vals <- unlist(cor_); vals <- vals[!is.na(vals)]
      res$qa_signo_L <- list(corr_L = cor_, ok = length(vals) > 0 && mean(vals) > 0)
    }
    guardar_json(res, file.path(out, "lli_resumen.json"))
    p$info("LLI por hilera:\n%s", paste(capture.output(print(as.data.frame(H[, intersect(c("id_fila", atr, "P_media", "P_h06_media", "LLI_estricta_media", "LLI_percentil_media", "LLI_estricta_mejor_media", "cob_MH_media"), names(H)), with = FALSE]), digits = 4)), collapse = "\n"))
    p$info("QA signo L: %s | avisos: %s", toJSON(res$qa_signo_L, auto_unbox = TRUE, digits = 3), paste(res$avisos_rango, collapse = "; "))
    for (k in unique(df$id_fila)) {
      g <- df[id_fila == k]
      png_lineas(list(list(g$s_centro, g$LLI_estricta, "LLI estricta"), list(g$s_centro, g$LLI_percentil, "LLI percentil")),
                 file.path(out, sprintf("lli_hilera%s.png", k)), sprintf("Hilera %s: LLI por rebanada", k), "s (m)", "m2")
    }
    p$n_despues <- nrow(df)
    p$salida <- out
    list(H = H, T = T)
  })
}
if (sys.nframe() == 0) main_10()
