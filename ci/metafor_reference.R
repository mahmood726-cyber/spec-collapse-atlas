#!/usr/bin/env Rscript
# Compute metafor reference values (REML/DL tau^2 and pooled estimate) for each
# built-in dataset, so CI can validate the pure-Python engine against metafor.
# Reads ci/datasets.json (single source of truth, exported from datasets.py),
# writes ci/metafor_reference.json.
suppressMessages({
  library(jsonlite)
  library(metafor)
})

ds <- fromJSON("ci/datasets.json", simplifyVector = FALSE)
out <- list()
for (name in names(ds)) {
  yi <- unlist(ds[[name]]$yi)
  vi <- unlist(ds[[name]]$vi)
  reml <- rma(yi = yi, vi = vi, method = "REML")
  dl   <- rma(yi = yi, vi = vi, method = "DL")
  out[[name]] <- list(
    reml_tau2 = as.numeric(reml$tau2),
    reml_est  = as.numeric(reml$beta[1]),
    dl_tau2   = as.numeric(dl$tau2)
  )
  cat(sprintf("%-10s REML tau2=%.4f est=%.4f  DL tau2=%.4f\n",
              name, reml$tau2, reml$beta[1], dl$tau2))
}
write_json(out, "ci/metafor_reference.json", auto_unbox = TRUE, digits = 8)
cat("wrote ci/metafor_reference.json\n")
