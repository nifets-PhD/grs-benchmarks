#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(GillespieSSA2); library(arrow)})
Sys.setenv(R_MAKEVARS_USER = normalizePath("runners/rcpp-nobounds.mk"))

args <- commandArgs(trailingOnly = TRUE)
out <- args[1]; model <- args[2]
getarg <- function(flag, default) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) default else args[i + 1]
}
method   <- getarg("--method", "etl")
tau      <- as.numeric(getarg("--tau", "1"))
duration <- as.numeric(getarg("--duration", "20000"))
samples  <- as.integer(getarg("--samples", "10"))
ntraj    <- as.integer(getarg("--trajectories", "3"))
seed     <- as.integer(getarg("--seed", "1"))
stopneg  <- as.logical(getarg("--stop-on-neg", "FALSE"))

source(model)

t_build <- system.time({
  compiled <- compile_reactions(reactions, state_ids = state_ids,
                                params = c(dummy = 0), buffer_ids = NULL)
})[["elapsed"]]

ssa_method <- switch(method,
  etl   = ssa_etl(tau = tau),
  exact = ssa_exact(),
  btl   = ssa_btl(mean_firings = as.numeric(getarg("--mean-firings", "10"))),
  stop("unknown method: ", method))

frames <- list(); aborted <- 0
t_sim <- system.time({
  for (r in seq_len(ntraj)) {
    set.seed(seed + r - 1L)
    res <- ssa(initial_state = initial_state, reactions = compiled,
               final_time = duration, params = c(dummy = 0), method = ssa_method,
               census_interval = duration / samples,
               stop_on_neg_state = stopneg, verbose = FALSE,
               sim_name = sprintf("traj%d", r))
    if (nrow(res$state) < samples + 1L) aborted <- aborted + 1
    frames[[r]] <- list(t = res$time, state = res$state, path = r)
    cat(sprintf("[%s] trajectory %d/%d\n", format(Sys.time(), "%H:%M:%S"), r, ntraj),
        file = stderr())
  }
})[["elapsed"]]

dir.create(out, recursive = TRUE, showWarnings = FALSE)
rows <- do.call(rbind, lapply(frames, function(f) {
  st <- as.data.frame(f$state)
  colnames(st) <- state_ids
  st$t <- f$t
  long <- reshape(st, direction = "long", varying = state_ids,
                  v.names = "value", timevar = "name", times = state_ids,
                  idvar = "t")
  long$name <- state_names[match(long$name, state_ids)]
  long$path <- as.character(f$path)
  long$sample <- as.integer(round(long$t / (duration / samples)))
  long[, c("path", "sample", "t", "name", "value")]
}))
write_feather(rows, file.path(out, "counts.arrow"))

rss_peak_bytes <- function() {
  if (!file.exists("/proc/self/status")) return("")
  v <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", v)) * 1024
}

writeLines(c("engine=gillespiessa2", paste0("method=", method),
             paste0("tau=", tau), paste0("gillespiessa2_stop_on_neg_state=", stopneg),
             paste0("trajectories=", ntraj), paste0("seed=", seed),
             paste0("model=", model), paste0("makevars=", Sys.getenv("R_MAKEVARS_USER"))),
           file.path(out, "config.txt"))

writeLines(c("metric,value",
             paste0("rss_peak_bytes,", rss_peak_bytes()),
             paste0("build_seconds,", t_build),
             paste0("simulate_seconds,", t_sim),
             paste0("simulate_seconds_per_trajectory,", t_sim / ntraj),
             paste0("aborted_trajectories,", aborted)),
           file.path(out, "metrics.csv"))

cat(sprintf("gillespie %s: %.1f s simulate (%d traj, %.2f s/traj), %.1f s build, %d aborted, %d count rows\n",
            out, t_sim, ntraj, t_sim / ntraj, t_build, aborted, nrow(rows)))
