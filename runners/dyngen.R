suppressMessages({
    library(GillespieSSA2)
    library(dplyr)
    library(arrow)
    library(optparse)
})

options <- parse_args(OptionParser(option_list = list(
  make_option("--out", type = "character", help = "output path prefix"),
  make_option("--model", type = "character", help = "path to the .model.rds"),
  make_option("--samples", type = "integer", default = 10L),
  make_option("--trajectories", type = "integer"),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--method", type = "character", default = "exact"),
  make_option("--tau", type = "double", default = 300 / 3600),
  make_option("--mean-firings", type = "double", default = 100)
)))

Sys.setenv(R_MAKEVARS_USER = normalizePath("runners/rcpp-nobounds.mk"))

started <- proc.time()[["elapsed"]]
report <- function(stage) {
    cat(sprintf("[%.0fs] %s\n", proc.time()[["elapsed"]] - started, stage))
    flush.console()
}

report("loading model")
model <- readRDS(options$model)
simulation_system <- model$simulation_system

duration <- model$simulation_params$burn_time + model$simulation_params$total_time
interval <- duration / options$samples
trajectories <- if (is.null(options$trajectories)) {
    model$simulation_params$experiment_params$num_simulations[[1]]
} else {
    options$trajectories
}

method <- switch(options$method,
    exact = ssa_exact(),
    etl = ssa_etl(tau = options$tau),
    btl = ssa_btl(mean_firings = options$`mean-firings`),
    stop("unknown method: ", options$method))

report("building")
buffer_ids <- unique(unlist(lapply(simulation_system$reactions, function(r) r$buffer_ids)))
build_seconds <- system.time(
    reactions <- compile_reactions(
        simulation_system$reactions,
        state_ids = simulation_system$molecule_ids,
        params = simulation_system$parameters,
        buffer_ids = buffer_ids,
        hardcode_params = FALSE,
        fun_by = 1000L
    )
)[["elapsed"]]

report("simulating")
set.seed(options$seed)

simulate_seconds <- 0
runs <- lapply(seq_len(trajectories), function(i) {
    elapsed <- system.time(
        out <- ssa(
            initial_state = simulation_system$initial_state,
            reactions = reactions,
            final_time = duration,
            params = simulation_system$parameters,
            method = method,
            census_interval = interval,
            stop_on_neg_state = FALSE
        )
    )[["elapsed"]]
    simulate_seconds <<- simulate_seconds + elapsed

    counts <- as.matrix(out$state)
    tibble(
        path = as.character(i),
        sample = rep(round(out$time / interval), times = ncol(counts)),
        t = rep(out$time, times = ncol(counts)),
        name = rep(colnames(counts), each = nrow(counts)),
        value = as.integer(as.vector(counts))
    ) |> filter(sample > 0)
})

report("writing")
dir.create(options$out, recursive = TRUE, showWarnings = FALSE)

rss_peak_bytes <- function() {
  if (!file.exists("/proc/self/status")) return("")
  v <- grep("^VmHWM:", readLines("/proc/self/status"), value = TRUE)
  as.numeric(gsub("[^0-9]", "", v)) * 1024
}

writeLines(
    c("metric,value",
      paste0("build_seconds,", build_seconds),
      paste0("simulate_seconds,", simulate_seconds),
      paste0("rss_peak_bytes,", rss_peak_bytes())),
    file.path(options$out, "metrics.csv"))

writeLines(
    c("engine=dyngen",
      paste0("method=", options$method),
      paste0("tau=", options$tau),
      paste0("trajectories=", trajectories),
      paste0("model=", options$model),
      paste0("makevars=", Sys.getenv("R_MAKEVARS_USER"))),
    file.path(options$out, "config.txt"))

write_feather(bind_rows(runs), file.path(options$out, "counts.arrow"))

cat(sprintf("%s: %.1f s simulate, %.1f s build, %d count rows\n",
    options$out, simulate_seconds, build_seconds, sum(sapply(runs, nrow))))
