suppressMessages({
    library(dyngen)
    library(dplyr)
    library(jsonlite)
    library(optparse)
})

POOLS <- list(polymerases = 5e5, ribosomes = 2e6, proteasomes = 1e6)
ELONGATION <- 3.6

options <- parse_args(OptionParser(option_list = list(
  make_option("--out", type = "character", help = "output directory"),
  make_option("--genes", type = "integer", default = 150L),
  make_option("--trajectories", type = "integer", default = 100L),
  make_option("--seed", type = "integer", default = 1L),
  make_option("--samples", type = "integer", default = 10L),
  make_option("--backbone", type = "character", default = "bifurcating_cycle"),
  make_option("--species", type = "character", default = "premrnas,mrnas,proteins"),
  make_option("--schedule", type = "character", default = "schedule.json"),
  make_option("--promoter", type = "character", default = "equilibrium"),
  make_option("--reuse-model", action = "store_true", default = FALSE)
)))

species <- strsplit(options$species, ",")[[1]]

set.seed(options$seed)

backbone <- get(paste0("backbone_", options$backbone))()
num_tfs <- nrow(backbone$module_info)
num_targets <- round((options$genes - num_tfs) / 2)
num_hks <- options$genes - num_tfs - num_targets

reused <- file.path(options$out, "model.rds")

model <- if (options$`reuse-model` && file.exists(reused)) readRDS(reused) else
initialise_model(
    backbone = backbone,
    num_tfs = num_tfs,
    num_targets = num_targets,
    num_hks = num_hks,
    num_cells = options$trajectories * options$samples,
    simulation_params = simulation_default(
        experiment_params = simulation_type_wild_type(
            num_simulations = options$trajectories
        )
    ),
    verbose = FALSE
) |>
generate_tf_network() |>
generate_feature_network() |>
generate_kinetics()

duration <- model$simulation_params$burn_time + model$simulation_params$total_time

dir.create(options$out, recursive = TRUE, showWarnings = FALSE)
if (!options$`reuse-model` || !file.exists(reused)) saveRDS(model, reused)


feature_info <- model$feature_info
feature_network <- model$feature_network |>
    left_join(select(feature_info, to = feature_id, basal), by = "to") |>
    mutate(
        clamped_basal = pmin(pmax(basal, 1e-6), 1 - 1e-6),
        at = if (options$promoter == "thermodynamic")
            dissociation / strength^(1 / hill)
        else dissociation * ifelse(
            effect > 0,
            clamped_basal / strength,
            (1 - clamped_basal) / strength
        )^(1 / hill),
        k = -hill
    )

slots <- function(rows) lapply(seq_len(nrow(rows)), function(i)
    list(from = unbox(rows$from[i]), at = unbox(rows$at[i]), k = unbox(rows$k[i])))

regulation <- function(rows) if (options$promoter == "thermodynamic")
    list(slots = slots(rows), aggregate = unbox("thermodynamic")) else slots(rows)

genes <- lapply(seq_len(nrow(feature_info)), function(i) {
    gene <- feature_info[i, ]
    basal <- if (options$promoter == "thermodynamic") max(gene$basal, 1e-12)
        else min(max(gene$basal, 1e-6), 1 - 1e-6)
    regulators <- filter(feature_network, to == gene$feature_id)
    out <- list(
        name = unbox(gene$feature_id),
        base_rates = lapply(list(
            activation    = 1,
            deactivation  = (1 - basal) / basal,
            initiation    = gene$transcription_rate / POOLS$polymerases,
            transcription = ELONGATION,
            processing    = gene$splicing_rate,
            translation   = gene$translation_rate / POOLS$ribosomes,
            abortion      = 0,
            premrna_decay = gene$mrna_decay_rate,
            mrna_decay    = gene$mrna_decay_rate,
            protein_decay = gene$protein_decay_rate / POOLS$proteasomes
        ), unbox),
        activation = regulation(filter(regulators, effect > 0)),
        repression = regulation(filter(regulators, effect < 0)))
    out$color <- unbox(
        if (!is.na(gene$color)) gene$color
        else if (gene$is_tf) "#000000"
        else if (gene$is_hk) "#D3D3D3"
        else "#A9A9A9"
    )
    out
})

schedule <- list(
    method = unbox("RSSACR"),
    num_reps = unbox(options$trajectories),
    duration = unbox(duration),
    sample_interval = unbox(duration / options$samples),
    step = list(
        do = list(`{regulation/v1}` = list(
            method = list(`$` = unbox("method")),
            species = species,
            genes = genes
        )),
        step = list(
            list(`{add}` = lapply(POOLS, unbox)),
            list(
                branch = unbox(TRUE),
                step = list(
                    each = list(start = unbox(1), stop = list(`$` = unbox("num_reps"))),
                    step = list(
                        to = list(`$` = unbox("duration")),
                        step = list(`$` = unbox("sample_interval"))
                    )
                )
            )
        )
    )
)

write_json(schedule, file.path(options$out, options$schedule), digits = 17, auto_unbox = FALSE, pretty = TRUE)

cat(sprintf("%s/%s: %d genes, %d edges, %d trajectories, to = %g, %d samples, species %s\n",
    options$out, options$schedule, nrow(feature_info), nrow(feature_network),
    options$trajectories, duration, options$samples, options$species))
