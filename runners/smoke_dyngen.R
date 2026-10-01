suppressMessages({library(dyngen); library(tibble)})
set.seed(1)

backbone <- backbone(
  module_info = tribble(
    ~module_id, ~basal, ~burn, ~independence, ~color,
    "M1",       1,      TRUE,  1,             "#FF0000",
    "M2",       0,      FALSE, 1,             "#0000FF"),
  module_network = tribble(
    ~from, ~to,  ~effect, ~strength, ~hill,
    "M1",  "M2", 1L,      4,         2),
  expression_patterns = tribble(
    ~from,   ~to,  ~module_progression, ~start, ~burn, ~time,
    "sBurn", "sA", "+M1",               TRUE,   TRUE,  100,
    "sA",    "sB", "+M2",               FALSE,  FALSE, 100)
)

model <- initialise_model(
  backbone = backbone, num_cells = 1, num_tfs = 1, num_targets = 0, num_hks = 0,
  verbose = FALSE,
  simulation_params = simulation_default(
    census_interval = 10, ssa_algorithm = ssa_etl(tau = 300 / 3600),
    experiment_params = simulation_type_wild_type(num_simulations = 1)))

model <- generate_tf_network(model)
model <- generate_feature_network(model)
model <- generate_kinetics(model)
model <- generate_gold_standard(model)
model <- generate_cells(model)

fi <- model$feature_info
counts <- model$simulations$counts
cat(sprintf("dyngen %s: %d genes (%s), mean mRNA/gene = %.1f\n",
  as.character(packageVersion("dyngen")), nrow(fi), paste(fi$feature_id, collapse = ","), mean(Matrix::colMeans(counts[, fi$mol_mrna, drop = FALSE]))))
