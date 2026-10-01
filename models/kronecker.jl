using ArgParse
using JSON

const MODELS = Dict(
    "kronecker" => (;
        base_rates = Dict{String, Float64}(
            "activation" => 3.0e-4,
            "deactivation" => 2.0e-3,
            "initiation" => 1.93e-7,
            "abortion" => 0.01,
            "transcription" => 0.001,
            "premrna_decay" => 0.001,
            "processing" => 0.02,
            "mrna_decay" => 6.4e-5,
            "translation" => 7.2e-9,
            "protein_decay" => 2.4e-11,
        ),
        at = 8.5, ribosomes = 2.0e6),
    "grs-defaults" => (;
        base_rates = Dict{String, Float64}(
            "activation" => 2.5,
            "deactivation" => 10.0,
            "initiation" => 6.6e-7,
            "abortion" => 0.01,
            "transcription" => 0.001,
            "premrna_decay" => 0.001,
            "processing" => 0.02,
            "mrna_decay" => 0.001,
            "translation" => 2.5e-9,
            "protein_decay" => 3e-10,
        ),
        at = 2.0, ribosomes = 2.0e6),
)

const SPECIES = Dict(
    "equilibrium" => ["elongations", "premrnas", "mrnas", "proteins"],
    "full" => ["active", "elongations", "premrnas", "mrnas", "proteins"],
)

settings = @add_arg_table! ArgParseSettings() begin
    "--out"
        required = true
    "--k"
        arg_type = Int
        default = 11
    "--switching-rate"
        dest_name = "switching_rate"
        arg_type = Float64
    "--variant"
        default = "full"
        range_tester = v -> haskey(SPECIES, v)
    "--method"
        default = "RSSACR"
    "--trajectories"
        arg_type = Int
        default = 100
    "--duration"
        arg_type = Float64
        default = 172800.0
    "--samples"
        arg_type = Int
        default = 10
    "--seed"
        default = "1"
    "--unique"
        arg_type = Bool
        default = true
    "--copies"
        arg_type = Int
        default = 1
    "--rates"
        default = "kronecker"
        range_tester = v -> haskey(MODELS, v)
    "--profile-reactions"
        dest_name = "profile_reactions"
        arg_type = Bool
        default = false
end

args = parse_args(settings, as_symbols = true)

genes = 2^args[:k]
model = MODELS[args[:rates]]

switching(rates, ::Nothing) = rates
switching(rates, d::Real) = let s = d / (rates["activation"] + rates["deactivation"])
    merge(rates, Dict("activation" => s * rates["activation"],
                      "deactivation" => s * rates["deactivation"]))
end

rates = switching(model.base_rates, args[:switching_rate])

adjacency(initiator, links) = Dict(
    "adjacency" => Dict("initiator" => initiator,
                        "power" => Dict("\$" => "k"),
                        "expected_links_per_gene" => links),
    "at" => ["LogNormal", model.at, 1.0],
)

initial = Dict{String, Any}("polymerases" => 5.0e5,
                          "ribosomes" => model.ribosomes,
                          "proteasomes" => 1.0e6)
if !args[:unique]
    width = length(string(genes))
    for i in 1:genes
        initial[string(lpad(i, width, '0'), ".inactive")] = args[:copies]
    end
end

schedule = Dict(
    "method" => args[:method],
    "k" => args[:k],
    "num_reps" => args[:trajectories],
    "duration" => args[:duration],
    "sample_interval" => args[:duration] / args[:samples],
    "step" => Dict(
        "do" => Dict("{regulation/kronecker}" => Dict(
            "seed" => "seed-$(args[:seed])",
            "method" => Dict("\$" => "method"),
            "species" => SPECIES[args[:variant]],
            "unique" => args[:unique],
            "base_rates" => rates,
            "profile_reactions" => args[:profile_reactions],
            "activation" => adjacency([[0.8, 0.3], [0.2, 0.5]], 0.8),
            "repression" => adjacency([[0.3, 0.6], [0.4, 0.1]], 1.2),
        )),
        "step" => [
            Dict("{add}" => initial),
            Dict("branch" => true, "step" => Dict(
                "each" => Dict("start" => 1, "stop" => Dict("\$" => "num_reps")),
                "step" => Dict("to" => Dict("\$" => "duration"),
                               "step" => Dict("\$" => "sample_interval")))),
        ],
    ),
)

mkpath(dirname(args[:out]))
open(io -> JSON.print(io, schedule, 2), args[:out], "w")

println(args[:out], ": 2^", args[:k], " = ", genes, " genes, ", args[:rates], " ", args[:variant],
        " (on ", round(1 / rates["deactivation"], digits = 1), "s, off ",
        round(1 / rates["activation"], digits = 1), "s), ",
        args[:trajectories], " traj, to=", args[:duration])
