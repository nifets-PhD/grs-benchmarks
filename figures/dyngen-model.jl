using CSV
using DataFrames
using Printf
using Statistics

const ROWS = [
    "dyngen-tauleap-tau0.0833" => ("dyngen", raw"tau-leaping, fixed $\tau$, clamped negatives"),
    "grs-tauleap-tau0.0833" => (raw"\thetool{}", raw"tau-leaping, fixed $\tau$, clamped negatives"),
    "grs-tauleap-tau0.0833-thin" => ("", raw"tau-leaping, fixed $\tau$, thinned negatives"),
    "grs-tauleap-eps0.05" => ("", raw"adaptive tau-leaping, $\epsilon = 0.05$"),
    "grs-hybrid-tau0.0833-tausplitting" => ("", raw"hybrid, fixed $\tau$, threshold $n_c = 10$"),
    "grs-tausplitting" => ("", "exact SSA (TauSplitting)"),
]

df = CSV.read("measurements/dyngen-model.csv", DataFrame)

build(row) = row.engine == "grs" ? row.warmup_build_seconds : row.build_seconds

groups = Vector{Vector{Pair{String, Tuple{String, String}}}}()
for entry in ROWS
    isempty(first(last(entry))) ? push!(last(groups), entry) : push!(groups, [entry])
end

open("figures/dyngen-model.tex", "w") do io
    println(io, raw"\begin{tabular}{lr|lr}")
    println(io, raw"\hline")
    println(io, raw"\textbf{Engine} & \textbf{Build (s)} & \textbf{Method} & \textbf{Runtime (s/traj.)} ", "\\\\")
    for group in groups
        rows = [only(eachrow(filter(:arm => ==(arm), df))) for (arm, _) in group]
        seconds = median(build.(rows))
        println(io, raw"\hline")
        for (i, ((_, (engine, method)), row)) in enumerate(zip(group, rows))
            cells = i == 1 ? (engine, @sprintf("%.1f", seconds)) : ("", "")
            @printf(io, "%s & %s & %s & %.1f \\\\\n", cells..., method, row.trajectory_runtime)
        end
    end
    println(io, raw"\hline")
    println(io, raw"\end{tabular}")
end
