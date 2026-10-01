using CSV, DataFrames

const RUN = r"genes-0*(\d+)/([^/]+)-n(\d+)-s(\d+)$"
const ARM = r"^(grs|gillespiessa2|copasi|dyngen)-(ssa|hybrid|tauleap|tausplitting)(?:-(tau|eps)([0-9.]+))?(?:-thin|-qss|-tausplitting)?$"

const LEADING = [:genes, :arm, :engine, :method, :engine_method, :tau, :epsilon,
                 :trajectories, :seed, :model, :path,
                 :total_runtime, :trajectory_runtime, :build_seconds, :warmup_seconds,
                 :rss_before_build_bytes, :rss_after_build_bytes, :rss_peak_bytes,
                 :aborted_trajectories]

entries(path, sep) = Dict{Symbol, Any}(Symbol(strip(k)) => String(strip(v)) for (k, v) in
    (split(line, sep, limit = 2) for line in eachline(path) if occursin(sep, line)))

function identify(dir)
    r = match(RUN, dir)
    r === nothing && error("unparseable run path: $dir")
    a = match(ARM, r[2])
    a === nothing && error("unparseable arm: $(r[2])")
    (; genes = parse(Int, r[1]), arm = r[2], engine = a[1], method = a[2],
       tau = a[3] == "tau" ? parse(Float64, a[4]) : missing,
       epsilon = a[3] == "eps" ? parse(Float64, a[4]) : missing,
       trajectories = parse(Int, r[3]), seed = parse(Int, r[4]), path = dir)
end

function collate(dirs)
    df = DataFrame()
    for dir in dirs
        isfile(joinpath(dir, "metrics.csv")) &&
            isfile(joinpath(dir, "config.txt")) || continue
        row = merge(entries(joinpath(dir, "config.txt"), "="),
                    entries(joinpath(dir, "metrics.csv"), ","))
        delete!(row, :metric)
        haskey(row, :method) && (row[:engine_method] = pop!(row, :method))
        id = identify(dir)
        for k in (:trajectories, :seed)
            !haskey(row, k) || parse(Int, row[k]) == id[k] ||
                error("$k mismatch in $dir: path $(id[k]), config $(row[k])")
        end
        total_runtime = parse(Float64, row[:simulate_seconds])
        trajectory_runtime = haskey(row, :simulate_seconds_per_trajectory) ?
            parse(Float64, row[:simulate_seconds_per_trajectory]) :
            total_runtime / id.trajectories
        push!(df, (; row..., id..., total_runtime, trajectory_runtime); cols = :union)
    end
    select!(df, intersect(LEADING, propertynames(df)), :)
    key = [:genes, :arm, :trajectories, :seed, :model]
    clashes = df[nonunique(df, key), key]
    nrow(clashes) == 0 || error("duplicate configs:\n$clashes")
    sort!(df, [:genes, :engine, :method, :tau])
end

if abspath(PROGRAM_FILE) == @__FILE__
    df = collate(ARGS[2:end])
    nrow(df) == 0 && error("no completed runs among $(length(ARGS) - 1) planned")
    mkpath(dirname(ARGS[1]))
    CSV.write(ARGS[1], df)
    println(stderr, "$(nrow(df))/$(length(ARGS) - 1) runs -> $(ARGS[1])")
end
