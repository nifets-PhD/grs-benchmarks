using Arrow
using CSV
using DataFrames
using Statistics

const PROBES = ("initiated", "aborted", "transcribed", "premrnas_decayed", "processed",
    "mrnas_decayed", "translated", "proteins_decayed", "activated", "deactivated")

kind(name) = last(split(name, '.'; limit = 2))
probe(name) = occursin('.', name) &&
    (kind(name) in PROBES || startswith(kind(name), "proteolyzed_by_"))

function events(dir)
    t = DataFrame(Arrow.Table(joinpath(dir, "counts.arrow")))
    t.name = String.(t.name)
    final = filter(r -> r.sample == maximum(t.sample) && probe(r.name), t)
    final.kind = kind.(final.name)
    genes = parse(Int, match(r"genes-0*(\d+)", dir)[1])
    seconds = maximum(final.t)
    totals = combine(groupby(final, [:path, :kind]), :value => sum => :total)
    combine(groupby(totals, :kind),
        :total => (v -> mean(v) / genes / seconds) => :rate,
        :path => (p -> length(unique(p))) => :trajectories)
end

rows = DataFrame()
for dir in ARGS[2:end]
    append!(rows, insertcols!(events(dir), 1, :path => dir); cols = :union)
end
mkpath(dirname(ARGS[1]))
CSV.write(ARGS[1], rows)
println(stderr, nrow(rows), " event rates from ", length(ARGS) - 1, " runs -> ", ARGS[1])
