using Arrow
using CSV
using DataFrames

rows = DataFrame()
for dir in ARGS[3:end]
    t = DataFrame(Arrow.Table(joinpath(dir, "counts.arrow")))
    t = t[(t.sample .== maximum(t.sample)) .& endswith.(String.(t.name), ARGS[2]), :]
    append!(rows, DataFrame(run = dir, name = String.(t.name), value = t.value))
end
mkpath(dirname(ARGS[1]))
CSV.write(ARGS[1], rows)
println(stderr, nrow(rows), " final counts from ", length(ARGS) - 2, " runs -> ", ARGS[1])
