using CairoMakie
using CSV
using DataFrames

const ARM = "grs-ssa"

df = CSV.read("measurements/kronecker-scaling.csv", DataFrame)
df = sort!(filter(:arm => ==(ARM), df), :genes)

plain(v) = [x >= 1 ? string(Int(round(x))) :
            rstrip(rstrip(string(round(x, sigdigits = 2)), '0'), '.') for x in v]

genes = df.genes
shared = (; xscale = log2, xticks = float.(genes), xtickformat = plain)

fig = Figure(size = (480, 580))

ax1 = Axis(fig[1, 1]; yscale = log10, yticks = 10.0 .^ (-2:4), ytickformat = plain,
    yminorticks = IntervalsBetween(9), yminorgridvisible = true,
    yminorgridcolor = (:black, 0.08),
    ylabel = "wall clock (s)", title = "runtime", shared...)

scatterlines!(ax1, genes, df.build_seconds;
    label = "build time", color = Makie.wong_colors()[1], markersize = 8)
scatterlines!(ax1, genes, df.trajectory_runtime;
    label = "time per trajectory", color = Makie.wong_colors()[2], markersize = 8)

w = filter(:warmup_seconds => !ismissing, df)
scatterlines!(ax1, w.genes, w.warmup_seconds;
    label = "warmup", color = Makie.wong_colors()[6], linestyle = :dash,
    markersize = 8)

axislegend(ax1; position = :lt)

ax2 = Axis(fig[2, 1]; xticklabelrotation = pi / 2,
    xlabel = "# genes", ylabel = "GiB", title = "peak memory",
    limits = (nothing, (0, nothing)), shared...)

scatterlines!(ax2, genes, df.rss_peak_bytes ./ 2^30;
    color = Makie.wong_colors()[3], markersize = 8)

linkxaxes!(ax1, ax2)
hidexdecorations!(ax1; grid = false)
rowgap!(fig.layout, 10)

out = "figures/scaling-kronecker-full.pdf"
save(out, fig)
save("/tmp/fig-scaling.png", fig; px_per_unit = 2)
println("wrote ", out)
