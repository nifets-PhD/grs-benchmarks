using CairoMakie
using CairoMakie.GeometryBasics
using CSV
using DataFrames
using Statistics

const DIR = "measurements/promoter-approximation"

rate_of(path) = parse(Float64, match(r"-d([0-9.e-]+)(?:-profiled)?/", path)[1])
variant_of(path) = occursin("/equilibrium-", path) ? "equilibrium" : "full"
label(rate) = rate >= 1 ? string(Int(rate)) : string(rate)

costs = CSV.read(joinpath(DIR, "costs.csv"), DataFrame)
events = CSV.read(joinpath(DIR, "events.csv"), DataFrame)
events.switching = rate_of.(events.path)
runtime(variant, rate) = only(costs.trajectory_runtime[
    (variant_of.(costs.path) .== variant) .& (rate_of.(costs.path) .== rate)])
fired(rate, kind) = only(events.rate[(events.switching .== rate) .& (events.kind .== kind)])

rows = DataFrame()
for file in filter(f -> startswith(f, "distances-d"), readdir(DIR))
    d = CSV.read(joinpath(DIR, file), DataFrame)
    rate = parse(Float64, match(r"distances-d(.+)\.csv", file)[1])
    push!(rows, (; rate, hist = only(d.hist[d.arm .== "grs-ssa-qss"]),
        floor = mean(d.hist[d.arm .== "null"]),
        switching_runtime = runtime("full", rate), qss_runtime = runtime("equilibrium", 1.0)))
end
sort!(rows, :rate; rev = true)

const PROMOTER = [("activated", "activation"), ("deactivated", "deactivation")]
const CASCADE = [("initiated", "initiation"), ("aborted", "abortion"), ("transcribed", "elongation"),
    ("premrnas_decayed", "pre-mRNA decay"), ("processed", "processing"),
    ("mrnas_decayed", "mRNA decay"), ("translated", "translation"),
    ("proteins_decayed", "protein decay")]
const REACTIONS = [PROMOTER; CASCADE]

const GENE = "014"
function example(species)
    d = filter(:name => ==("$GENE.$species"), CSV.read(joinpath(DIR, "$species.csv"), DataFrame))
    d.arm = replace.(d.run, r".*/((?:full|equilibrium)-d[0-9.]+)/.*" => s"\1")
    d
end

set_theme!(Theme(
    font = "TeX Gyre Heros Makie", fontsize = 14, backgroundcolor = :white,
    Axis = (backgroundcolor = :white,
            xgridcolor = (:black, 0.12), ygridcolor = (:black, 0.12),
            xgridwidth = 0.8, ygridwidth = 0.8, spinewidth = 1.1,
            xtickalign = 1, ytickalign = 1, xticksize = 4, yticksize = 4,
            xlabelfont = :bold, ylabelfont = :bold, xlabelsize = 16, ylabelsize = 16,
            xticklabelsize = 14, yticklabelsize = 14,
            xticklabelfont = :bold, yticklabelfont = :bold,
            width = 380, height = 200),
    Legend = (framevisible = false, labelsize = 13, labelfont = :bold, patchsize = (16, 14))))

function roundedrect(x, y, w, h, r; n = 16)
    pts = Point2f[]
    for (cx, cy, a0) in ((x + w - r, y + h - r, 0.0), (x + r, y + h - r, pi / 2),
                         (x + r, y + r, pi), (x + w - r, y + r, 3pi / 2))
        for t in range(a0, a0 + pi / 2; length = n)
            push!(pts, Point2f(cx + r * cos(t), cy + r * sin(t)))
        end
    end
    Polygon(pts)
end

const QSS, TELEGRAPH = colorant"#4063D8", colorant"#E6A817"
const FAST, SLOW = colorant"#F5CD55", colorant"#E08E00"
const FLOOR = RGBf(0.45, 0.45, 0.45)
superscript(e) = join(Char(Dict('-' => 0x207B, '0' => 0x2070, '1' => 0x00B9, '2' => 0x00B2,
    '3' => 0x00B3)[c]) for c in string(e))
ticklabel(rate) = rate >= 1 ? label(rate) : "10" * superscript(round(Int, log10(rate)))
ticks = (rows.rate, ticklabel.(rows.rate))

fig = Figure(figure_padding = 24)
fast, slow = maximum(rows.rate), minimum(rows.rate)

a = Axis(fig[1, 1]; height = 130, ylabelsize = 14, xscale = log10, xticks = ticks, xticklabelsvisible = false,
    ylabel = "Histogram distance")
scatterlines!(a, rows.rate, rows.hist; color = QSS, label = "QSS vs. telegraph")
scatterlines!(a, rows.rate, rows.floor; color = FLOOR, linestyle = :dash, label = "Noise floor")
ylims!(a, 0, 2)
axislegend(a; position = :lt)
b = Axis(fig[2, 1]; height = 130, ylabelsize = 14, xscale = log10, xticks = ticks,
    xlabel = "Promoter switching rate (s⁻¹)", ylabel = "Runtime per trajectory (s)")
scatterlines!(b, rows.rate, rows.switching_runtime; color = TELEGRAPH, label = "Telegraph")
scatterlines!(b, rows.rate, rows.qss_runtime; color = QSS, label = "QSS")
axislegend(b; position = :lt)
linkxaxes!(a, b)
ylims!(b, 0, nothing)

c = Axis(fig[3, 1]; height = 150, xscale = log10, yreversed = true, ygridvisible = false,
    yticks = (eachindex(REACTIONS), last.(REACTIONS)), yticklabelsize = 11,
    xlabel = "Event rate per gene (s⁻¹)")
for (y, (k, _)) in enumerate(PROMOTER)
    lines!(c, [fired(r, k) for r in rows.rate], fill(y, nrow(rows)); color = TELEGRAPH)
    scatter!(c, [fired(r, k) for r in rows.rate if r != fast && r != slow],
        fill(y, nrow(rows) - 2); color = TELEGRAPH, markersize = 5)
end
for (rate, colour) in ((fast, FAST), (slow, SLOW))
    scatter!(c, [fired(rate, k) for (k, _) in PROMOTER], eachindex(PROMOTER);
        color = colour, markersize = 9, strokewidth = 0.8, strokecolor = :black)
end
present = [k in events.kind[events.switching .== fast] for (k, _) in CASCADE]
scatter!(c, [fired(fast, k) for (k, _) in CASCADE[present]], length(PROMOTER) .+ findall(present);
    color = FLOOR, markersize = 9, strokewidth = 0.8, strokecolor = :black)

cascade = fig[1:3, 2] = GridLayout()
for (i, (species, title)) in enumerate((("mrnas", "mRNA"), ("proteins", "Protein")))
    d = example(species)
    counts(arm) = d.value[d.arm .== arm]
    edges = range(0, quantile(d.value, 0.99); length = 31)
    ax = Axis(cascade[i, 1]; width = 220, height = Auto(), title = "Gene $GENE, $title",
        titlesize = 14, xlabel = "Counts", ylabel = "Fraction", xlabelsize = 14, ylabelsize = 14,
        xticks = WilkinsonTicks(3; k_min = 3, k_max = 4),
        xtickformat = v -> string.(round.(Int, v)))
    for (rate, colour) in ((fast, FAST), (slow, SLOW))
        hist!(ax, counts("full-d$(label(rate))"); bins = edges, normalization = :probability,
            color = (colour, 0.5), strokewidth = 0.6, strokecolor = colour,
            label = "Telegraph, $(ticklabel(rate)) s⁻¹")
    end
    stephist!(ax, counts("equilibrium-d1"); bins = edges, normalization = :probability,
        color = QSS, linewidth = 2, label = "QSS")
    xlims!(ax, 0, last(edges))
    ylims!(ax, 0, nothing)
end

Legend(cascade[3, 1], cascade.content[end].content; framevisible = false, labelsize = 12,
    patchsize = (12, 10), tellwidth = false, tellheight = true)

for (l, pos) in zip("ABCD", (fig[1, 1, TopLeft()], fig[2, 1, TopLeft()], fig[3, 1, TopLeft()],
        cascade[1, 1, TopLeft()]))
    Label(pos, string(l); font = :bold, fontsize = 20, padding = (0, 20, 8, 0))
end
resize_to_layout!(fig)
save("figures/promoter-approximation.pdf", fig)
