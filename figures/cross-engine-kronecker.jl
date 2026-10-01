using CairoMakie
using CairoMakie.GeometryBasics
using CSV
using DataFrames
using Statistics

load(f) = CSV.read("measurements/$f.csv", DataFrame)
cost, accuracy, tradeoff = load("kronecker-scaling"), load("kronecker-accuracy"),
                           load("kronecker-accuracy-cost")

curve(arm) = sort!(filter(:arm => ==(arm), cost), :genes)
runtime = Dict(tradeoff.arm .=> tradeoff.trajectory_runtime)
bias = Dict(accuracy.arm .=> accuracy.signed_smd)
nulls = filter(:arm => ==("null"), accuracy)

const GRS, GILLESPIE2, COPASI = Makie.wong_colors()[[1, 3, 6]]
const EXACT, HYBRID, LEAP = :rect, :diamond, :utriangle
const STYLE = Dict(EXACT => :solid, HYBRID => :dash, LEAP => :dot)
markersize(m) = m === HYBRID ? Vec2f(8, 10) : 10

ladders = [(GRS, EXACT, "grs-ssa"), (GRS, HYBRID, "grs-hybrid-tau3"),
           (GRS, LEAP, "grs-tauleap-tau3"),
           (GILLESPIE2, EXACT, "gillespiessa2-ssa"),
           (GILLESPIE2, LEAP, "gillespiessa2-tauleap-tau3"),
           (COPASI, EXACT, "copasi-ssa"), (COPASI, HYBRID, "copasi-hybrid-eps0.05"),
           (COPASI, LEAP, "copasi-tauleap-eps0.05")]
engines = [(GRS, "grs-ssa"), (GILLESPIE2, "gillespiessa2-tauleap-tau3"),
           (COPASI, "copasi-tauleap-eps0.05")]

plain(v) = [x >= 1 ? string(Int(round(x))) :
            rstrip(rstrip(string(round(x, sigdigits = 2)), '0'), '.') for x in v]

set_theme!(Theme(
    font = "TeX Gyre Heros Makie", fontsize = 14, backgroundcolor = :white,
    Axis = (backgroundcolor = :white, xtickformat = plain, ytickformat = plain,
            xgridcolor = (:black, 0.12), ygridcolor = (:black, 0.12),
            xgridwidth = 0.8, ygridwidth = 0.8,
            xminorgridvisible = false, yminorgridvisible = false,
            topspinevisible = true, rightspinevisible = true, spinewidth = 1.1,
            xtickalign = 1, ytickalign = 1, xticksize = 4, yticksize = 4,
            xlabelfont = :bold, ylabelfont = :bold, xlabelsize = 16, ylabelsize = 16,
            xticklabelsize = 14, yticklabelsize = 14,
            xticklabelfont = :bold, yticklabelfont = :bold),
    Legend = (framevisible = false, labelsize = 13, labelfont = :bold,
              patchsize = (16, 14))))

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

const WIDTH, HEIGHT = 980, 820
fig = Figure(size = (WIDTH, HEIGHT), figure_padding = 12)
card = poly!(fig.scene, roundedrect(6, 6, WIDTH - 12, HEIGHT - 12, 20);
             color = RGBf(0.94, 0.94, 0.94), space = :pixel)
translate!(card, 0, 0, -100)

shared = (; xscale = log2, xticks = 2.0 .^ (3:2:15), xlabel = "# Genes",
          backgroundcolor = :white)

ax1 = Axis(fig[1, 1]; yscale = log10, yticks = 10.0 .^ (-2:4),
           ylabel = "Runtime per trajectory (s)", shared...)
ax2 = Axis(fig[1, 2]; yscale = log10, yticks = 10.0 .^ (-2:4),
           ylabel = "Build time (s)", shared...)
ax3 = Axis(fig[2, 1]; yscale = log10, yticks = vec([m * 10.0^e for m in (1, 2, 5), e in -1:0]),
           ylabel = "Peak Memory (GiB)", shared...)

for (colour, marker, arm) in ladders
    d = curve(arm)
    scatterlines!(ax1, d.genes, d.trajectory_runtime;
        color = colour, marker, markersize = markersize(marker), linewidth = 1.6,
        linestyle = STYLE[marker])
end
build(row) = row.engine == "grs" ? row.warmup_build_seconds : row.build_seconds
for (colour, arm) in engines
    d = curve(arm)
    b = filter(:build_seconds => !ismissing, d)
    scatterlines!(ax2, b.genes, build.(eachrow(b));
        color = colour, markersize = 7, linewidth = 1.6)
    scatterlines!(ax3, d.genes, d.rss_peak_bytes ./ 2^30;
        color = colour, markersize = 7, linewidth = 1.6)
end

floor = quantile(abs.(nulls.signed_smd), 0.95)
const TAUS = ("10", "3", "1", "0.3", "0.1")
point(arm) = (bias[arm], runtime[arm])
ladderof(prefix) = [(point("$prefix$t")..., t) for t in TAUS]
trails = [
  (GRS, LEAP, ladderof("grs-tauleap-tau"), :gray15),
  (GRS, HYBRID, ladderof("grs-hybrid-tau"), GRS),
  (GILLESPIE2, LEAP, ladderof("gillespiessa2-tauleap-tau"), nothing),
  (GILLESPIE2, EXACT, [(point("gillespiessa2-ssa")..., "")], nothing),
  (COPASI, HYBRID, [(point("copasi-hybrid-eps0.05")..., "")], nothing),
  (COPASI, EXACT, [(point("copasi-ssa")..., "")], nothing),
  (COPASI, LEAP, [(point("copasi-tauleap-eps0.05")..., "")], nothing),
  (GRS, EXACT, [(0.0, runtime["grs-ssa"], "")], nothing),
]

ax4 = Axis(fig[2, 2]; xscale = Makie.Symlog10(floor), yscale = log10,
    backgroundcolor = :white, xlabel = "Mean SMD",
    ylabel = "Runtime per trajectory (s)", xtickformat = Makie.automatic,
    xticks = ([-0.01, 0, 0.01, 0.1, 1],
              ["−0.01", "0", "0.01", "0.1", "1"]))
vspan!(ax4, -floor, floor; color = (:black, 0.12))
vlines!(ax4, 0; color = (:black, 0.35), linewidth = 1, linestyle = :dash)
for (colour, marker, points, _) in trails
    x, y = [p[1] for p in points], [p[2] for p in points]
    length(points) > 1 && lines!(ax4, x, y;
        color = (colour, 0.6), linewidth = 1.6, linestyle = STYLE[marker])
    scatter!(ax4, x, y; color = colour, marker, markersize = markersize(marker) .* 1.2)
end
const TAULABELS = [
  ("grs-hybrid-tau0.1", GRS, (:right, :center), (-9, 0)),
  ("grs-hybrid-tau10", GRS, (:left, :bottom), (-2, 8)),
  ("grs-tauleap-tau0.1", GRS, (:left, :top), (7, -3)),
  ("grs-tauleap-tau10", GRS, (:right, :center), (-8, -1)),
  ("gillespiessa2-tauleap-tau0.1", GILLESPIE2, (:right, :bottom), (-6, 4)),
  ("gillespiessa2-tauleap-tau10", GILLESPIE2, (:right, :bottom), (-6, 3)),
]
for (arm, align, offset) in (("copasi-tauleap-eps0.05", (:center, :top), (0, -9)),
                             ("copasi-hybrid-eps0.05", (:left, :center), (8, 0)))
    text!(ax4, point(arm)...; text = "ε=0.05", align, offset,
          fontsize = 11, font = :bold, color = COPASI)
end
for (arm, colour, align, offset) in TAULABELS
    text!(ax4, point(arm)...; text = "τ=" * last(split(arm, "tau")), align, offset,
          fontsize = 11, font = :bold, color = colour)
end
text!(ax4, 0.0, runtime["grs-ssa"]; text = "reference", align = (:right, :center),
      offset = (-10, 0), fontsize = 11, font = :bold, color = GRS)

for (i, (row, col)) in enumerate([(1, 1), (1, 2), (2, 1), (2, 2)])
    Label(fig[row, col, TopLeft()], string('A' + i - 1);
          font = :bold, fontsize = 22, padding = (0, 8, 4, 0), halign = :right)
end

engine_key = [MarkerElement(color = c, marker = :circle, markersize = 16)
              for c in (GRS, GILLESPIE2, COPASI)]
method_key = [[MarkerElement(color = :gray35, marker = m, markersize = markersize(m)),
               LineElement(color = :gray35, linestyle = STYLE[m])]
              for m in (EXACT, HYBRID, LEAP)]
Legend(fig[3, 1:2], engine_key,
       ["GeneRegulatorySystems.jl", "GillespieSSA2 (dyngen)", "COPASI"];
       orientation = :horizontal)
Legend(fig[4, 1:2], method_key,
       ["SSA", "Hybrid SSA/Tau-leaping", "Tau-leaping"];
       orientation = :horizontal)

linkyaxes!(ax1, ax2)

rowgap!(fig.layout, 12)
rowgap!(fig.layout, 3, 0)
colgap!(fig.layout, 18)
rowsize!(fig.layout, 1, Auto(1.25))
save("figures/cross-engine-kronecker.pdf", fig)
