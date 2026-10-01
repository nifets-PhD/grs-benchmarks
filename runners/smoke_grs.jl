using GeneRegulatorySystems
const GRS = GeneRegulatorySystems

example = joinpath(pkgdir(GeneRegulatorySystems),
                   "examples", "specification", "minimal.schedule.json")

schedule! = GRS.Models.load(example; seed="1234")
x = GRS.Models.FlatState(schedule!())

counts = x.counts
get_of(suffix) = sum(v for (k, v) in counts if endswith(String(k), suffix); init = 0)

println("grs $(pkgversion(GeneRegulatorySystems)): minimal.schedule.json, t = $(x.t)")
println("  elongations = $(counts[Symbol("1.elongations")])\texpected ~6")
println("  premrnas    = $(counts[Symbol("1.premrnas")])\texpected ~0.3")
println("  mrnas       = $(counts[Symbol("1.mrnas")])\texpected ~6")
println("  proteins    = $(counts[Symbol("1.proteins")])\texpected ~95")
