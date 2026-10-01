include(joinpath(@__DIR__, "CatalystSBML.jl"))
include(joinpath(@__DIR__, "CatalystGillespie.jl"))

using GeneRegulatorySystems, Catalyst, ModelingToolkit, JSON
const CS = CatalystSBML
const CG = CatalystGillespie

const M = GeneRegulatorySystems.Models
const SC = M.SciML
const MT = ModelingToolkit

const POOLS = ["polymerases", "ribosomes", "proteasomes"]

v1model(x) = v1model(x, x.definition)
v1model(x, ::M.V1.Definition) = x
v1model(x, ::Any) = v1model(x.model)

findkey(x, key) = nothing
function findkey(x::AbstractVector, key)
    for e in x
        r = findkey(e, key)
        isnothing(r) || return r
    end
    nothing
end
function findkey(x::AbstractDict, key)
    haskey(x, key) && return x[key]
    for (_, v) in x
        r = findkey(v, key)
        isnothing(r) || return r
    end
    nothing
end

write_model(::Val{:sbml}, definition, out; tag, kwargs...) =
    CS.writesbml(definition, out; modelid = tag, kwargs...)

write_model(::Val{:gillespie}, definition, out; tag, kwargs...) =
    CG.writegillespie(definition, out; kwargs...)

function convert_schedule(target, schedule, out; tag = first(split(basename(out), '.')),
                          freeze_pools = true, seed = "1", path = "+.do")
    sched = M.load(schedule; seed = seed)
    add = something(findkey(JSON.parsefile(schedule, dicttype = Dict{Symbol, Any}),
                            Symbol("{add}")), Dict{Symbol, Any}())

    f! = v1model(M.Scheduling.reify(sched, path))
    rs = Catalyst.flatten(f!.model.definition)
    jm = f!.model.model

    pv = Dict{Any, Float64}(MT.unwrap(p) => float(jm[p])
                            for p in SC.parameter_symbols(jm))
    boot = Dict{Symbol, Int}(k => Int(v) for (k, v) in M.load_defaults()[:bootstrap])
    merge!(boot, Dict{Symbol, Int}(k => Int(v) for (k, v) in add))
    ia = Dict{Any, Float64}(MT.unwrap(s) => float(get(boot, SC.normalize_name(s), 0))
                            for s in MT.unknowns(rs))

    mkpath(dirname(out))
    write_model(Val(Symbol(target)), f!.model.definition, out;
                tag = tag,
                parameter_values = pv,
                initial_amounts = ia,
                constant_species = freeze_pools ? POOLS : String[],
                ratelaw = CS.foldedratelaw(pv, tag))
    out
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) >= 4 ||
        error("usage: julia --project=. models/export.jl <sbml|gillespie> " *
              "<schedule.json> <out> <seed> [freeze_pools]")
    freeze = length(ARGS) >= 5 ? parse(Bool, ARGS[5]) : true
    out = convert_schedule(ARGS[1], ARGS[2], ARGS[3]; seed = ARGS[4], freeze_pools = freeze)
    println("wrote $out  $(round(filesize(out) / 1e6, digits = 2)) MB")
end
