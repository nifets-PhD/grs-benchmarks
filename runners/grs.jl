using ArgParse
using JSON
using Logging
using GeneRegulatorySystems: Models, Specifications
const Scheduling = Models.Scheduling
import Arrow

settings = @add_arg_table! ArgParseSettings() begin
    "out"
        required = true
        help = "output path prefix."
    "schedule"
        required = true
        help = "schedule spec to run"
    "--seed"
        default = "1"
    "--model"
        default = "+.do"
        help = "path to the regulation model to prebuild"
    "--parallel"
        action = :store_true
    "--method"
        help = "override the aggregator defined in the schedule"
    "--tau"
        arg_type = Float64
    "--epsilon"
        arg_type = Float64
    "--quantile"
        arg_type = Float64
        default = 0.03
    "--tau-cap"
        dest_name = "tau_cap"
        arg_type = Float64
        default = Inf
        help = "upper bound on AdaptiveTau step"
    "--policy"
        default = "default"
        help = "blending policy, or \"default\" to let the definition choose"
    "--nc"
        arg_type = Int
        help = "criticality threshold for CriticalBlend; omit to use the policy default"
    "--exact"
        default = "RSSACR"
        help = "inner exact aggregator for HybridTau"
    "--thin-negative"
        dest_name = "thin_negative"
        arg_type = Bool
        default = true
        help = "false clamps negatives to zero like dyngen's ssa_etl"
    "--trajectories"
        arg_type = Int
        help = "override num_reps from the schedule"
end

args = parse_args(settings, as_symbols = true)

paths, times, names, values = String[], Float64[], Symbol[], Int64[]

guard = ReentrantLock()

function record(x; path, into = nothing, _...)
    into === nothing && return nothing
    t = Models.t(x)
    counts = Models.FlatState(x).counts
    lock(guard) do
        for (name, value) in counts
            push!(paths, string(path))
            push!(times, t)
            push!(names, name)
            push!(values, Int64(value))
        end
    end
    nothing
end

started = time()
report(stage) = (println("[", round(Int, time() - started), "s] ", stage);
                 flush(stdout))

specification = JSON.parsefile(args[:schedule])
interval = specification["sample_interval"]
trajectories = something(args[:trajectories], specification["num_reps"])

policy!(spec) = let
    s = args[:policy] == "default" ? spec :
        merge(spec, Dict{Symbol, Any}(:policy => args[:policy]))
    args[:nc] === nothing ? s : merge(s, Dict{Symbol, Any}(:nc => args[:nc]))
end

required(flag, value, method) = value === nothing ?
    error("--$flag is required for $method") : value

aggregator(::Val{:AdaptiveTau}) = policy!(Dict{Symbol, Any}(
    :name => "HybridTau", :exact => args[:exact],
    :dt => args[:tau_cap],
    :epsilon => required("epsilon", args[:epsilon], "AdaptiveTau"),
    :quantile => args[:quantile], :thin_negative => args[:thin_negative]))

aggregator(::Val{:FixedTau}) = policy!(Dict{Symbol, Any}(
    :name => "HybridTau", :exact => args[:exact],
    :dt => required("tau", args[:tau], "FixedTau"), :epsilon => nothing,
    :thin_negative => args[:thin_negative]))

aggregator(::Val{T}) where {T} = String(T)

spec = args[:method] === nothing ? nothing : aggregator(Val(Symbol(args[:method])))

override(f!, name, ::Nothing) = f!
override(f!, name, value) =
    Specifications.set(f!, ".$name", Specifications.Template(value))

function build(; num_reps = trajectories, duration = nothing)
    f! = Models.load(args[:schedule]; seed = args[:seed])
    f! = override(f!, :method, spec)
    f! = override(f!, :num_reps, num_reps)
    f! = override(f!, :duration, duration)
    Scheduling.prebuild(f!, args[:model])
end

report("warming up")
t_warmup_build = @elapsed warmup! = build(num_reps = 1, duration = interval)
t_warmup_simulate = @elapsed warmup!(Models.FlatState(), interval;
    trace = record, parallel = args[:parallel])
t_warmup = t_warmup_build + t_warmup_simulate
foreach(empty!, (paths, times, names, values))

rss_before_build = Sys.maxrss()
report("building")
t_build = @elapsed schedule! = build()
rss_built = Sys.maxrss()
report("simulating")
t_simulate = @elapsed schedule!(; trace = record, parallel = args[:parallel])
rss_peak = Sys.maxrss()

report("writing")
mkpath(args[:out])

open(joinpath(args[:out], "config.txt"), "w") do io
    println(io, "engine=grs")
    println(io, "method=", args[:method])
    println(io, "grs_policy=", args[:policy])
    println(io, "trajectories=", trajectories)
    println(io, "seed=", args[:seed])
    println(io, "model=", args[:schedule])
    spec isa AbstractDict || return
    for key in sort(collect(keys(spec)))
        spec[key] === nothing || println(io, "grs_", key, "=", spec[key])
    end
end

Arrow.write(joinpath(args[:out], "counts.arrow"),
    (; path = paths, sample = round.(Int32, times ./ interval),
       t = times, name = names, value = values))

open(joinpath(args[:out], "metrics.csv"), "w") do io
    println(io, "metric,value")
    println(io, "build_seconds,", t_build)
    println(io, "warmup_seconds,", t_warmup)
    println(io, "warmup_build_seconds,", t_warmup_build)
    println(io, "simulate_seconds,", t_simulate)
    println(io, "simulate_seconds_per_trajectory,", t_simulate / trajectories)
    println(io, "rss_before_build_bytes,", rss_before_build)
    println(io, "rss_after_build_bytes,", rss_built)
    println(io, "rss_peak_bytes,", rss_peak)
end

println("grs $(args[:out]): ", round(t_simulate, digits = 1), " s simulate, ",
        round(t_build, digits = 1), " s build, ",
        round(rss_peak / 2^30, digits = 2), " GiB peak rss, ",
        length(values), " count rows")
