module CatalystSBML

using Catalyst, ModelingToolkit, Symbolics, SymbolicUtils
using SymbolicUtils: iscall, operation, arguments, maketerm
const MT = ModelingToolkit

export writesbml, foldedratelaw

const SBML_NS = "http://www.sbml.org/sbml/level3/version2/core"
const MATHML_NS = "http://www.w3.org/1998/Math/MathML"

subs_params(ex, rules) = MT.substitute(ex, rules)

function constant_value(ex)
    ex isa Number && return float(ex)
    isempty(Symbolics.get_variables(ex)) || return nothing
    v = try
        Base.eval(@__MODULE__, Symbolics.toexpr(ex))
    catch
        return nothing
    end
    v isa Number ? float(v) : nothing
end

function constant_bool(ex)
    ex isa Bool && return ex
    isempty(Symbolics.get_variables(ex)) || return nothing
    v = try
        Base.eval(@__MODULE__, Symbolics.toexpr(ex))
    catch
        return nothing
    end
    v isa Bool ? v : nothing
end

function no_negative_coefficient(ex)
    c = constant_value(ex)
    isnothing(c) || return c >= 0
    iscall(ex) || return true
    all(no_negative_coefficient, arguments(ex))
end

function positive_on_nonnegative(ex)
    terms = iscall(ex) && operation(ex) === (+) ? arguments(ex) : [ex]
    constant = 0.0
    for t in terms
        c = constant_value(t)
        isnothing(c) || (constant += c)
    end
    constant > 0 && all(no_negative_coefficient, terms)
end

function decidably_true(cond)
    iscall(cond) && operation(cond) === (<) || return false
    lo = constant_value(arguments(cond)[1])
    !isnothing(lo) && lo <= 0 && positive_on_nonnegative(arguments(cond)[2])
end

function fold_guards(ex)
    iscall(ex) || return ex
    args = map(fold_guards, arguments(ex))
    if operation(ex) === ifelse
        b = constant_bool(args[1])
        isnothing(b) || return b ? args[2] : args[3]
        decidably_true(args[1]) && return args[2]
    end
    maketerm(typeof(ex), operation(ex), args, SymbolicUtils.metadata(ex))
end

function assert_continuous(ex, tag)
    iscall(ex) || return
    operation(ex) === ifelse &&
        error("undecidable ifelse survives in $tag: $ex")
    foreach(a -> assert_continuous(a, tag), arguments(ex))
end

function foldedratelaw(parameter_values, tag = "model")
    function (rx)
        law = fold_guards(subs_params(Catalyst.oderatelaw(rx), parameter_values))
        assert_continuous(law, tag)
        law
    end
end

function sanitize(name::AbstractString)
    s = replace(name, r"\(t\)$" => "")
    s = replace(s, '₊' => '_')
    s = replace(s, r"[^A-Za-z0-9_]" => "_")
    isempty(s) && (s = "id")
    occursin(r"^[A-Za-z_]", s) || (s = "_" * s)
    s
end

function idmap(syms)
    out = Dict{Any,String}()
    seen = Dict{String,Int}()
    for s in syms
        base = sanitize(string(s))
        n = get(seen, base, 0)
        seen[base] = n + 1
        out[MT.unwrap(s)] = n == 0 ? base : "$(base)_$n"
    end
    out
end

xmlescape(s) = replace(string(s), "&" => "&amp;", "<" => "&lt;", ">" => "&gt;",
                       "\"" => "&quot;")

num(x::Integer) = "<cn type=\"integer\"> $x </cn>"
function num(x::Real)
    f = float(x)
    isinteger(f) && abs(f) < 1e15 && return "<cn type=\"integer\"> $(Int(f)) </cn>"
    s = string(f)
    if occursin('e', s) || occursin('E', s)
        m, e = split(replace(s, 'E' => 'e'), 'e')
        return "<cn type=\"e-notation\"> $m <sep/> $(parse(Int, e)) </cn>"
    end
    "<cn> $s </cn>"
end

const OPS = Dict(:+ => "plus", :- => "minus", :* => "times", :/ => "divide",
                 :^ => "power", :abs => "abs", :min => "min", :max => "max",
                 :exp => "exp", :log => "ln", :sqrt => "root", :sin => "sin",
                 :cos => "cos", :tan => "tan",
                 :< => "lt", :> => "gt", :<= => "leq", :>= => "geq",
                 :(==) => "eq", :!= => "neq", :& => "and", :| => "or", :! => "not")

function mathml(ex, ids)
    ex isa Number && return num(ex)
    u = MT.unwrap(ex)
    haskey(ids, u) && return "<ci> $(ids[u]) </ci>"
    u isa Number && return num(u)
    if !iscall(u)
        v = try
            Symbolics.toexpr(u)
        catch
            nothing
        end
        v isa Number && return num(v)
        error("no SBML id for symbol $u")
    end
    op = operation(u)
    args = arguments(u)
    if op === getindex || (op isa MT.Differential)
        error("unsupported operation in kinetic law: $op")
    end
    if nameof(op) === :ifelse
        return string("<piecewise><piece>", mathml(args[2], ids),
                      mathml(args[1], ids), "</piece><otherwise>",
                      mathml(args[3], ids), "</otherwise></piecewise>")
    end
    name = get(OPS, nameof(op), nothing)
    isnothing(name) && error("unsupported function in kinetic law: $(nameof(op))")
    if name == "root"
        return string("<apply><root/><degree>", num(2), "</degree>",
                      mathml(args[1], ids), "</apply>")
    end
    string("<apply><", name, "/>", join((mathml(a, ids) for a in args)), "</apply>")
end

function writesbml(rs::Catalyst.ReactionSystem, path::AbstractString;
                   modelid = "model", compartment = "c", substance_units = "item",
                   parameter_values = Dict(), initial_amounts = Dict(),
                   constant_species = String[], ratelaw = Catalyst.oderatelaw)
    sys = Catalyst.flatten(rs)
    MT.iscomplete(sys) || (sys = MT.complete(sys))
    sp = MT.unknowns(sys)
    ps = MT.parameters(sys)
    ids = merge(idmap(sp), idmap(ps))

    lookup(tbl, sym) = begin
        u = MT.unwrap(sym)
        for k in (u, sym, Symbol(sym), string(sym))
            haskey(tbl, k) && return float(tbl[k])
        end
        nothing
    end

    value_of(p) = begin
        v = lookup(parameter_values, p)
        isnothing(v) || return v
        MT.hasdefault(p) && return float(Symbolics.value(MT.getdefault(p)))
        error("no value for parameter $p; pass it in parameter_values")
    end
    initial_of(s) = begin
        v = lookup(initial_amounts, s)
        isnothing(v) || return v
        MT.hasdefault(s) ? float(Symbolics.value(MT.getdefault(s))) : 0.0
    end

    io = IOBuffer()
    println(io, "<?xml version=\"1.0\" encoding=\"UTF-8\"?>")
    println(io, "<sbml xmlns=\"", SBML_NS, "\" level=\"3\" version=\"2\">")
    println(io, "<model id=\"", xmlescape(modelid), "\" substanceUnits=\"",
            substance_units, "\" timeUnits=\"dimensionless\" extentUnits=\"",
            substance_units, "\">")

    println(io, "<listOfCompartments>")
    println(io, "<compartment id=\"", compartment,
            "\" spatialDimensions=\"3\" size=\"1\" constant=\"true\"/>")
    println(io, "</listOfCompartments>")

    frozen = Set(constant_species)
    println(io, "<listOfSpecies>")
    for s in sp
        fz = string(s) in frozen || ids[MT.unwrap(s)] in frozen
        println(io, "<species id=\"", ids[MT.unwrap(s)], "\" name=\"",
                xmlescape(string(s)), "\" compartment=\"", compartment,
                "\" initialAmount=\"", initial_of(s),
                "\" hasOnlySubstanceUnits=\"true\" boundaryCondition=\"", fz,
                "\" constant=\"", fz, "\"/>")
    end
    println(io, "</listOfSpecies>")

    if !isempty(ps)
        println(io, "<listOfParameters>")
        for p in ps
            println(io, "<parameter id=\"", ids[MT.unwrap(p)], "\" name=\"",
                    xmlescape(string(p)), "\" value=\"", value_of(p),
                    "\" constant=\"true\"/>")
        end
        println(io, "</listOfParameters>")
    end

    println(io, "<listOfReactions>")
    for (i, rx) in enumerate(Catalyst.reactions(sys))
        println(io, "<reaction id=\"r", i, "\" reversible=\"false\">")
        if !isempty(rx.substrates)
            println(io, "<listOfReactants>")
            for (s, n) in zip(rx.substrates, rx.substoich)
                println(io, "<speciesReference species=\"", ids[MT.unwrap(s)],
                        "\" stoichiometry=\"", n, "\" constant=\"true\"/>")
            end
            println(io, "</listOfReactants>")
        end
        if !isempty(rx.products)
            println(io, "<listOfProducts>")
            for (s, n) in zip(rx.products, rx.prodstoich)
                println(io, "<speciesReference species=\"", ids[MT.unwrap(s)],
                        "\" stoichiometry=\"", n, "\" constant=\"true\"/>")
            end
            println(io, "</listOfProducts>")
        end
        law = ratelaw(rx)
        inrx = Set(MT.unwrap(s) for s in vcat(rx.substrates, rx.products))
        mods = [v for v in Symbolics.get_variables(law)
                if !(v in inrx) && haskey(ids, v) && any(isequal(v), MT.unwrap.(sp))]
        if !isempty(mods)
            println(io, "<listOfModifiers>")
            for m in mods
                println(io, "<modifierSpeciesReference species=\"", ids[m], "\"/>")
            end
            println(io, "</listOfModifiers>")
        end
        println(io, "<kineticLaw><math xmlns=\"", MATHML_NS, "\">")
        println(io, mathml(law, ids))
        println(io, "</math></kineticLaw>")
        println(io, "</reaction>")
    end
    println(io, "</listOfReactions>")
    println(io, "</model>")
    println(io, "</sbml>")

    open(path, "w") do f
        write(f, String(take!(io)))
    end
    path
end

end
