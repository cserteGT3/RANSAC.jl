## export to JSON

"""
    toDict(s::FittedShape)

Convert `s` to a `Dict{String,Any}`. It's "type" is defined by [`strt`](@ref).
All fields of the struct is saved to the dict.
"""
function toDict(s::FittedShape)
    d = Dict{String,Any}("type"=>strt(s))
    for f in fieldnames(typeof(s))
        sf = String(f)
        merge!(d, Dict(sf=>getfield(s, f)))
    end
    d
end

toDict(sc::ShapeCandidate) = toDict(sc.shape)

toDict(scs::ScoredShape) = toDict(scs.candidate.shape)

"""
    toDict(a::Vector{T}) where {T<:Union{FittedShape,ShapeCandidate,ScoredShape}}

Convert a vector of shapes to a `Dict{String,Any}`. The top key is a "primitive", whose value is the array of the shapes.
See the documentation for examples.
"""
function toDict(a::Vector{T}) where {T<:Union{FittedShape,ShapeCandidate,ScoredShape}}
    ad = toDict.(a)
    return Dict("primitives"=>ad)
end

"""
    printJSON(io::IO, s, indent)

Print a `FittedShape`, `ShapeCandidate`, `ScoredShape` or a vector of them to `io` as a JSON string.
With `indent` given, it prints a representation with newlines and indents.

# Arguments:
- `io::IO`: must be specified, use `stdout` for interactive purposes.
- `s`: a `FittedShape`, `ShapeCandidate`, `ScoredShape` or a vector of one of them.
- `indent::Int`: indentation level.
"""
function exportJSON(io::IO, s, indent)
    JSON.print(io, toDict(s), indent)
end

"""
    printJSON(io::IO, s)

Print a `FittedShape`, `ShapeCandidate`, `ScoredShape` or a vector of them to `io` as a compact JSON string.

# Arguments:
- `io::IO`: must be specified, use `stdout` for interactive purposes.
- `s`: a `FittedShape`, `ShapeCandidate`, `ScoredShape` or a vector of one of them.
"""
function exportJSON(io::IO, s)
    JSON.print(io, toDict(s))
end


## import yaml config file

"""
    readconfig(fname, paramtype=RANSACParameters)

Read a config file to a parameter sruct type of `paramtype` (`RANSACParameters` is the default).

# Arguments:
- `fname`: file name.
- `paramtype=RANSACParameters`: type of the parameter struct. One must be able to instantiate it with keyword arguments only (eg. `paramtype()` should work).
"""
function readconfig(fname, paramtype=RANSACParameters)
    fdict = YAML.load(open(fname))
    # get element type; defaults to Float64
    #ftype = get(fdict, "eltype", "Float64") == "Float32" ? Float32 : Float64

    # field names of paramtype
    field_names = fieldnames(paramtype)
    vals = []
    def_struct = paramtype()

    for f_n in field_names
        # `shape_types` is an array of symbols, must be handled differently
        if f_n === :shape_types
            if ! (get(fdict, "shape_types", nothing) === nothing)
                push!(vals, Symbol.(fdict["shape_types"]))
                continue
            end
        end
        # default value
        def_val = getproperty(def_struct, f_n)
        # fdict uses string keys
        param_val = get(fdict, String(f_n), def_val)
        push!(vals, param_val)
    end
    return paramtype(;zip(field_names, Tuple(vals))...)
end