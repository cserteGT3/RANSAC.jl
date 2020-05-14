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

toDict(sc::ExtractedShape) = toDict(sc.shape)

"""
    toDict(a::Vector{T}) where {T<:Union{FittedShape,ExtractedShape}}

Convert a vector of shapes to a `Dict{String,Any}`.
The top key is a "primitive", whose value is the array of the shapes.
See the documentation for examples.
"""
function toDict(a::Vector{T}) where {T<:Union{FittedShape,ExtractedShape}}
    ad = toDict.(a)
    return Dict("primitives"=>ad)
end

"""
    printJSON(io::IO, s, indent)

Print a `FittedShape`, `ExtractedShape` or a vector of them to `io` as a JSON string.
With `indent` given, it prints a representation with newlines and indents.

# Arguments:
- `io::IO`: must be specified, use `stdout` for interactive purposes.
- `s`: a `FittedShape`, `ExtractedShape` or a vector of one of them.
- `indent::Int`: indentation level.
"""
function exportJSON(io::IO, s, indent)
    JSON.print(io, toDict(s), indent)
end

"""
    printJSON(io::IO, s)

Print a `FittedShape`, `ExtractedShape`
or a vector of them to `io` as a compact JSON string.

# Arguments:
- `io::IO`: must be specified, use `stdout` for interactive purposes.
- `s`: a `FittedShape`, `ExtractedShape` or a vector of one of them.
"""
function exportJSON(io::IO, s)
    JSON.print(io, toDict(s))
end


## import yaml config file

"""
    readconfig(fname; toextend=DEFAULT_PARAMETERS, shapedict=DEFAULT_SHAPE_DICT)

Read a config file to a `NamedTuple`.
A "base" ntuple is expected, that gets overwritten/extended
with the values in the config file.
A `Dict{String,FittedShape}` dictionary is also expected,
that translates the string primitive types to julia types.

# Arguments

- `fname`: name of the config file.
- `toextend=DEFAULT_PARAMETERS`: a named tuple,
    that will be overwritten/extended with the values from the config file.
- `shapedict=DEFAULT_SHAPE_DICT`: a dictionary that
    translates the string primitive names to julia types.
"""
function readconfig(fname; toextend=DEFAULT_PARAMETERS, shapedict=DEFAULT_SHAPE_DICT)
    fio = open(fname, "r")
    fdict = YAML.load(fio)
    close(fio)

    newnd = ((;Symbol(k)=>dict2nt(v)) for (k,v) in pairs(fdict))
    for v in newnd
        if haskey(v[1], :shape_types)
            nv = merge(v[1], (shape_types=[shapedict[k] for k in v[1].shape_types],))
            toextend = ransacparameters(toextend; iteration=nv)
        else
            toextend = ransacparameters(toextend; v...)
        end
    end
    return toextend
end

function dict2nt(dict)
    dd = [NamedTuple{Tuple(Symbol.(keys(d)))}(values(d)) for d in dict]
    if length(dd) < 2
        return dd[1]
    end
    newd = dd[1]
    for i in eachindex(dd)
        i == 1 && continue
        newd = merge(newd, dd[i])
    end
    return newd
end
