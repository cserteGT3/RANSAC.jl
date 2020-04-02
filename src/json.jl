# export/import as/from JSON

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
