"""
isntequal(nt1::T1, nt2::T2) where {T1<:NamedTuple, T2<:NamedTuple}

Kind of `isequal()` for `NamedTuple`. Order of the keys doesn't matter.
"""
function isntequal(nt1::T1, nt2::T2) where {T1<:NamedTuple, T2<:NamedTuple}

# Collect the keys of a `NamedTuple` into an array, then sort it.
sortkeys(nt) = nt |> keys |> collect |> sort

sfn1 = sortkeys(nt1)
sfn2 = sortkeys(nt2)
fnameseq = sfn1 == sfn2

# keys are not equal -> they are not equal
fnameseq || return false
for key in sfn1
    v1 = getindex(nt1, key)
    v2 = getindex(nt2, key)
    # recursively call 
    isntequal(v1, v2) || return false
end
return true
end

"""
    isntequal(a, b)

Define for any type, using default `isequal`.
It is used when a value is not `NamedTuple` in recursive calls.
"""
isntequal(a, b) = a == b
