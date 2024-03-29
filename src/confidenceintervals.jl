struct ConfidenceInterval
    min::Float64
    max::Float64
    E::Float64
    ConfidenceInterval(x, y) = x > y ? error("out of order") : new(x, y, (x+y)/2)
end

Base.show(io::IO, x::ConfidenceInterval) =
    print(io, "CI: [$(x.min), $(x.max)]")

Base.show(io::IO, ::MIME"text/plain", x::ConfidenceInterval) =
    print(io, "ConfidenceInterval: [$(x.min), $(x.max)]")


"""
    notsoconfident(x, y)

Return with a new `ConfidenceInterval` parameterized correctly with x and y.
"""
function notsoconfident(x, y)
    ConfidenceInterval(min(x, y), max(x, y))
end

"""
    isoverlap(i1::A, i2::A) where {A<:ConfidenceInterval}

Test if the two [ConfidenceInterval](@ref)s are overlapping.
"""
function isoverlap(i1::A, i2::A) where {A<:ConfidenceInterval}
    i1.min == i2.min && return true
    if i1.min < i2.min
        return i2.min <= i1.max
    else
        return isoverlap(i2, i1)
    end
end

"""
    E(x::ConfidenceInterval)

Expected value of a [ConfidenceInterval](@ref).
"""
E(x::ConfidenceInterval) = x.E


"""
    hypergeomdev(N, x, n)

Standard deviation of the hypergeometric distribution.

Just copied from the article (Efficient RANSAC).
"""
function hypergeomdev(N, x, n)
    sq_ = (x*n*(N-x)*(N-n))/(N-1)
    # handle floating point errors:
    # small negative numbers, that should be zero
    sq = sq_ < 0 ? zero(sq_) : sqrt(sq_)
    (min=(x*n+sq)/N, max=(x*n-sq)/N)
end

"""
    estimatescore(S1length, Plength, σS1)

Give an estimate score for the whole pointcloud.

# Arguments
- `S1length`: number of points that are searched for compatible points (size of a subset).
- `Plength`: size of the whole point cloud.
- `σS1`: number of compatible points.
"""
function estimatescore(S1length, Plength, σS1)
    gd = hypergeomdev(-2-S1length, -2-Plength, -1-σS1)
    return notsoconfident(-1-gd.min, -1-gd.max)
end
