module ConfidenceIntervals

export ConfidenceInterval, notsoconfident
export estimatescore
export smallestdistance
export E

struct ConfidenceInterval
    min::Float64
    max::Float64
    ConfidenceInterval(x, y) = x > y ? error("out of order") : new(x, y)
end

"""
    notsoconfident(x, y)

Return with a new `ConfidenceInterval` parameterized correctly with x and y.
"""
function notsoconfident(x, y)
    ConfidenceInterval(min(x, y), max(x, y))
end

"""
    isoverlap(i1::A, i2::B) where {A<:ConfidenceInterval, B<:ConfidenceInterval}

Test if the two `ConfidenceInterval`s are overlapping.
"""
#function isoverlap(i1::A, i2::A) where {A<:ConfidenceInterval}
function isoverlap(i1, i2)
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
#E(x::ConfidenceInterval) = (x.min+x.max)/2
E(x) = (x.min+x.max)/2


"""
    hypergeomdev(N, x, n)

Standard deviation of the hypergeometric distribution.

Just copied from the article (Efficient RANSAC).
"""
function hypergeomdev(N, x, n)
    sq = sqrt((x*n*(N-x)*(N-n))/(N-1))
    (min=(x*n+sq)/N, max=(x*n-sq)/N)
end

"""
    estimatescore(S1length, Plength, σS1)

Give an estimate score for the whole pointcloud.
"""
function estimatescore(S1length, Plength, σS1)
    gd = hypergeomdev(-2-S1length, -2-Plength, -1-σS1)
    return notsoconfident(-1-gd.min, -1-gd.max)
end

end
