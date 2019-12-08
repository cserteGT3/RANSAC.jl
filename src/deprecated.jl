# lines that are not used but kept for legacy reasons
# what does that even mean???

"""
    impldistance2segment3(point, shape)

Distance: nearest point of the segment is chosen,
then the distance from the linesegments joining to that point is computed,
and the smaller is returned.
"""
function impldistance2segment3(point, shape)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    _, i = nearestpoint(point, shape.contour)
    # if the first is the nearest, then the i-1 segment is the last
    i_1 = i==1 ? lastindex(shape.contour) : i-1
    di_1 = distance2onesegment(point, shape.contour, i_1)
    di = distance2onesegment(point, shape.contour, i)
    # absolute distance is computed
    j = argmin([abs(di_1), abs(di)])

    d = [di_1, di][j]
    return (shape.outwards*d, j)
end

"""
    contournormal(shape, i)

Compute the normal of the i-th segment of a shape.
`shape` should be an `ExtractedTranslational` or `ImplicitTranslational`.
DEPRECATED
"""
function contournormal(shape, i)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    return shape.flipnormal .* segmentnormal(shape.contour, i)
end

"""
    outwardsnormal(shape, i)

Compute the normal of the i-th segment of a shape.
This normal points always outwards.
`shape` should be an `ExtractedTranslational` or `ImplicitTranslational`.
"""
function outwardsnormal(shape, i)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    return shape.outwards .* segmentnormal(shape.contour, i)
end

"""
    distance2onesegment(point, A, i)

Compute distance of `point` to the i-th segment in A.
This distance is computed along the normal vector. No further +- considered.
"""
function distance2onesegment(point, A, i)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    nv = segmentnormal(A, i)
    topoint = point-A[i]
    return dot(topoint, nv)
end

"""
    dist2segment(point, A)

Compute the shortest distance from `point` to the linesegments `A`.
Also return the index of that segment.
This distance is computed along the normal vector. No further +- considered.
"""
function dist2segment(point, A)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    leastd = distance2onesegment(point, A, 1)
    size(A,1) == 1 && return (leastd, 1)
    best = 1
    for i in 2:size(A,1)
        d = distance2onesegment(point, A, i)
        if abs(d) < abs(leastd)
            leastd = d
            best = i
        end
    end
    return (leastd, best)
end

"""
    impldistance2segment(point, shape)

Compute the shortest signed distance from `point` to the linesegments `shape`.
Sign is decided so, that the normal of the surface points outwards.
"""
function impldistance2segment(point, shape)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    d, i = dist2segment(point, shape.contour)
    return (shape.outwards*d, i)
end

"""
    impldistance2segment2(point, shape)

Distance: distance from the nearest point of the shape.contour point.
"""
function impldistance2segment2(point, shape)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    d, i = nearestpoint(point, shape.contour)
    dotp = dot(shape.contour[i]-shape.center, shape.contour[i]-point)
    signi = dotp < 0 ? 1 : -1
    return (signi*d, i)
end



"""
    dandn2contour(point, contour)

Distance and normal from `point` to segment, which is just two points.
Distance and normals signs are consistently computed from the contour.
Also return the number of THE segment.
Inwards-outwards signs must be dealt separately (use `contournormal()` or `outwardsnormal()` with `Ai` ).
"""
function _dandn2segment(point, contour)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    @assert size(contour, 1) > 1 "Contour should contain at least 2 points."
    _, i = nearestpoint(point, contour)
    # is i the first of the contour? if so i-1 is the last
    i_1 = i == firstindex(contour) ? size(contour,1) : i-1
    di_1 = contourdistance(point, contour, i_1)
    di = contourdistance(point, contour, i)
    # selected segment:
    Ai, d = di_1 < di ? (i_1,di_1) : (i,di)
    # normal in the point: normal of the segment
    n = segmentnormal(contour, Ai)
    # sign of the distance: is point-A[1]
    #tp = contour[i]
    tp = midpoint(contour, Ai)
    dsigned = dot(point-tp, n) > 0 ? d : -d
    return (dsigned, n, Ai)
end

function _dn2shape_outwards(point, shape)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    d, n, _ = dandn2segment(point, shape.contour)
    outw = shape.outwards
    return outw*d, outw*n
end

function _dn2shape_contour(point, shape)
    @warn "depricated file: $(@__FILE__), line:$(@__LINE__)"
    d, n, _ = dandn2segment(point, shape.contour)
    flips = shape.flipnormal
    return flips*d, flips*n
end
