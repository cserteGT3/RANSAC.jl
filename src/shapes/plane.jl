# fitting

struct FittedPlane{A<:AbstractArray} <: FittedShape
    isplane::Bool
    point::A
    normal::A
end

function isshape(shape::FittedPlane)
    return shape.isplane
end

"""
    fitplane(p, n, params)

Fit a plane to 3 points. Their and additional point's normals are used to validate the fit.

A collinearity check is used to not filter out points on one line.
# Arguments:
- `alpharad::Real`: maximum difference between the normals (in radians).
- `collin_threshold::Real=0.2`: threshold for the collinearity check (lower is stricter).
"""
function fitplane(p, n, params)
    @unpack α_plane, collin_threshold = params
    lp = length(p)
    @assert length(p) > 2 "At least 3 point is needed."
    @assert lp == length(n) "Size must be the same."
    crossv = normalize(cross(p[2]-p[1], p[3]-p[1]))
    # how many point's normal must be checked
    if norm(crossv) < collin_threshold
        # The points are collinear
        # solution: use 1 more point
        if length(p) < 4
            # return with false, cause no more points can be used
            return FittedPlane(false, NaNVec, NaNVec)
        end
        crossv = normalize(cross(p[2]-p[4], p[3]-p[1]))
        if norm(crossv) < collin_threshold
            # return false, cause these are definitely on one line
            return FittedPlane(false, NaNVec, NaNVec)
        end
    end
    # here we have the normal of the theoretical plane
    norm_ok = falses(lp)
    invnorm_ok = falses(lp)

    thr = cos(α_plane)
    for i in 1:lp
        dotp = dot(crossv, normalize(n[i]))
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    norm_ok == trues(lp) && return FittedPlane(true, p[1], crossv)
    invnorm_ok == trues(lp) && return FittedPlane(true, p[1], -1*crossv)
    return FittedPlane(false, NaNVec, NaNVec)
end

# bitmapping

function scorecandidate(pc, candidate::ShapeCandidate{T}, subsetID, params) where {T<:FittedPlane}
    ps = @view pc.vertices[pc.subsets[subsetID]]
    ns = @view pc.normals[pc.subsets[subsetID]]
    ens = @view pc.isenabled[pc.subsets[subsetID]]

    cp, pp = compatiblesPlane(candidate.shape, ps, ns, params)
    inder = cp.&ens
    inpoints = (pc.subsets[subsetID])[inder]
    score = estimatescore(length(pc.subsets[subsetID]), pc.size, length(inpoints))
    pc.levelscore[candidate.octree_lev] += E(score)
    return ScoredShape(candidate, score, inpoints)
end

"""
    project2plane(plane, points)

Project `points` on to the `plane`.

# Arguments:
- 'plane::FittedPlane': a plane.
- 'points::AbstractArray': an array of points (array-of-arrays).
"""
function project2plane(plane, points)
    # create a coordinate frame
    # z is the plane's normal
    o_z = normalize(plane.normal)
    # x is a random orthogonal vector (in the plane)
    o_x = normalize(arbitrary_orthogonal(o_z))
    # y is created so it's a right hand coord. frame
    o_y = normalize(cross(o_z, o_x))
    answer = similar(points)
    # get the coordinates of the points in the prev. created coord. frame
    for i in eachindex(points)
        v = points[i]-plane.point
        answer[i] = eltype(answer)(dot(o_x,v), dot(o_y,v), dot(o_z,v))
    end
    answer
end

"""
    compatiblesPlane(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the plane.
Give back the projected points too for parameter space magic.

Compatibility is measured with an `eps` distance to the plane and an `alpharad` angle to it's normal.
"""
function compatiblesPlane(plane, points, normals, params)
    @unpack ϵ_plane, α_plane = params
    @assert length(points) == length(normals) "Size must be the same."
    projecteds = project2plane(plane, points)
    # eps check
    c1 = [abs(a[3]) < ϵ_plane for a in projecteds]
    # alpha check
    c2 = [isparallel(plane.normal, normals[i], α_plane) && c1[i] for i in eachindex(normals)]
    # projecteds[c2] are the compatible points
    return c2, projecteds
end

"""
    refit(s, pc, ϵ, α)

Refit plane. Only s.inpoints is updated.
"""
function refitplane(s, pc, params)
    # TODO: use octree for that
    cp, _ = compatiblesPlane(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], params)
    s.inpoints = ((1:pc.size)[pc.isenabled])[cp]
    s
end
