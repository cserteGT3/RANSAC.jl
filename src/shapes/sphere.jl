# fitting

struct FittedSphere{A<:AbstractArray, R<:Real} <: FittedShape
    issphere::Bool
    center::A
    radius::R
    outwards::Bool
end

function isshape(shape::FittedSphere)
    return shape.issphere
end

function fit2pointsphere(v, n, params)
    @unpack parallelthrdeg, sphere_par = params
    n1n = normalize(n[1])
    n2n = normalize(n[2])

    b = false
    if abs(dot(n1n, n2n)) > cosd(parallelthrdeg)
        # parallel normals
        # fit it, if bad, will fall out later
        centerp = (v[1]+v[2])/2
        return FittedSphere(true, centerp, norm(centerp-v[1]), b)
    else
        # not parallel normals
        # first check if they intersect
        g = v[2]-v[1]
        h = cross(n2n, g)
        k = cross(n2n, n1n)
        nk = norm(k)
        nh = norm(h)

        if nk < sphere_par || nh < sphere_par
            # if nk == 0 || nh == 0, but it's numeric
            # no intersection
            n2 = n2n × (n1n × n2n)
            n1 = n1n × (n2n × n1n)
            c1 = v[1] + dot( (v[2]-v[1]), n2 )/dot( n[1], n2 ) * n[1]
            c2 = v[2] + dot( (v[1]-v[2]), n1 )/dot( n[2], n1 ) * n[2]
            centerp = (c1+c2)/2
            r = (norm(v[1]-centerp) + norm(v[1]-centerp))/2
            return FittedSphere(true, centerp, r, b)
        else
            # intersection
            if dot(h, k) > 0
                # point the same direction -> +
                M = v[1] + nh/nk * n1n
                return FittedSphere(true, M, norm(M-v[1]), b)
            else
                # point in different direction -> -
                M = v[1] - nh/nk * n1n
                return FittedSphere(true, M, norm(M-v[1]), b)
            end
        end
    end
end

function setsphereOuterity(sp, b)
    FittedSphere(sp.issphere, sp.center, sp.radius, b)
end

"""
    fitsphere(p, n, params)

Fit a sphere to 2 points. Additional points and their normals are used to validate the fit.

# Arguments:
- `epsilon::Real`: maximum distance-difference between the fitted and measured spheres.
- `alpharad::Real`: maximum difference between the normals (in radians).
"""
function fitsphere(p, n, params)
    @unpack ϵ_sphere, α_sphere = params
    pl = length(p)
    @assert pl == length(n) "Size must be the same."
    @assert pl > 2 "Size must be at least 3."
    # "forcefit" a sphere
    sp = fit2pointsphere(p, n, params)
    # check if real sphere
    sp.issphere || return sp
    thr = cos(α_sphere)
    vert_ok = falses(pl)
    norm_ok = falses(pl)
    invnorm_ok = falses(pl)
    for i in 1:pl
        # vertice check
        vert_ok[i] = abs(norm(p[i]-sp.center)-sp.radius) < ϵ_sphere
        # normal check
        dotp = dot( normalize(p[i]-sp.center), normalize(n[i]) )
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    vert_ok == trues(pl) || return FittedSphere(false, NaNVec, 0, false)
    norm_ok == trues(pl) && return setsphereOuterity(sp, true)
    invnorm_ok == trues(pl) && return setsphereOuterity(sp, false)
    return FittedSphere(false, NaNVec, 0, false)
end

# bitmapping

function scorecandidate(pc, candidate::ShapeCandidate{T}, subsetID, params) where {T<:FittedSphere}
    ps = @view pc.vertices[pc.subsets[subsetID]]
    ns = @view pc.normals[pc.subsets[subsetID]]
    ens = @view pc.isenabled[pc.subsets[subsetID]]

    cpl, uo, sp = compatiblesSphere(candidate.shape, ps, ns, params)
    # verti: összes pont indexe, ami enabled és kompatibilis
    # lenne, ha működne, de inkább a boolean indexelést machináljuk
    verti = pc.subsets[subsetID]
    underEn = uo.under .& cpl
    overEn = uo.over .& cpl

    inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
    score = estimatescore(length(pc.subsets[subsetID]), pc.size, length(inpoints))
    pc.levelscore[candidate.octree_lev] += E(score)
    return ScoredShape(candidate, score, inpoints)
end

"""
    compatiblesSphere(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the sphere.
Give back the projected points too for parameter space magic.
Return a bool indexer for (under,over) too.

Compatibility is measured with an `eps` distance to the sphere and an `alpharad` angle to it's normal.
"""
function compatiblesSphere(sphere, points, normals, params)
    @unpack ϵ_sphere, α_sphere = params
    @assert length(points) == length(normals) "Size must be the same."
    # eps check
    o = sphere.center
    R = sphere.radius
    c1 = [abs(norm(a-o)-R) < ϵ_sphere for a in points]
    # alpha check
    if sphere.outwards
        c2 = [isparallel(normalize(points[i]-o), normals[i], α_sphere) && c1[i] for i in eachindex(points)]
    else
        c2 = [isparallel(normalize(o-points[i]), normals[i], α_sphere) && c1[i] for i in eachindex(points)]
    end

    under = falses(length(points))
    over = falses(length(points))
    for i in eachindex(points)
        if points[i][3] <= o[3]
            under[i] = true
        else
            over[i] = true
        end
    end
    # TODO
    # normalize by it's own length or maximum length??? r+3ϵ
    param_points = [unitdisk2square(normalize(a[1:2]-o[1:2])) for a in points]
    return c2, (under=under, over=over), param_points
end

"""
    refit(s, pc, ϵ, α)

Refit sphere. Only s.inpoints is updated.
"""
function refitsphere(s, pc, params)
    # TODO: use octree for that
    cpl, uo, sp = compatiblesSphere(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], params)
    # verti: összes pont indexe, ami enabled és kompatibilis
    verti = (1:pc.size)[pc.isenabled]
    underEn = uo.under .& cpl
    overEn = uo.over .& cpl
    s.inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
    s
end
