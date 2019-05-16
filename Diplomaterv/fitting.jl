module Fitting

include("utilities.jl")
include("ConfidenceIntervals.jl")

using StaticArrays: SVector, MVector
using LinearAlgebra: cross, ×, dot, normalize, normalize!, norm, det
using ZChop: zchop, zchop!

using .Utilities
using .ConfidenceIntervals: ConfidenceInterval, E, isoverlap

export FittedShape, isshape
export FittedPlane, isplane
export FittedSphere, issphere
export FittedCylinder, iscylinder
export ShapeCandidate, findhighestscore
export ScoredShape
export largestshape
#export refit

"""
An abstract type that wraps the fitted shapes.
"""
abstract type FittedShape end

struct FittedPlane{A<:AbstractArray} <: FittedShape
    isplane::Bool
    point::A
    normal::A
end

function isshape(shape::FittedPlane)
    return shape.isplane
end

mutable struct ShapeCandidate{S<:FittedShape}
    shape::S
    octree_lev::Int
end

# mutable struct ScoredShape{S<:FittedShape, A<:AbstractArray}

mutable struct ScoredShape{A<:AbstractArray}
    candidate::ShapeCandidate
    score
    inpoints::A
end

# TODO: delete these
# ShapeCandidate(shape, score, octlev) = ShapeCandidate(shape, score, [], false, octlev)
# ShapeCandidate(shape, octlev) = ShapeCandidate(shape, ConfidenceInterval(0,0), [], false, octlev)

"""
    findhighestscore(A)

Find the largest expected value in an array of `ScoredShape`s.

Indicate if there's an overlap.
"""
function findhighestscore(A)
    (length(A) > 0) || return (index = 0, overlap = false)
    ind = 1
    highest = E(A[1].score)
    for i in eachindex(A)
        esc = E(A[i].score)
        if esc > highest
            highest = esc
            ind = i
        end
    end
    overlaps = [isoverlap(A[ind].score, a.score) for a in A]
    overlaps[ind] = false
    # no overlap:
    overlaps == falses(length(A)) && return (index = ind, overlap = false)
    # overlap:
    return (index = ind, overlap = true)
end

"""
    largestshape(A)

Find the largest shape in an array of `ScoredShape`s.
"""
function largestshape(A)
    length(A) < 1 && return (index = 0, size = 0)
    bestscore = length(A[1].inpoints)
    bestind = 1
    for i in eachindex(A)
        length(A[i].inpoints) <= bestscore && continue
        bestscore = length(A[i].inpoints)
        bestind = i
    end
    return (index = bestind, size = bestscore)
end

"""
This is "dummy" vector for type stability.
"""
const NaNVec = SVector(0.0,0.0,0.0)

"""
    isplane(p, n, alpharad, collin_threshold = 0.2)

Fit a plane to 3 points. Their and additional point's normals are used to validate the fit.

A collinearity check is used to not filter out points on one line.
# Arguments:
- `alpharad::Real`: maximum difference between the normals (in radians).
- `collin_threshold::Real=0.2`: threshold for the collinearity check (lower is stricter).
"""
function isplane(p, n, alpharad, collin_threshold = 0.2)
    lp = length(p)
    @assert length(p) > 2 "At least 3 point is needed."
    @assert lp == length(n) "Size must be the same."
    crossv = MVector{3}(normalize(cross(p[2]-p[1], p[3]-p[1])))
    zchop!(crossv)
    # how many point's normal must be checked
    if norm(crossv) < collin_threshold
        # The points are collinear
        # solution: use 1 more point
        if length(p) < 4
            # return with false, cause no more points can be used
            return FittedPlane(false, NaNVec, NaNVec)
        end
        crossv = MVector{3}(normalize(cross(p[2]-p[4], p[3]-p[1])))
        if norm(crossv) < collin_threshold
            # return false, cause these are definitely on one line
            return FittedPlane(false, NaNVec, NaNVec)
        end
    end
    # here we have the normal of the theoretical plane
    norm_ok = falses(lp)
    invnorm_ok = falses(lp)

    thr = cos(alpharad)
    for i in 1:lp
        dotp = dot(crossv, zchop!(MVector{3}(normalize(n[i]))))
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    norm_ok == trues(lp) && return FittedPlane(true, p[1], SVector{3}(crossv))
    invnorm_ok == trues(lp) && return FittedPlane(true, p[1], SVector{3}(-1*crossv))
    return FittedPlane(false, NaNVec, NaNVec)
end

struct FittedSphere{A<:AbstractArray, R<:Real} <: FittedShape
    issphere::Bool
    center::A
    radius::R
    outwards::Bool
end

function isshape(shape::FittedSphere)
    return shape.issphere
end

function fitsphere(v, n)
    n1n = normalize(n[1])
    n2n = normalize(n[2])

    b = false
    if abs(dot(n1n, n2n)) > 0.98
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

        if nk < 0.02 || nh < 0.02
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
                return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]), b)
            else
                # point in different direction -> -
                M = v[1] - nh/nk * n1n
                return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]), b)
            end
        end
    end
end

function setsphereOuterity(sp, b)
    FittedSphere(sp.issphere, sp.center, sp.radius, b)
end

"""
    issphere(p, n, epsilon, alpharad)

Fit a sphere to 2 points. Additional points and their normals are used to validate the fit.

# Arguments:
- `epsilon::Real`: maximum distance-difference between the fitted and measured spheres.
- `alpharad::Real`: maximum difference between the normals (in radians).
"""
function issphere(p, n, epsilon, alpharad)
    pl = length(p)
    @assert pl == length(n) "Size must be the same."
    @assert pl > 2 "Size must be at least 3."
    # "forcefit" a sphere
    sp = fitsphere(p, n)
    # check if real sphere
    sp.issphere || return sp
    thr = cos(alpharad)
    vert_ok = falses(pl)
    norm_ok = falses(pl)
    invnorm_ok = falses(pl)
    for i in 1:pl
        # vertice check
        vert_ok[i] = abs(norm(p[i]-sp.center)-sp.radius) < epsilon
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

#=
"""
    refit(s, pc, ϵ, α)

Refit plane. Only s.inpoints is updated.
"""
function refit(s, pc, ϵ, α)
    # TODO: use octree for that
    cp, _ = compatiblesPlane(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
    s.inpoints = ((1:pc.size)[pc.isenabled])[cp]
    s
end

"""
    refit(s, pc, ϵ, α)

Refit sphere. Only s.inpoints is updated.
"""
function refit(s, pc, ϵ, α)
    # TODO: use octree for that
    cpl, uo, sp = compatiblesSphere(s.candidate.shape, pc.vertices[pc.isenabled], pc.normals[pc.isenabled], ϵ, α)
    # verti: összes pont indexe, ami enabled és kompatibilis
    verti = (1:pc.size)[pc.isenabled]
    underEn = uo.under .& cpl
    overEn = uo.over .& cpl
    s.inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
    s
end
=#

struct FittedCylinder{A<:AbstractArray, R<:Real} <: FittedShape
    iscylinder::Bool
    axis::A
    center::A
    radius::R
    outwards::Bool
end

function isshape(shape::FittedCylinder)
    return shape.iscylinder
end

function setcylinderOuterity(fc, b)
    FittedCylinder(fc.iscylinder, fc.axis, fc.center, fc.radius, b)
end

function fitcylinder(p, n, epsilon, alpharad; parallel_threshold_deg = 1)
    # the normals are too parallel so cannot fit cylinder to it
    if abs(dot(n[1], n[2])) > cos(deg2rad(parallel_threshold_deg))
        return FittedCylinder(false, NaNVec, NaNVec, 0, false)
    end

    an = normalize(n[1] × n[2])

    function lineplaneintersect(n, u, w)
        # http://geomalgorithms.com/a05-_intersect-1.html
        # n=plane normal
        # u= direction of the line
        # w=

        return dot(-n, w)/dot(n,u)
    end

    function project2plane(n, w)
        #n: normal
        #w: helyvektor
        return w+n*lineplaneintersect(n,n,w)
    end

    function projectto2d(xaxis, yaxis, zaxis, p1)
        # xax: pl a ponban muataó vektor
        xx = xaxis[1];
        xy = xaxis[2];
        xz = xaxis[3];

        yx = yaxis[1];
        yy = yaxis[2];
        yz = yaxis[3];

        zx = zaxis[1];
        zy = zaxis[2];
        zz = zaxis[3];

        px = p1[1];
        py = p1[2];
        pz = p1[3];


        r1 = -((-(pz*yy*zx) + py*yz*zx + pz*yx*zy - px*yz*zy - py*yx*zz + px*yy*zz)/(xz*yy*zx - xy*yz*zx - xz*yx*zy + xx*yz*zy + xy*yx*zz - xx*yy*zz))
        r2 = -((pz*xy*zx - py*xz*zx - pz*xx*zy + px*xz*zy + py*xx*zz - px*xy*zz)/(xz*yy*zx - xy*yz*zx - xz*yx*zy + xx*yz*zy + xy*yx*zz - xx*yy*zz))
        r3 = -((pz*xy*yx - py*xz*yx - pz*xx*yy + px*xz*yy + py*xx*yz - px*xy*yz)/(-(xz*yy*zx) + xy*yz*zx + xz*yx*zy - xx*yz*zy - xy*yx*zz + xx*yy*zz))

        return SVector(r1, r2)
    end

    function lineintersectionpoint(l1, l2)
        # l1 = [[p1x, p1y],[p2x,p2y]]
        a = l1[1]
        b = l1[2]
        c = l2[1]
        d = l2[2]

        amb_ = a-b
        cmd_ = c-d

        d1 = det(vcat(a',b'))
        d2 = det(vcat(c',d'))
        d3 = det(vcat(amb_',cmd_'))
        return (d1*cmd_-d2*amb_)/d3
    end

    xax = normalize(project2plane(an, p[1]))
    yax = normalize(an × xax)

    p11proj = projectto2d(xax, yax, an, project2plane(an, p[1]))
    p12proj = projectto2d(xax, yax, an, project2plane(an, p[1]+n[1]))

    p21proj = projectto2d(xax, yax, an, project2plane(an, p[2]))
    p22proj = projectto2d(xax, yax, an, project2plane(an, p[2]+n[2]))

    interc = lineintersectionpoint([p11proj,p12proj], [p21proj,p22proj])

    c = interc[1]*xax+interc[2]*yax

    nnormies = [norm(pt - c - an*dot( an, pt-c )) for pt in p[1:2]]
    R = (nnormies[1] + nnormies[2])/2
    # R = (norm(interc-p11proj) + norm(interc-p21proj))/2

    #p[1] - an*dot( an, p[1]-c )

    outw = dot(p12proj-p11proj, p11proj-interc) > 0

    return FittedCylinder(true, an, c, R, outw)


#=
    p1proj = p[1] - an*dot(p[1], an)
    p2proj = p[2] - an*dot(p[2], an)

    v = [p1proj, p2proj]

    # first check if they intersect
    g = v[2]-v[1]
    h = cross(n[2], g)
    k = cross(n[2], n[1])
    nk = norm(k)
    nh = norm(h)

    if abs(nk) < 0.02 || abs(nh) < 0.02
        error("That is impossible with p: $p, n: $n and nk=$nk, nh=$nh")
    else
        # intersection
        if dot(h, k) > 0
            # point the same direction -> +
            c = v[1] + nh/nk * n[1]
            #return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]), b)
        else
            # point in different direction -> -
            c = v[1] - nh/nk * n[1]
            #return FittedSphere(true, SVector{3}(zchop!(MVector{3}(M))), norm(M-v[1]), b)
        end



        # radius is the average of: center-point
        R = (norm(p1proj-c) + norm(p2proj-c))/2
        # outwards?
        outw = isparallel(n[1], p1proj-c, alpharad)

        return FittedCylinder(true, an, c, R, outw)

    end=#
end

"""
    iscylinder(p, n, epsilon, alpharad)

Fit a cylinder to 2 points. Additional points and their normals are used to validate the fit.
Normals are expected to be normalized.

# Arguments:
- `epsilon::Real`: maximum distance-difference between the fitted and measured spheres.
- `alpharad::Real`: maximum difference between the normals (in radians).
"""
function iscylinder(p, n, epsilon, alpharad)
    pl = length(p)
    @assert pl == length(n) "Size must be the same."
    @assert pl > 2 "Size must be at least 3."

    # "forcefit" a cylinder
    fc = fitcylinder(p, n, epsilon, alpharad, parallel_threshold_deg=3)

    fc.iscylinder || return fc


    thr = cos(alpharad)
    vert_ok = falses(pl)
    norm_ok = falses(pl)
    invnorm_ok = falses(pl)
    for i in 1:pl
        # current normal
        curr_norm = p[i] - fc.axis*dot( fc.axis, p[i]-fc.center ) - fc.center
        # vertice check
        vert_ok[i] = abs(norm(curr_norm)-fc.radius) < epsilon
        # normal check
        dotp = dot( normalize(curr_norm), n[i] )
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end
    vert_ok == trues(pl) || return FittedCylinder(false, NaNVec, NaNVec, 0, false)
    norm_ok == trues(pl) && return setcylinderOuterity(fc, true)
    invnorm_ok == trues(pl) && return setcylinderOuterity(fc, false)
    return FittedCylinder(false, NaNVec, NaNVec, 0, false)
end


end #module
