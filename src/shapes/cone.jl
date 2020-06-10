# fitting

"""
    struct FittedCone{A<:AbstractArray, R<:Real} <: FittedShape

Cone primitive, defined by its apex,
its axis (that points from the apex towards the opening),
and its opening angle in radians.
Also stored, if the normals point outwards of the shape.
"""
struct FittedCone{A<:AbstractArray, R<:Real} <: FittedShape
    apex::A
    # points from apex towards the opening
    # unit length
    axis::A
    # opening angle in radians
    opang::R
    outwards::Bool
end

Base.show(io::IO, x::FittedCone) =
    print(io, """cone, ω: $(x.opang)""")

Base.show(io::IO, ::MIME"text/plain", x::FittedCone{A, R}) where {A, R} =
    print(io, """FittedCone{$A, $R}\napex: $(x.apex), axis: $(x.axis), ω: $(x.opang), $(x.outwards ? "outwards" : "inwards" )""")

strt(x::FittedCone) = "cone"

function defaultshapeparameters(::Type{FittedCone})
    return (cone=(ϵ=0.3, α=deg2rad(5), minconeopang=deg2rad(2)),)
end

function setconeOuterity(fc, b)
    FittedCone(fc.apex, fc.axis, fc.opang, b)
end

## fitting

function fit3pointcone(p, n)
    # rank of the coefficient matrix
    ftype = eltype(eltype(p))
    r = Array{ftype,2}(undef, (3,3))
    for i in 1:3; for j in 1:3; r[i,j] = n[i][j]; end; end
    rank(r) == 3 || return nothing
    ds = [dot(p[i], n[i]) for i in 1:3]
    # rank of the augmented matrix
    rv = hcat(r, -1 .* ds)
    rank(rv) == 3 || return nothing
    # apex
    ap = SVector{3, ftype}(r\ds)
    # axis
    axis3p = [ap+((p[i]-ap)/norm(p[i]-ap)) for i in 1:3]
    ax = normalize(cross(axis3p[2]-axis3p[1], axis3p[3]-axis3p[1]))
    midp = sum(axis3p)/3
    dirv = normalize(midp-ap)
    ax = dot(ax, dirv) < 0 ? -1*ax : ax
    # opening angle
    angles = [acos(clamp(dot(normalize(p[i]-ap), ax), -1, 1)) for i in 1:3]
    opangle = 2*sum(angles)/3
    return FittedCone(ap, ax, opangle, true)
end

"""
    project2cone(cone, p)

Return the distance from the cone's surface and normal of the surface at that point.
"""
function project2cone(cone, p)
    # vector to the point from the apex
    to_point = cone.apex-p
    to_pointn = normalize(to_point)
    # rotation axis defined by the above vector and the axis of the cone
    rot_ax = normalize(cross(cone.axis, to_pointn))
    # normal: from the rotation axis and the axis fo the cone
    # need to be rotated to get the normal of the cone
    comp_n = normalize(cross(cone.axis, rot_ax))
    # rotation matrix from the rotation axis and opening angle
    rM = rodriguesrad(rot_ax, -cone.opang/2)
    # normal of the cone at the point
    current_normal = normalize(rM*comp_n)
    # distance from the surface of the cone to the point
    dist = dot(-current_normal, -to_point)
    # return the distance and the normal
    return (dist, current_normal)
end

function validatecone(cone, ps, ns, params)
    #@unpack α_cone, ϵ_cone, minconeopang = params
    @extract params : params_cone=cone
    @extract params_cone : α_cone=α ϵ_cone=ϵ minconeopang
    calcs = [project2cone(cone, ps[i]) for i in eachindex(ps)]
    for i in eachindex(calcs)
        if calcs[i][1] > ϵ_cone
            return nothing
        end
    end

    #cone with small opening angle is filtered
    cone.opang < minconeopang && return nothing

    lp = size(ps, 1)
    norm_ok = falses(lp)
    invnorm_ok = falses(lp)

    thr = cos(α_cone)
    for i in eachindex(calcs)
        dotp = dot(calcs[i][2], ns[i])
        norm_ok[i] = dotp > thr
        invnorm_ok[i] = dotp < -thr
    end

    norm_ok == trues(lp) && return setconeOuterity(cone, true)
    invnorm_ok == trues(lp) && return setconeOuterity(cone, false)
    return nothing
end

"""
    fit(::Type{FittedCone}, p, n, pc, params)

Fit a cone to 3 points. Normals are expected to be normalized.
Return `nothing` if points do not fit to a cone.
"""
function fit(::Type{FittedCone}, p, n, pc, params)
    fcone = fit3pointcone(p, n)
    fcone === nothing && return nothing
    valid_cone = validatecone(fcone, p, n, params)
    return valid_cone
end

## scoring

function compatiblesCone(cone, points, normals, params)
    #@unpack α_cone, ϵ_cone = params
    @extract params : params_cone=cone
    @extract params_cone : α_cone=α ϵ_cone=ϵ
    calcs = (project2cone(cone, points[i]) for i in eachindex(points))

    # eps check
    # c1 = [abs(calcs[i][1]) < ϵ_cone for i in eachindex(calcs)]
    #if cone.outwards
    #    c2=[isparallel(calcs[i][2], normals[i], α_cone) && c1[i] for i in eachindex(calcs)]
    #else
    #    c2=[isparallel(-calcs[i][2], normals[i], α_cone) && c1[i] for i in eachindex(calcs)]
    #end

    zcn = zip(calcs, normals)
    if cone.outwards
        c2=[isparallel(c[2], n, α_cone) && (abs(c[1]) < ϵ_cone) for (c,n) in zcn]
    else
        c2=[isparallel(-c[2], n, α_cone) && (abs(c[1]) < ϵ_cone) for (c,n) in zcn]
    end
    return c2
end

function scorecandidate(pc, candidate::FittedCone, subsetID, params)
    ps = @view pc.vertices[pc.subsets[subsetID]]
    ns = @view pc.normals[pc.subsets[subsetID]]
    ens = @view pc.isenabled[pc.subsets[subsetID]]

    cp = compatiblesCone(candidate, ps, ns, params)
    # verti: összes pont indexe, ami enabled és kompatibilis
    # lenne, ha működne, de inkább a boolean indexelést machináljuk
    inder = cp.&ens
    inpoints = (pc.subsets[subsetID])[inder]
    score = estimatescore(length(pc.subsets[subsetID]), pc.size, length(inpoints))
    return (score, inpoints)
end

## refit

"""
    refit(s::T, pc, params) where {T<:FittedCone}

Refit cone.
"""
function refit(s::T, pc, params) where {T<:FittedCone}
    # TODO: use octree for that
    p = @view pc.vertices[pc.isenabled]
    n = @view pc.normals[pc.isenabled]
    cp = compatiblesCone(s, p, n, params)
    return ExtractedShape(s, ((1:pc.size)[pc.isenabled])[cp])
end
