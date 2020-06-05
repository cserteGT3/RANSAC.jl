# fitting

"""
    struct FittedSphere{A<:AbstractArray, R<:Real} <: FittedShape

Sphere primitive, defined by its center and radius.
Also stored, if the normals point outwards of the shape.
"""
struct FittedSphere{A<:AbstractArray, R<:Real} <: FittedShape
    center::A
    radius::R
    outwards::Bool
end

Base.show(io::IO, x::FittedSphere) =
    print(io, """sphere, R: $(x.radius)""")

Base.show(io::IO, ::MIME"text/plain", x::FittedSphere{A, R}) where {A, R} =
    print(io, "FittedSphere{$A, $R}\n", """center: $(x.center), R: $(x.radius), $(x.outwards ? "outwards" : "inwards" )""")

strt(x::FittedSphere) = "sphere"

function defaultshapeparameters(::Type{FittedSphere})
    # `sphere_par`: parameter in sphere fitting
    sp = (ϵ=0.3, α=deg2rad(5), sphere_par=0.02)
    return (sphere=sp,)
end

function fit2pointsphere(v, n, params)
    #@unpack parallelthrdeg, sphere_par = params
    @extract params : params_sphere=sphere
    @extract params_sphere : sphere_par
    @extract params.common : parallelthrdeg
    n1n = normalize(n[1])
    n2n = normalize(n[2])

    b = false
    if abs(dot(n1n, n2n)) > cosd(parallelthrdeg)
        # parallel normals
        # fit it, if bad, will fall out later
        centerp = (v[1]+v[2])/2
        return FittedSphere(centerp, norm(centerp-v[1]), b)
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
            return FittedSphere(centerp, r, b)
        else
            # intersection
            if dot(h, k) > 0
                # point the same direction -> +
                M = v[1] + nh/nk * n1n
                return FittedSphere(M, norm(M-v[1]), b)
            else
                # point in different direction -> -
                M = v[1] - nh/nk * n1n
                return FittedSphere(M, norm(M-v[1]), b)
            end
        end
    end
end

function setsphereOuterity(sp, b)
    FittedSphere(sp.center, sp.radius, b)
end

"""
    fit(::Type{FittedSphere}, p, n, pc, params)

Fit a sphere to 2 points. Additional points and their normals are used to validate the fit.
Return `nothing` if points do not fit to a sphere.
"""
function fit(::Type{FittedSphere}, p, n, pc, params)
    #@unpack ϵ_sphere, α_sphere = params
    @extract params : params_sphere=sphere
    @extract params_sphere : α_sphere=α ϵ_sphere=ϵ
    pl = length(p)
    @assert pl == length(n) "Size must be the same."
    @assert pl > 2 "Size must be at least 3."
    # "forcefit" a sphere
    sp = fit2pointsphere(p, n, params)
    # check if real sphere
    sp === nothing && return nothing
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
    vert_ok == trues(pl) || return nothing
    norm_ok == trues(pl) && return setsphereOuterity(sp, true)
    invnorm_ok == trues(pl) && return setsphereOuterity(sp, false)
    return nothing
end

# bitmapping

function scorecandidate(pc, candidate::FittedSphere, subsetID, params)
    ps = @view pc.vertices[pc.subsets[subsetID]]
    ns = @view pc.normals[pc.subsets[subsetID]]
    ens = @view pc.isenabled[pc.subsets[subsetID]]

    cpl = compatiblesSphere(candidate, ps, ns, params)
    # verti: összes pont indexe, ami enabled és kompatibilis
    # lenne, ha működne, de inkább a boolean indexelést machináljuk
    verti = pc.subsets[subsetID]
    #underEn = uo.under .& cpl
    #overEn = uo.over .& cpl

    #inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
    inpoints = verti[cpl]
    score = estimatescore(length(pc.subsets[subsetID]), pc.size, length(inpoints))
    return (score, inpoints)
end

"""
    compatiblesSphere(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the sphere.

Compatibility is measured with an `eps` distance to the sphere
and an `alpharad` angle to it's normal.
"""
function compatiblesSphere(sphere, points, normals, params)
    #@unpack ϵ_sphere, α_sphere = params
    @extract params : params_sphere=sphere
    @extract params_sphere : α ϵ
    @assert length(points) == length(normals) "Size must be the same."
    # eps check
    o = sphere.center
    R = sphere.radius
    # c1 = [abs(norm(a-o)-R) < ϵ_sphere for a in points]
    # alpha check
    #if sphere.outwards
    #    c2=[isparallel(normalize(points[i]-o), normals[i], α_sphere) && c1[i] for i in eachindex(points)]
    #else
    #    c2=[isparallel(normalize(o-points[i]), normals[i], α_sphere) && c1[i] for i in eachindex(points)]
    #end

    zpn = zip(points, normals)

    if sphere.outwards
        c2=[isparallel(normalize(p-o), n, α) && (abs(norm(p-o)-R) < ϵ) for (p,n) in zpn]
    else
        c2=[isparallel(normalize(o-p), n, α) && (abs(norm(p-o)-R) < ϵ) for (p,n) in zpn]
    end

    # TODO
    # normalize by it's own length or maximum length??? r+3ϵ
    # param_points = [unitdisk2square(normalize(a[1:2]-o[1:2])) for a in points]
    return c2
end

"""
    refit(s::T, pc, params) where {T<:FittedSphere}

Refit sphere.
"""
function refit(s::T, pc, params) where {T<:FittedSphere}
    # TODO: use octree for that
    pcv = @view pc.vertices[pc.isenabled]
    pcn = @view pc.normals[pc.isenabled]
    cpl = compatiblesSphere(s, pcv, pcn, params)
    # verti: összes pont indexe, ami enabled és kompatibilis
    verti = (1:pc.size)[pc.isenabled]
    #underEn = uo.under .& cpl
    #overEn = uo.over .& cpl
    #s.inpoints = count(underEn) >= count(overEn) ? verti[underEn] : verti[overEn]
    return ExtractedShape(s, verti[cpl])
end
