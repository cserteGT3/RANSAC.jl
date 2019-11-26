# Method:
# 1. van-e közös merőleges? nincs -> break
# 2. összes pont levetítése erre a síkra (egyik pont és a közös normális)
# 3. legkisebb és legnagyobb távolság megnézése
# 4. AABB
# 5. ha az AABB területe nagyon kicsi -> break
# 6. legnagyobb összefüggő terület
# 7. kontúr kiszedése: kell-e, hogy zárt görbe legyen? - szerintem kell -> 2 végpont összekötése
# 8. kör/egyenes illesztése
# 9. visszaellenőrzés?

abstract type AbstractTranslationalSurface <: FittedShape end

struct FittedTranslational <: AbstractTranslationalSurface
    istranslational::Bool
    coordframe
    contourindexes
    subsetnum::Int
end

struct ExtractedTranslational <: AbstractTranslationalSurface
    istranslational::Bool
    coordframe
    contour
    # center of gravity
    center
    # normal of the contour is parallel to the direction
    # towards the center of the contour?
    # == should flip the computed normals to direct outwards?
    # this is used in e.g. CSGBuilding
    # true means, that the computed normals must be turned to direct outside
    outwards::Int
    # should the computed normal be flipped to match the measured points
    # this is used in this package to ensure that in/outwards is correct
    # true means that computed normals must be turned to match the measured points
    flipnormal::Int
end

Base.show(io::IO, x::AbstractTranslationalSurface) =
    print(io, """$(x.istranslational ? "o" : "x") translational""")

Base.show(io::IO, ::MIME"text/plain", x::FittedTranslational) =
    print(io, """FittedTranslational\n$(x.istranslational ? "o" : "x") translational""")

Base.show(io::IO, ::MIME"text/plain", x::ExtractedTranslational) =
    print(io, """ExtractedTranslational\n$(x.istranslational ? "o" : "x") translational""")

strt(x::AbstractTranslationalSurface) = "translational"

isshape(shape::AbstractTranslationalSurface) = shape.istranslational

function transldir(p, n, params)
    # 1. van-e közös merőleges? nincs -> break
    @unpack α_perpend = params
    is = [[1,2, 3], [2,3, 1], [3,1, 2]]
    ds = [abs(dot(n[i[1]], n[i[2]])) for i in is]
    sel_i = is[argmin(ds)]
    direction = normalize(cross(n[sel_i[1]], n[sel_i[2]]))
    isapprox(norm(direction), 1) || return (false, direction)
    abs(dot(direction, n[sel_i[3]])) < α_perpend && return (true, direction)
    return (false, direction)
end

"""
    project2sketchplane(pcr, indexes, transl_frame, params)

Project points defined by `indexes` (`pcr.vertices[indexes]`) to the plane defined by `transl_frame`.
All the points should be enabled.
Only those points are considered, whose normal is perpendicular to the plane.
"""
function project2sketchplane(pcr, indexes, transl_frame, params)
    @unpack α_perpend = params
    xv = transl_frame[1]
    yv = transl_frame[2]
    z = transl_frame[3]
    # enabled points and normals
    p = @view pcr.vertices[indexes]
    n = @view pcr.normals[indexes]

    projected = Array{SVector{2,Float64},1}(undef, 0)
    inds = Int[]
    for i in eachindex(p)
        # point is not part of the translational surface
        abs(dot(z, n[i])) > α_perpend && continue
        nv1 = dot(xv, p[i])
        nv2 = dot(yv, p[i])
        nvv = SVector{2, Float64}(nv1, nv2)
        push!(projected, nvv)
        push!(inds, indexes[i])
    end
    return projected, inds
end

"""
    project2sketchplane(points, transl_frame)

Just project the points to the plane.
"""
function project2sketchplane(points, transl_frame)
    xv = transl_frame[1]
    yv = transl_frame[2]
    projd = [SVector{2,Float64}(dot(xv, p), dot(yv, p)) for p in points]
    return projd
end

"""
    segmentpatches(points, ϵ_inrange)

Return connected patches of a pointcloud
"""
function segmentpatches(points, ϵ_inrange)
    btree = BallTree(points)
    uf = UnionFinder(size(points,1))

    for i in eachindex(points)
        inr = inrange(btree, points[i], ϵ_inrange, true)
        for j in eachindex(inr)
            # the point itself is always in range
            i == inr[j] && continue
            union!(uf, i, inr[j])
        end
    end
    return CompressedFinder(uf)
end

## distance and normal computation of linesegments

"""
    midpoint(A, i)

Compute the midpoint of the i-th segment of a linesegment-list.
`A`: list of points.
"""
function midpoint(A, i)
    i == lastindex(A) && return (A[1]+A[end])/2
    return (A[i]+A[i+1])/2
end

"""
    nearestpoint(point, A)

Return the index of the nearest point from A to `point`.
`A`: list of points.
"""
function nearestpoint(point, A)
    di = norm(point-A[1])
    i = 1
    for j in eachindex(A)
        dj = norm(point-A[j])
        if dj < di
            di = dj
            i = j
        end
    end
    return di, i
end

"""
    twopointnormal(a)

Compute the normal of a segment. Inwards/outwards is not considered here.
The normal points towards "left".
`A`: list of points.
"""
function twopointnormal(a)
    dirv = normalize(a[2]-a[1])
    return normalize(convert(eltype(a), [-dirv[2], dirv[1]]))
end

"""
    segmentnormal(A, i)

Compute the normal of the i-th segment of a linesegment-list.
`A`: list of points.
"""
function segmentnormal(A, i)
    if i == lastindex(A)
        a = [A[end], A[1]]
        return twopointnormal(a)
    end
    b = @view A[i:i+1]
    return twopointnormal(b)
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

# THIS SHOULD BE USED!!!!!!!!!!!!!!!!
"""
    dn2contour(point, contour)

Distance: distance from the nearest point of the contour point.
"""
function dn2contour(point, contour)
    d, i = nearestpoint(point, contour)
    i_1 = i==1 ? lastindex(contour) : i-1
    # original:
    #pn_ = (segmentnormal(contour, i_1)+segmentnormal(shape.contour, i))/2
    dss = [contourdistance(point, contour, is) for is in [i_1, i]]
    # 1. this leaves unwanted things
    #pn_ = segmentnormal(contour, [i_1, i][argmin(dss)])
    # 2.
    sdss = sum(dss)
    if isapprox(sdss, 0)
        # can't divide by 0
        pn_ = segmentnormal(contour, i)
    else
        pn_ = (dss[2]*segmentnormal(contour, i_1)+dss[1]*segmentnormal(contour, i))/sdss
    end
    return (d, pn_, i)
end

# THIS SHOULD BE USED!!!!!!!!!!!!!!!!
"""
    dn2shape_outw(point, shape)

Distance: distance from the nearest point of the shape.contour point.
Normal points always outwards. Use in CSG.
"""
function dn2shape_outw(point, shape)
    d, pn_, i = dn2contour(point, shape.contour)
    pn = shape.outwards * pn_
    dotp = dot(pn, point-shape.contour[i])
    signi = dotp < 0 ? -1 : 1
    return (signi*d, pn, i)
end

# THIS SHOULD BE USED!!!!!!!!!!!!!!!!
"""
    dn2shape_contour(point, shape)

Distance: distance from the nearest point of the shape.contour point.
Use in RANSAC.
"""
function dn2shape_contour(point, shape)
    d, pn_, i = dn2contour(point, shape.contour)
    pn = shape.flipnormal * pn_
    dotp = dot(pn, point-shape.contour[i])
    signi = dotp < 0 ? -1 : 1
    return (signi*d, pn, i)
end

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

## not line but segment distances

"""
    segmentdistance(q, ab)

Distance from `q` to segment, which is just two points.
This distance is an absolute value.
(Distance to segment, not line!)
"""
function segmentdistance(q, ab)
    a = ab[1]
    b = ab[2]
    if isapprox(norm(b-a), 0)
        @warn "$ab is just a point, not a linesegment."
        return norm(q-a)
    end
    v = normalize(b-a)
    a2qv = dot(q-a,v)*v
    qv = a+a2qv
    dv = abs.(b-a)
    # indc = abs(bx-ax) > abs(by-ay) ? 1 : 2
    i = dv[1] > dv[2] ? 1 : 2
    t = (q[i]-a[i])/(b[i]-a[i])
    if t < 0
        return norm(q-a)
    elseif t > 1
        return norm(q-b)
    else
        return norm(q-qv)
    end
end

"""
    contourdistance(p, contour, i)

Distance of `p` from the `i`-th segment of `contour`.
Uses `segmentdistance`.
"""
function contourdistance(p, contour, i)
    if i == lastindex(contour)
        a = [contour[end], contour[1]]
        return segmentdistance(p, a)
    end
    b = @view contour[i:i+1]
    return segmentdistance(p, b)
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

"""
    normaldirs(segments, points, normals, center, params)

Compute the normal directions.
This is based on the points, which are closer than `ϵ_transl`.
Return a boolean first that indicates that the normals direct to the "same direction".
"""
function normaldirs(segments, points, normals, center, params)
    @assert size(points) == size(normals)
    fff() = (false, false, false)
    @unpack ϵ_transl, min_normal_num = params
    #calcs = [dist2segment(p, segments) for p in points]
    calcs = [dn2contour(p, segments) for p in points]
    compats = [abs(calcs[i][1]) < ϵ_transl for i in eachindex(calcs)]
    compatsize = count(compats)
    compatsize == 0 && return fff()
    # later working with points[compats]
    psize = size(points,1)
    flipnormal = Vector{Bool}(undef, psize)
    outwards = Vector{Bool}(undef, psize)
    for i in 1:psize
        # continue if not compatible
        compats[i] || continue
        # this is the fitted normal
        #contour_n = segmentnormal(segments, calcs[i][2])
        contour_n = calcs[i][2]
        tocenter = normalize(center-midpoint(segments, calcs[i][3]))

        # normal of the contour is parallel to the direction
        # towards the center of the contour?
        # == should flip the computed normals to direct outwards?
        # this is used in e.g. CSGBuilding
        outwards[i] = dot(contour_n, tocenter) < 0.0 ? false : true

        # should the computed normal be flipped to match the measured points
        # this is used in this package to ensure that in/outwards is correct
        flipnormal[i] = dot(contour_n, normals[i]) < 0.0 ? true : false
    end
    thisoutw = @view outwards[compats]
    outwr = count(thisoutw)/compatsize
    # can't agree on outwards
    (outwr > min_normal_num) || (outwr <= 1-min_normal_num) || return fff()

    # can't agree on flipnormals
    thisflip = @view flipnormal[compats]
    flipr = count(thisflip)/compatsize
    (flipr > min_normal_num) || (flipr <= 1-min_normal_num) || return fff()
    # this means, that the computed normals must be turned to direct outside
    outwb = outwr > min_normal_num ? -1 : 1

    # this means that computed normals must be turned to match the measured points
    flipn = flipr > min_normal_num ? -1 : 1

    return (true, outwb, flipn)
end

function fittranslationalsurface(pcr, p, n, params)
    @unpack α_perpend, diagthr, max_group_num = params
    @unpack max_contour_it, thinning_par, ϵ_transl = params
    # Method:
    # 1. van-e közös merőleges? nincs -> break
    ok, dir = transldir(p, n, params)
    ok || return nothing
    # 2. összes pont levetítése erre a síkra (egyik pont és a közös normális)
    o = p[1]
    xv = n[1]
    zv = dir
    yv = normalize(cross(zv, xv))
    coordframe = [xv, yv, zv]

    #TODO: okosabb kéne ide
    # these are the enabled points from the first subset
    ien = pcr.isenabled
    subsnum = 1
    sbs = pcr.subsets[subsnum]
    # index in isenabled with the subset, to get those who are enabled in the subset
    # then use it to index into the subset
    used_i = sbs[ien[sbs]]
    projected, proj_ind = project2sketchplane(pcr, used_i, coordframe, params)
    size(projected, 1) < 2 && return nothing
    # 4. AABB
    aabb = findAABB(projected)
    sidelength = aabb[2]-aabb[1]
    #TODO: azt kéne inkább nézni, hogy az egyik oldal nagyon kicsi a másikhoz képest=sík
    # if one of the side's length is <<< then the other -> nothing
    sidelength[1] < 0.01*sidelength[2] && return nothing
    sidelength[2] < 0.01*sidelength[1] && return nothing
    # 6. összefüggő kontúrok
    thr = ϵ_transl
    maxit = max_contour_it
    spatchs = segmentpatches(projected, thr)
    while spatchs.groups > max_group_num
        maxit < 1 && return nothing
        maxit -= 1
        thr = 0.9*thr
        spatchs = segmentpatches(projected, thr)
        #spatchs.groups <= max_group_num && break
    end
    # hereby spatchs should contain maximum max_group_num of patches
    fitresults = Array{FittedTranslational,1}(undef, 0)
    for i in 1:spatchs.groups
        # i if part of the i-th group
        cur_group = findall(x->x==i, spatchs.ids)
        patch_indexes = proj_ind[cur_group]
        ft = FittedTranslational(true, coordframe, patch_indexes, subsnum)
        push!(fitresults, ft)
    end
    return fitresults
end

#=
function fittranslationalsurface_orig(pcr, p, n, params)
    @unpack α_perpend, diagthr, max_group_num = params
    @unpack max_contour_it, thinning_par, ϵ_transl = params
    # Method:
    # 1. van-e közös merőleges? nincs -> break
    ok, dir = transldir(p, n, params)
    ok || return nothing
    # 2. összes pont levetítése erre a síkra (egyik pont és a közös normális)
    o = p[1]
    xv = n[1]
    zv = dir
    yv = normalize(cross(zv, xv))
    coordframe = [xv, yv, zv]

    #TODO: okosabb kéne ide
    # these are the enabled points from the first subset
    ien = pcr.isenabled
    sbs = pcr.subsets[1]
    # index in isenabled with the subset, to get those who are enabled in the subset
    # then use it to index into the subset
    used_i = sbs[ien[sbs]]
    projected, proj_ind = project2sketchplane(pcr, used_i, coordframe, params)

    # 3. legkisebb és legnagyobb távolság megnézése
    @unpack mind, maxd, avgd = minmaxdistance(projected)
    # 4. AABB
    aabb = findAABB(projected)
    # 5. ha az AABB területe nagyon kicsi -> break
    #TODO: azt kéne inkább nézni, hogy az egyik oldal nagyon kicsi a másikhoz képest=sík
    #diagd = norm(aabb[2]-aabb[1])
    #diagd < diagthr && return nothing
    # 6. összefüggő kontúrok
    thr = ϵ_transl
    maxit = max_contour_it
    spatchs = segmentpatches(projected, thr)
    while spatchs.groups > max_group_num
        maxit < 1 && return nothing
        maxit -= 1
        thr = 0.9*thr
        spatchs = segmentpatches(projected, thr)
        #spatchs.groups <= max_group_num && break
    end
    # hereby spatchs should contain maximum max_group_num of patches
    # 7. kontúr kiszedése: kell-e, hogy zárt görbe legyen? - szerintem kell -> 2 végpont összekötése
    fitresults = Array{FittedTranslational,1}(undef, 0)
    for i in 1:spatchs.groups
        # if spatchs.groups==1 -> .ids is not Array of Array, only Array

        # 1 if part of the current group
        cur_group = findall(x->x==i, spatchs.ids)
        patch_p = projected[cur_group]
        #tris = delaunay(patch_p)
        #all_edges = to_edges!(tris)
        #weights = [norm(patch_p[e[1]] - patch_p[e[2]]) for e in all_edges]
        #tree = spanning_tree(all_edges, weights)
        #thinned, chunks = thinning(p, tree, 2.0)
        #println("pathc")
        #@show patch_p
        println("time slow:")
        thinned, chunks = @time thinning_slow(patch_p, ϵ_transl/2)
        # put them into a FittedTranslational
        closed = [SVector{2,Float64}(th) for th in thinned]
        c = centroid(closed)

        used_n = @view pcr.normals[proj_ind[cur_group]]
        proj_n = project2sketchplane(used_n, coordframe)

        isok, outw, flips = normaldirs(closed, patch_p, proj_n, c, params)
        if !isok
            @debug "normals not ok"
            continue
        end
        ft = FittedTranslational(true, coordframe, closed, c, outw, flips)
        push!(fitresults, ft)
    end

    # 8. kör/egyenes illesztése
    # - skipp
    # 9. visszaellenőrzés?
    # validation is currently skipped
    isempty(fitresults) && return nothing
    return fitresults
end
=#

## scoring

function scorecandidate(pc, candidate::ShapeCandidate{T}, subsetID, params) where {T<:FittedTranslational}
    inpoints = candidate.shape.contourindexes
    subsID = candidate.shape.subsetnum
    score = estimatescore(length(pc.subsets[subsID]), pc.size, length(inpoints))
    pc.levelscore[candidate.octree_lev] += E(score)
    return ScoredShape(candidate, score, inpoints)
end

## refit

function compatiblesTranslational(shape, points, normals, params)
    @unpack α_transl, ϵ_transl = params

    # project to plane
    ps = project2sketchplane(points, shape.coordframe)
    ns = project2sketchplane(normals, shape.coordframe)

    calcs = [dn2shape_contour(p, shape) for p in ps]

    #eps check
    c1 = [abs(calcs[i][1]) < ϵ_transl for i in eachindex(points)]

    #alpha check
    c2 = Vector{Bool}(undef, size(points))

    for i in eachindex(c2)
        #comp_n = contournormal(shape, calcs[i][2])
        comp_n = calcs[i][2]
        c2[i] = isparallel(comp_n, ns[i], α_transl)
    end
    return c1.&c2
end

"""
    refittransl(s, pc, params)

Refit translational. Only s.inpoints is updated.
"""
function refittransl(s, pc, params)
    @unpack ϵ_transl, force_transl = params
    transl = s.candidate.shape
    cf = transl.coordframe
    cidxs = transl.contourindexes

    # 1. project points & normals
    old_p = @view pc.vertices[cidxs]
    old_n = @view pc.normals[cidxs]
    o_pp = project2sketchplane(old_p, cf)
    o_np = project2sketchplane(old_n, cf)

    # 2. thinning
    thinned, _ = thinning_slow(o_pp, ϵ_transl/2)
    closed = [SVector{2,Float64}(th) for th in thinned]
    c = centroid(closed)

    # 3. normaldirs()
    isok, outw, flips = normaldirs(closed, o_pp, o_np, c, params)
    (isok || force_transl) || return nothing
    et = ExtractedTranslational(true, cf, closed, c, outw, flips)

    # 4. search for all enabled and compatibel points
    # TODO: use octree for that
    p = @view pc.vertices[pc.isenabled]
    n = @view pc.normals[pc.isenabled]
    cp = compatiblesTranslational(et, p, n, params)
    ip = ((1:pc.size)[pc.isenabled])[cp]
    cand = ShapeCandidate(et, s.candidate.octree_lev)
    return ScoredShape(cand, s.score, ip)
end
