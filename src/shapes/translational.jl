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

struct FittedTranslational <: FittedShape
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

Base.show(io::IO, x::FittedTranslational) =
    print(io, """$(x.istranslational ? "o" : "x") extruded""")

Base.show(io::IO, ::MIME"text/plain", x::FittedTranslational) =
    print(io, """FittedTranslational\n$(x.istranslational ? "o" : "x") extruded""")

strt(x::FittedTranslational) = "extruded"

isshape(shape::FittedTranslational) = shape.istranslational

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

function project2sketchplane(pcr, indexes, transl_frame, params)
    @unpack α_perpend = params
    # enabled indices
    en_ind = pcr.isenabled[indexes]
    xv = transl_frame[1]
    yv = transl_frame[2]
    z = transl_frame[3]
    # enabled points and normals
    p = @view pcr.vertices[en_ind]
    n = @view pcr.normals[en_ind]

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
    largestpatch_depr(p, pind, resolution, bb, params)

Indexes of the points in the largest patch.
This should not be used.

# Arguments:
- `p`: points in 2D.
- `pind`: index of the points in the pointcloud.
- `resolution`: resolution of the bitmap.
- `bb`: axis aligned bounding box.
- `params::RANSACParameters`: parameters.
"""
function largestpatch_depr(p, pind, resolution, bb, params)
    @unpack transl_conn = params
    minv = bb[1]
    maxv = bb[2]

    # size of the bitmap
    xs = trunc(Int, (maxv[1]-minv[1])/resolution)
    ys = trunc(Int, (maxv[2]-minv[2])/resolution)
    @assert xs > 0 && ys > 0 "max-min must be > 0. Check the code! xs: $xs, ys:$ys"
    # size of pixels
    βx = (maxv[1]-minv[1])/xs
    βy = (maxv[2]-minv[2])/ys

    bitmap = falses(xs, ys)
    indexmap = [Int[] for i in 1:xs, j in 1:ys]
    for i in eachindex(p)
        # coordinate frame is transposed
        xplace = ceil(Int, abs(p[i][1]-minv[1])/βx)
        yplace = ceil(Int, abs(p[i][2]-minv[2])/βy)
        xplace = xplace == 0 ? 1 : xplace
        yplace = yplace == 0 ? 1 : yplace
        bitmap[xplace, yplace] = true
        append!(indexmap[xplace, yplace], pind[i])
    end
    return p, bitmap
    # return largestconncomp(bitmap, indexmap, transl_conn)
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
"""
function midpoint(A, i)
    i == lastindex(A) && return (A[1]+A[end])/2
    return (A[i]+A[i+1])/2
end

"""
    twopointnormal(a)

Compute the normal of a segment. Inwards/outwards is not considered here.
The normal points towards "left".
"""
function twopointnormal(a)
    dirv = normalize(a[2]-a[1])
    return normalize(convert(eltype(a), [-dirv[2], dirv[1]]))
end

"""
    segmentnormal(A, i)

Compute the normal of the i-th segment of a linesegment-list.
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
    contournormal(shape::FittedTranslational, i)

Compute the normal of the i-th segment of a FittedTranslational shape.
"""
function contournormal(shape::FittedTranslational, i)
    return shape.flipnormal .* segmentnormal(shape.contour, i)
end

"""
    outwardsnormal(shape::FittedTranslational, i)

Compute the normal of the i-th segment of a FittedTranslational shape.
This normal points always outwards.
"""
function outwardsnormal(shape::FittedTranslational, i)
    return shape.outwards .* segmentnormal(shape.contour, i)
end

"""
    distance2onesegment(point, A, i)

Compute distance of `point` to the i-th segment in A.
"""
function distance2onesegment(point, A, i)
    nv = segmentnormal(A, i)
    topoint = point-A[i]
    return dot(topoint, nv)
end

"""
    dist2segment(point, A)

Compute the shortest distance from `point` to the linesegments `A`.
Also return the index of that segment.
"""
function dist2segment(point, A)
    leastd = distance2onesegment(point, A, 1)
    size(A,1) == 1 && return leastd
    best = 1
    for i in 2:size(A,1)
        d = distance2onesegment(point, A, i)
        if d < leastd
            leastd = d
            best = i
        end
    end
    return (leastd, best)
end

function validatetrans(candidate, ps, ns, params)
    @unpack α_transl , ϵ_transl = params
    calcs = [distandnormal2segment(p, candidate.contour) for p in ps]
    #TODO:
    return nothing
end

"""
    normaldirs(segments, points, normals, center, params)

Compute the normal directions.
This is based on the points, which are closer than `ϵ_transl`.
Return a boolean first that indicates that the normals direct to the "same direction".
"""
function normaldirs(segments, points, normals, center, params)
    fff() = (false, false, false)
    @unpack ϵ_transl, min_normal_num = params
    calcs = [dist2segment(p, segments) for p in points]
    compats = [abs(calcs[i][1]) < ϵ_transl for i in eachindex(calcs)]
    compatsize = size(compats, 1)
    compatsize == 0 && return fff()
    # later working with points[compats]
    psize = size(points,1)
    flipnormal = Vector{Bool}(undef, psize)
    outwards = Vector{Bool}(undef, psize)
    for i in 1:psize
        # continue if not compatible
        compats[i] || continue
        # this is the fitted normal
        contour_n = segmentnormal(segments, calcs[i][2])
        tocenter = normalize(center-midpoint(segments, calcs[i][2]))

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
    thisoutwn = count(thisoutw)
    # can't agree on outwards
    (outwr > min_normal_num) || (outwr <= 1-min_normal_num) || return fff()
    @show thisoutwn
    @show thisoutwn/compatsize

    # can't agree on flipnormals
    thisflip = @view flipnormal[compats]
    flipr = count(thisflip)/compatsize
    thisflipn = count(thisflip)
    (flipr > min_normal_num) || (flipr <= 1-min_normal_num) || return fff()

    # this means, that the computed normals must be turned to direct outside
    outwb = outwr > min_normal_num ? 1 : -1

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
    # most: az első subsetet nézem és abból is csak azt ami enabled
    #projected, proj_ind = project2sketchplane(pcr, pcr.subsets[1], coordframe, params)
    projected, proj_ind = project2sketchplane(pcr, 1:size(pcr.vertices,1), coordframe, params)

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
        patch_p = projected[findall(x->x==i, spatchs.ids)]
        #return patch_p, nothing
        #tris = delaunay(patch_p)
        #all_edges = to_edges!(tris)
        #weights = [norm(patch_p[e[1]] - patch_p[e[2]]) for e in all_edges]
        #tree = spanning_tree(all_edges, weights)
        #thinned, chunks = thinning(p, tree, 2.0)
        thinned, chunks = thinning_slow(patch_p, ϵ_transl/2)
        # put them into a FittedTranslational
        closed = [SVector{2,Float64}(th) for th in thinned]
        c = centroid(closed)

        used_n = @view pcr.normals[proj_ind]
        proj_n = project2sketchplane(used_n, coordframe)

        isok, outw, flips = normaldirs(closed, projected, proj_n, c, params)
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

## scoring

function compatiblesTranslational(shape, points, normals, params)
    @unpack α_transl, ϵ_transl = params

    # project to plane
    ps = project2sketchplane(points, shape.coordframe)
    ns = project2sketchplane(normals, shape.coordframe)

    calcs = [dist2segment(p, shape.contour) for p in ps]

    #eps check
    c1 = [abs(calcs[i][1]) < ϵ_transl for i in eachindex(points)]

    #alpha check
    c2 = Vector{Bool}(undef, size(points))

    flips = shape.flipnormal ? -1 : 1
    for i in eachindex(c2)
        comp_n = -1*segmentnormal(shape.contour, calcs[i][2])
        c2[i] = isparallel(comp_n, ns[i], α_transl)
    end
    return c1.&c2
end


function scorecandidate(pc, candidate::ShapeCandidate{T}, subsetID, params) where {T<:FittedTranslational}
    ps = @view pc.vertices[pc.subsets[subsetID]]
    ns = @view pc.normals[pc.subsets[subsetID]]
    ens = @view pc.isenabled[pc.subsets[subsetID]]
    # verti: összes pont indexe, ami enabled és kompatibilis
    # lenne, ha működne, de inkább a boolean indexelést machináljuk
    cp = compatiblesCone(candidate.shape, ps, ns, params)
    inder = cp.&ens
    inpoints = (pc.subsets[subsetID])[inder]
    score = estimatescore(length(pc.subsets[subsetID]), pc.size, length(inpoints))
    pc.levelscore[candidate.octree_lev] += E(score)
    return ScoredShape(candidate, score, inpoints)
end
