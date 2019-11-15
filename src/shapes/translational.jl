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

struct FittedTranslational{A<:AbstractArray} <: FittedShape
    istranslational::Bool
    coordframe::A
    contour
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

function fittranslationalsurface(pcr, p, n, params)
    @unpack α_perpend, diagthr, max_group_num = params
    @unpack max_contour_it, thinning_par = params
    # Method:
    # 1. van-e közös merőleges? nincs -> break
    ok, dir = transldir(p, n, params)
    ok || return FittedTranslational(false, NaNVec, nothing)
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
    diagd = norm(aabb[2]-aabb[1])
    diagd < diagthr && return FittedTranslational(false, NaNVec, nothing)
    # 6. összefüggő kontúrok
    thr = 1.5*avgd
    maxit = max_contour_it
    spatchs = CompressedFinder([], 0)
    while true
        spatchs = segmentpatches(projected, thr)
        spatchs.groups <= max_group_num && break
        maxit < 1 && return FittedTranslational(false, NaNVec, nothing)
        maxi -= 1
    end
    # hereby spatchs should contain maximum max_group_num of patches
    # 7. kontúr kiszedése: kell-e, hogy zárt görbe legyen? - szerintem kell -> 2 végpont összekötése
    results = Array{FittedTranslational,1}(undef, spatchs.groups)
    for i in 1:spatchs.groups
        patch_p = projected[spatchs.ids[i]]
        tris = delaunay(patch_p)
        all_edges = to_edges!(tris)
        weights = [norm(patch_p[e[1]] - patch_p[e[2]]) for e in all_edges]
        tree = spanning_tree(all_edges, weights)
        thinned, chunks = thinning(p, tree, 2.0)
        # put them into a FittedTranslational
    end

    # 8. kör/egyenes illesztése
    # - skipp
    # 9. visszaellenőrzés?
    # -skipp? :D 
end
