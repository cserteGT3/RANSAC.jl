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
    @unpack approxnull, α_perpend = params
    is = [[1,2, 3], [2,3, 1], [3,1, 2]]
    ds = [abs(dot(n[i[1]], n[i[2]])) for i in is]
    i = argmin(ds)
    direction = normalize(cross(n[i[1]], n[i[2]]))
    isapprox(norm(direction), 1) || return (false, direction)
    abs(dot(direction, n[i[3]])) < α_perpend && return (true, direction)
    return (false, direction)
end

function project2sketchplane(pcr, indexes, transl_frame, params)
    @unpack α_perpend = params
    # enabled indices
    en_ind = pcr.enabled[indexes]
    z = transl_frame[3]
    # enabled points and normals
    p = @view pcr.vertices[en_ind]
    n = @view pcr.vertices[en_ind]

    projected = Array{SVector{2,Float64},1}(undef, 0)
    inds = Int[]
    for i in eachindex(p)
        # point is not part of the translational surface
        abs(dot(z, n[i])) > α_perpend && continue
        push!(projected, SVector{2}(dot(xv, p[i]), dot(yv, p[i])))
        push!(inds, indexes[i])
    end
    return projected, inds
end

"""
    largestpatch(p, pind, resolution, bb, params)

Indexes of the points in the largest patch.

# Arguments:
- `p`: points in 2D.
- `pind`: index of the points in the pointcloud.
- `resolution`: resolution of the bitmap.
- `bb`: axis aligned bounding box.
- `params::RANSACParameters`: parameters.
"""
function largestpatch(p, pind, resolution, bb, params)
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
        yplace = ceil(Int, abs(p[i][1]-minv[1])/βx)
        xplace = ceil(Int, abs(p[i][2]-minv[2])/βy)
        bitmap[xplace, yplace] = true
        append!(indexmap[xplace, yplace], pind[i])
    end
    largestconncomp(bitmap, indexmap, transl_conn)
end

function fittranslationalsurface(pcr, p, n, params)
    @unpack α_perpend = params
    # Method:
    # 1. van-e közös merőleges? nincs -> break
    ok, dir = transldir(p, n, params)
    ok || return FittedTranslation(false, NaNVec, nothing)
    # 2. összes pont levetítése erre a síkra (egyik pont és a közös normális)
    o = p[1]
    xv = n[1]
    zv = dir
    yv = normalize(cross(zv, xv))
    coordframe = [xv, yv, zv]

    #TODO: okosabb kéne ide
    # most: az első subsetet nézem és abból is csak azt ami enabled
    projected, proj_ind = (pcr, subsets[1], coordframe, params)

    # 3. legkisebb és legnagyobb távolság megnézése
    @unpack mind, maxd = minmaxdistance(projected)
    # 4. AABB
    aabb = findAABB(projected)
    # 5. ha az AABB területe nagyon kicsi -> break
    #TODO: azt kéne inkább nézni, hogy az egyik oldal nagyon kicsi a másikhoz képest=sík
    diagd = norm(aabb[2]-aabb[1])
    diagd < diagthr && return FittedTranslation(false, NaNVec, nothing)
    # 6. legnagyobb összefüggő terület

    # 7. kontúr kiszedése: kell-e, hogy zárt görbe legyen? - szerintem kell -> 2 végpont összekötése
    # 8. kör/egyenes illesztése
    # 9. visszaellenőrzés?
end
