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
    @unpack approxnull, perpendnull = params
    is = [[1,2, 3], [2,3, 1], [3,1, 2]]
    ds = [abs(dot(n[i[1]], n[i[2]])) for i in is]
    i = argmin(ds)
    direction = normalize(cross(n[i[1]], n[i[2]]))
    isapprox(norm(direction), 1) || return (false, direction)
    dot(direction, n[i[3]]) < perpendnull && return (true, direction)
    return (false, direction)
end

function fittranslationalsurface(pcr, p, n, params)
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

    # okosabb kéne ide
    # most: az első subsetet nézem
    points = @view pcr.vertices[subsets[1]]
    normals = @view pcr.vertices[subsets[1]]

    projected = Array{SVector{2,Float64},1}(undef, size(points))
    for i in eachindex(points)
        v = points[i]-o
        projected[i] = SVector{2}(dot(xv, v), dot(yv, v))
    end
    # 3. legkisebb és legnagyobb távolság megnézése
    @unpack mind, maxd = minmaxdistance(projected)
    # 4. AABB
    aabb = findAABB(projected)
    # 5. ha az AABB területe nagyon kicsi -> break
    diagd = norm(aabb[1]-aabb[2])
    diagd < diagthr && return FittedTranslation(false, NaNVec, nothing)
    # 6. legnagyobb összefüggő terület
    # 7. kontúr kiszedése: kell-e, hogy zárt görbe legyen? - szerintem kell -> 2 végpont összekötése
    # 8. kör/egyenes illesztése
    # 9. visszaellenőrzés?
end
