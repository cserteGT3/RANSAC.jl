# fitting

struct FittedCone{A<:AbstractArray, R<:Real} <: FittedShape
    iscone::Bool
    apex::A
    axis::A
    # opening angle in radians
    opang::R
    outwards::Bool
end

Base.show(io::IO, x::FittedCone) =
    print(io, """$(x.iscylinder ? "o" : "x") cone, ω: $(x.opang)""")

Base.show(io::IO, ::MIME"text/plain", x::FittedCone{A, R}) where {A, R} =
    print(io, """FittedCone{$A, $R}\n$(x.iscone ? "o" : "x") cone, center: $(x.apex), axis: $(x.axis), ω: $(x.opang), $(x.outwards ? "outwards" : "inwards" )""")

strt(x::FittedCone) = "cone"

function isshape(shape::FittedCone)
    return shape.iscone
end

function setconeOuterity(fc, b)
    FittedCone(fc.iscone, fc.apex, fc.axis, fc.opang, b)
end

function fit3pointcone(psok, nsok, params)
        # rank of the coefficient matrix
        p = @view psok[1:3]
        n = @view nsok[1:3]
        r = hcat(n...)
        rank(r) == 3 || return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        ds = [-1*dot(p[i], n[i]) for i in 1:3]
        # rank of the augmented matrix
        rv = hcat(r, ds)
        rank(rv) == 3 || return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        # apex
        ap = r\ds
        # axis
        axis3p = [ap+((v-ap)/norm(v-ap)) for v in p]
        ax = normalize(cross(axis3p[2]-axis3p[1], axis3p[3]-axis3p[1]))
        # opening angle
        angles = [acos(dot(normalize(p[i]-ap), ax)) for i in 1:3]
        opangle = sum(angles)/3
        return FittedCone(true, ap, ax, opangle, true)
end
