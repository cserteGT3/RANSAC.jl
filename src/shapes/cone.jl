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
    print(io, """FittedCylinder{$A, $R}\n$(x.iscylinder ? "o" : "x") cone, center: $(x.center), axis: $(x.axis), ω: $(x.opang), $(x.outwards ? "outwards" : "inwards" )""")

strt(x::FittedCone) = "cone"

function isshape(shape::FittedCone)
    return shape.iscone
end

function setconeOuterity(fc, b)
    FittedCone(fc.iscone, fc.apex, fc.axis, fc.opang, b)
end

function fit3pointcone(p, n, params)
        # rank of the coefficient matrix
        r = hcat(n...)
        rank(r) == 3 || return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        ds = [-1*dot(p[i], n[i]) for i in 1:3]
        # rank of the augmented matrix
        rv = hcat(r, ds)
        rank(rv) == 3 || return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        # apex
        ap = r\ds
        # axis
        threeplane = fitplane(p, n, params)
        if ! isshape(threeplane)
            @debug "Fitted plane is not a plane!"
            return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        end
        ax = threeplane.normal
        # opening angle
        angles = [acos(dot(normalize(p[i]-ap), ax)) for i in 1:3]
        opangle = sum(angles)/3
        return FittedCone(true, ap, ax, opangle, true)
end
