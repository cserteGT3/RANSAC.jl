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
    print(io, """$(x.iscone ? "o" : "x") cone, ω: $(x.opang)""")

Base.show(io::IO, ::MIME"text/plain", x::FittedCone{A, R}) where {A, R} =
    print(io, """FittedCone{$A, $R}\n$(x.iscone ? "o" : "x") cone, apex: $(x.apex), axis: $(x.axis), ω: $(x.opang), $(x.outwards ? "outwards" : "inwards" )""")

strt(x::FittedCone) = "cone"

function isshape(shape::FittedCone)
    return shape.iscone
end

function setconeOuterity(fc, b)
    FittedCone(fc.iscone, fc.apex, fc.axis, fc.opang, b)
end

function fit3pointcone(psok, nsok)
        # rank of the coefficient matrix
        p = @view psok[1:3]
        n = @view nsok[1:3]
        r = Array{Float64,2}(undef, (3,3))
        for i in 1:3; for j in 1:3; r[i,j] = n[i][j]; end; end
        rank(r) == 3 || return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        ds = [dot(p[i], n[i]) for i in 1:3]
        # rank of the augmented matrix
        rv = hcat(r, -1 .* ds)
        rank(rv) == 3 || return FittedCone(false, NaNVec, NaNVec, 0.0, false)
        # apex
        ap = SVector{3, Float64}(r\ds)
        # axis
        axis3p = [ap+((v-ap)/norm(v-ap)) for v in p]
        ax = normalize(cross(axis3p[2]-axis3p[1], axis3p[3]-axis3p[1]))
        midp = sum(axis3p)/3
        dirv = normalize(midp-ap)
        ax = dot(ax, dirv) < 0 ? -1*ax : ax
        # opening angle
        angles = [acos(dot(normalize(p[i]-ap), ax)) for i in 1:3]
        opangle = sum(angles)/3
        return FittedCone(true, ap, ax, opangle, true)
end

function fitcone(p, n, params)
        @unpack α_cone, ϵ_cone = params
        fcone = fit3pointcone(p, n)

end
