module RodriguesRotations

using LinearAlgebra
using ZChop: zchop!

export rodriguesdeg, rodriguesrad

function rodrigues(nv, Θ)
    et = eltype(nv)
    R = zeros(et,3,3)
    R = nv*nv' + cos(Θ).*(Matrix{et}(I, 3,3) - nv*nv') + sin(Θ).*crossprodtensor(nv)
    return zchop!(R)
end

function crossprodtensor(v)
    [0 -v[3] v[2];
    v[3] 0 -v[1];
    -v[2] v[1] 0]
end

function rodriguesrad(nv, ϑ)
    nvn = normalize(nv)
    return rodrigues(nv, ϑ)
end

function rodriguesdeg(nv, ϑ)
    nvn = normalize(nv)
    return rodrigues(nv, deg2rad(ϑ))
end

end
