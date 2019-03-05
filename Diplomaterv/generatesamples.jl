module samples

using LinearAlgebra: normalize, cross, norm, normalize

export sampleplane

"""
    sampleplane(vp, v1, v2, sizet, lengtht)

Create a tuple of points and their normals based on a point (`vp`) and
two vectors (`v1`, `v2`).

# Arguments:
- `vp`: a point of the plane.
- `v1`, `v2`: two vectors that define a plane.
- `sizet`: tuple containing the number of samples along each side.
- `lengtht`: tuple containing the length of the plane along each side.
"""
function sampleplane(vp, v1, v2, sizet, lengtht)
    @assert size(v1) == size(v2) == size(vp) "Vector's size must be the same!"
    v1n = normalize(v1)
    v2n = normalize(v2)
    n = cross(v1n, v2n)
    @assert norm(n) != 0 "Vectors can't be collinear!"
    s1, s2 = sizet
    @assert s1 > 0 && s2 > 0 "Should sample more than 0."
    s1l, s2l = lengtht
    pn = repeat(n, 1, s1*s2)
    ps = similar(pn)
    linind = LinearIndices((1:s1, 1:s2))
    vp = vp - (s1l*v1n + s2l*v2n)/2
    for u in 1:s1
        for v in 1:s2
            ps[:,linind[u,v]] = vp + v1n*u/s1l + v2n*v/s2l
        end
    end
    return (ps, pn)
end

end
