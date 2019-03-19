module samples

include("rodrigues.jl")

using LinearAlgebra: normalize, normalize!, cross, norm
using .RodriguesRotations: rodriguesrad, rodriguesdeg
using AbstractPlotting: Point3f0
using StaticArrays: SVector
using Debugger

export sampleplane, samplecylinder, normalsforplot
export noisifyvertices, noisifynormals
export makemeanexample

"""
    arbitrary_orthogonal(vec)

Create an arbitrary orthogonal vector to `vec`.
"""
function arbitrary_orthogonal(vec)
    @assert size(vec,1) == 3 "Implemented only for 3 dimensional."
    v = normalize(vec)
    b0 = (v[1] <  v[2]) && (v[1] <  v[3]);
    b1 = (v[2] <= v[1]) && (v[2] <  v[3]);
    b2 = (v[3] <= v[1]) && (v[3] <= v[2]);
    return cross(v, normalize!([convert(Int, b0), convert(Int, b1), convert(Int, b2)]))
end

"""
    sampleplane(vp, v1, v2, lengtht, sizet)

Create a tuple of points on a plane and their normals based on a point (`vp`) and
two vectors (`v1`, `v2`).

# Arguments:
- `vp`: a point of the plane.
- `v1`, `v2`: two vectors that define a plane.
- `sizet`: tuple containing the number of samples along each side.
- `lengtht`: tuple containing the length of the plane along each side.
"""
function sampleplane(vp, v1, v2, lengtht, sizet)
    @assert size(v1) == size(v2) == size(vp) "Vector's size must be the same!"
    v1n = normalize(v1)
    v2n = normalize(v2)
    n = cross(v1n, v2n)
    @assert norm(n) != 0 "Vectors can't be collinear!"
    s1, s2 = sizet
    @assert s1 > 0 && s2 > 0 "Should sample more than 0."
    s1l, s2l = lengtht
    pn = [n for i in 1:s1*s2]
    vpn = vp - (s1l*v1n + s2l*v2n)/2
    ps = [vpn + v1n*u/s1l + v2n*v/s2l for u in 1:s1 for v in 1:s2]
    return (ps, pn)
end

"""
    samplecylinder(ax, vp, R, h, sizet)

Create a tuple of points on a cylinder and their normals based on an axis (`ax`),
a nullpoint (`vp`), radius and height.

# Arguments:
- `ax`: axis of the cylinder.
- `vp`: nullpoint of the cylinder.
- `R`: radius of the cylinder.
- `h`: height of the cylinder.
- `sizet`: tuple containing the number of samples along the circumference and the height.
"""
function samplecylinder(ax, vp, R, h, sizet)
    @assert size(ax,1) == size(vp,1) == 3 "Axis and nullpoint must be 3 dimensional."
    @assert R > 0 && h > 0 "Can't construct a line or circle."
    s1, s2 = sizet
    @assert s1 > 1 && s2 > 1 "Should sample more than 1."

    axn = normalize(ax)
    aort = arbitrary_orthogonal(axn)

    ps = Array{eltype(vp)}(undef, 3, s1*s2)
    ns = similar(ps)

    rotMat = rodriguesrad(axn, 2*π/s1)

    ps[:,1] = vp + R.*aort
    ns[:,1] = aort
    for i in 2:s1
        ps[:,i] = rotMat^(i-1) * ps[:,1]
        ns[:,i] = rotMat^(i-1) * ns[:,1]
    end

    linind = LinearIndices((1:s1, 1:s2))
    for i in 2:s2
        for j in 1:s1
            ps[:, linind[j,i]] = ps[:,j] + ((i-1)*h/(s2-1)).*axn
            ns[:, linind[j,i]] = ns[:,j]
        end
    end
    return (ps, ns)
end

"""
    normalsforplot(verts, norms, arrowsize = 0.5)

Create an array of pair of points for plotting the normals with Makie.
Only the direction of the normals is presented, their size not.

# Arguments:
- `arrowsize`: scaling factor for the size of the lines.
"""
function normalsforplot(verts, norms, arrowsize = 0.5)
    @assert size(verts) == size(norms) "They should have the same size."
    as = arrowsize
    return [verts[i] => verts[i] + as .*normalize(norms[i]) for i in 1:length(verts) ]
end

"""
    noisifyvertices(verts, allvs, scalef = 1)

Add gaussian noise to vertices.
Random subset or all vertices can be selected.

# Arguments:
- `allvs::Bool`: adds noise to every vertice if true.
- `scalef::Real`: scale the noise from the [-1,1] interval.
"""
function noisifyvertices(verts, allvs, scalef = 1)
    randis = similar(verts)
    for i in eachindex(randis)
        randis[i] = (2*scalef) .*SVector(rand(eltype(eltype(verts)),3)...) .-1
    end
    allvs && return verts+randis
    bools = rand(Bool, size(verts))
    verts[bools] + randis[bools]
end

"""
    noisifynormals(norms, maxrot)

Add gaussian noise to normals.
Rotates the normals around a random axis with `maxrot` degrees.

# Arguments:
- `maxrot::Real`: maximum rotation in degrees.
"""
function noisifynormals(norms, maxrot)
    retn = similar(norms)
    randis = similar(norms)
    for i in eachindex(randis)
        randis[i] = SVector(rand(eltype(eltype(norms)),3)...)
    end
    randrots = maxrot.*rand(eltype(eltype(norms)), length(norms)).-maxrot/2
    for i in 1:length(retn)
        nr = normalize(randis[i])
        while norm(cross(norms[i], nr)) < 0.1
            nr = normalize(rand(eltype(norms),3))
        end
        crv = cross(nr, norms[i])
        rM = rodriguesdeg(crv, randrots[i])
        retn[i] = rM*norms[i]
    end
    return retn
end

"""
    makemeanexample(nois = false; all = false)

Generate a definitely not random example.
"""
function makemeanexample(nois = false; all = false)
    pp_, pn_ = sampleplane([0,0,0], [0,0,1], [1,3,0], (13.7, 9.58), (96, 112))
    cp_, cn_ = samplecylinder([0,0,1], zeros(3), 15, 15, (100, 150))
    cp2_, cn2_ = samplecylinder([-3,2,1], [0.9, 1.8, -2.3], 7, 9, (150, 100))
    ps_ = hcat(pp_, cp_, cp2_)
    ns_ = hcat(pn_, cn_, cn2_)
    if nois
        noisifyvertices!(ps_, all)
        ns2_ = noisifynormals(ns_, 35)
        nsfp2_ = normalsforplot(ps_, ns2_)
        return ps_, ns_, nsfp2_
    end
    nsfp_ = normalsforplot(ps_, ns_)
    return ps_, ns_, nsfp_
end

end
