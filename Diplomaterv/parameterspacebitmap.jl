module ParameterspaceBitmap

include("utilities.jl")
include("fitting.jl")

using StaticArrays: SVector, MVector
using LinearAlgebra: cross, dot, normalize, normalize!, norm
using Logging

using .Fitting: FittedPlane, FittedSphere
using .Utilities

export project2plane, compatiblesPlane
export bitmapparameters
export compatiblesSphere

"""
    project2plane(plane, points)

Project `points` on to the `plane`.

# Arguments:
- 'plane::FittedPlane': a plane.
- 'points::AbstractArray': an array of points (array-of-arrays).
"""
function project2plane(plane, points)
    # create a coordinate frame
    # z is the plane's normal
    o_z = normalize(plane.normal)
    # x is a random orthogonal vector (in the plane)
    o_x = normalize(arbitrary_orthogonal(o_z))
    # y is created so it's a right hand coord. frame
    o_y = normalize(cross(o_z, o_x))
    answer = similar(points)
    # get the coordinates of the points in the prev. created coord. frame
    for i in eachindex(points)
        v = points[i]-plane.point
        answer[i] = eltype(answer)(dot(o_x,v), dot(o_y,v), dot(o_z,v))
    end
    answer
end

"""
    compatiblesPlane(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the plane.
Give back the projected points too for parameter space magic.

Compatibility is measured with an `eps` distance to the plane and an `alpharad` angle to it's normal.
"""
function compatiblesPlane(plane, points, normals, eps, alpharad)
    @assert length(points) == length(normals) "Size must be the same."
    projecteds = project2plane(plane, points)
    # eps check
    c1 = [abs(a[3]) < eps for a in projecteds]
    # alpha check
    c2 = [isparallel(plane.normal, normals[i], alpharad) && c1[i] for i in eachindex(normals)]
    # projecteds[c2] are the compatible points
    return c2, projecteds
end

"""
    bitmapparameters(parameters, compatibility, beta, idsource)

Bitmap the compatible parameters. An `idsource` is used to create the indexer map.

# Arguments:
- `parameters::AbstractArray`: contains the parameter pairs.
- `compatibility::AbstractArray`: indicates if the i-th parameter-pair is compatible.
- `beta<:Real`: resolution of the bitmap.
- `idsource::AbstractArray`: source of the indexer map.
"""
function bitmapparameters(parameters, compatibility, beta, idsource)
    @assert length(parameters) == length(compatibility) == length(idsource) "Everything must have the same length."
    minv, maxv = findAABB(parameters)
    xs = round(Int, (maxv[1]-minv[1])/beta)
    ys = round(Int, (maxv[2]-minv[2])/beta)
    @assert xs > 0 && ys > 0 "max-min should be positive. Check the code!"
    βx = (maxv[1]-minv[1])/xs
    βy = (maxv[2]-minv[2])/ys

    bitmap = falses(xs,ys)
    idxmap = zeros(Int,xs,ys)
    for i in eachindex(parameters)
        if compatibility[i]
            xplace = ceil(Int,(parameters[i][1]-minv[1])/βx)
            yplace = ceil(Int,(parameters[i][2]-minv[2])/βy)
            if xplace!=0 && yplace!=0
                if !bitmap[xplace, yplace]
                    bitmap[xplace, yplace] = true
                    idxmap[xplace, yplace] = idsource[i]
                else
                    @warn "Should not project 2 or more points into one pixel! $i-th iteration."
                end
            else
                @warn "that would be boundserror at $i !"
            end
        end
    end
    return bitmap, idxmap, (βx, βy)
end

"""
    bitmapparameters(parameters, compatibility, beta)

Bitmap the compatible parameters. The `idsource` is: `1:length(parameters)`.

# Arguments:
- `parameters::AbstractArray`: contains the parameter pairs.
- `compatibility::AbstractArray`: indicates if the i-th parameter-pair is compatible.
- `beta<:Real`: resolution of the bitmap.
"""
function bitmapparameters(parameters, compatibility, beta)
    return bitmapparameters(parameters, compatibility, beta, 1:length(parameters))
end

"""
    compatiblesSphere(plane, points, normals, eps, alpharad)

Create a bool-indexer array for those points that are compatible to the sphere.
Give back the projected points too for parameter space magic.

Compatibility is measured with an `eps` distance to the sphere and an `alpharad` angle to it's normal.
"""
function compatiblesSphere(sphere, points, normals, eps, alpharad)
    @assert length(points) == length(normals) "Size must be the same."
    # eps check
    o = sphere.center
    R = sphere.radius
    c1 = [abs(norm(a-o)-R) < eps for a in points]
    # alpha check
    α = cos(alpharad)
    if sphere.outwards
        c2 = [isparallel(normalize(points[i]-o), normals[i], α) && c1[i] for i in eachindex(points)]
    else
        c2 = [isparallel(normalize(o-points[i], normals[i]), α) && c1[i] for i in eachindex(points)]
    end

    under = Array{Int}(undef, 0)
    over = Array{Int}(undef, 0)
    for i in eachindex(points)
        if points[i][3] <= o[3]
            c2[i] && push!(under, i)
        else
            c2[i] && push!(over, i)
        end
    end
    # TODO
    # normalize by it's own length or maximum length??? r+3ϵ
    proj_points = [SVector{2}(normalize(a[1:2]-o[1:2])) for a in points]
    param_points = unitdisk2square.(proj_points)
    return c2, (under=under, over=over), param_points
end

end # module
