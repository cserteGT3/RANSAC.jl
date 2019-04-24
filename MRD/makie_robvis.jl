using Pkg
pkg"activate"

using StaticArrays: SVector
using Makie

"""
    drawcf(s, o, x, y, z; llength = 0.5, lwidth = 0.5)

Draw a coordinate frame. XYZ is RGB.

# Arguments:
- `s::Scene`: Makie scene.
- `o::AbstractArray`: origo of coordinate frame.
- `x::AbstractArray`: x axis of coordinate frame.
- `y::AbstractArray`: y axis of coordinate frame.
- `z::AbstractArray`: z axis of coordinate frame.
- `line_size = 0.5`: length of the lines.
- `lwidth = 0.5`: width of the lines.
"""
function drawcf(s, o, x, y, z; llength = 0.5, lwidth = 0.5)
    ox = [o => o+llength*x]
    oy = [o => o+llength*y]
    oz = [o => o+llength*z]
    linesegments!(s, ox, color = :red, linewidth = lwidth)
    linesegments!(s, oy, color = :green, linewidth = lwidth)
    linesegments!(s, oz, color = :blue, linewidth = lwidth)
    s
end

"""

Only ring tori.
"""
function tori(R, r, lR, lr)
    @assert R > r "Ring tori, I told you."
    x(Θ, ϕ) = (R + r*cos(Θ))*cos(ϕ)
    y(Θ, ϕ) = (R + r*cos(Θ))*sin(ϕ)
    z(Θ, ϕ) = r*sin(Θ)
    ps = [SVector(x(i,j), y(i,j), z(i,j)) for i in range(0, 2π, length=lR) for j in range(0, 2π, length=lr)]
end

function array2arrays(A)
    xout = Array{eltype(A[1])}(undef, length(A))
    yout = Array{eltype(A[1])}(undef, length(A))
    zout = Array{eltype(A[1])}(undef, length(A))
    for i in eachindex(A)
        xout[i] = A[i][1]
        yout[i] = A[i][2]
        zout[i] = A[i][3]
    end
    xout, yout, zout
end

o = SVector(0.0,0,0);
ax = SVector(1,0,0.0);
ay = SVector(0,1,0.0);
az = SVector(0,0,1.0);
sri = [SVector{3}(rand(3)) for _ in 1:100];
stori = tori(10,5, 100, 50);
x,y,z = array2arrays(stori);
meshscatter(x, y, z)

s = Scene();
scatter!(s, sri)
drawcf(s, o, ax, ay, az, lwidth = 5)

function makietori(R, r, lR, lr)
    @assert R > r "Ring tori, I told you."
    theta = range(0, 2π, length=lR)
    phi = range(0, 2π, length=lr)
    x_(Θ, ϕ) = (R + r*cos(Θ))*cos(ϕ)
    y_(Θ, ϕ) = (R + r*cos(Θ))*sin(ϕ)
    z_(Θ, ϕ) = r*sin(Θ)
    x = [x_(θ, ϕ) for θ in theta, ϕ in phi]
    y = [y_(θ, ϕ) for θ in theta, ϕ in phi]
    z = [z_(θ, ϕ) for θ in theta, ϕ in phi]
    return x, y, z
end

mx, my, mz = makietori(10,2, 500, 500);
surface(mx, my, mz, colormap=[:green], transparency = false, shading = true)
