using Pkg
pkg"activate"

using StaticArrays: SVector
using Makie
using Colors

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

function makietori(R, r, lR, lr; z0 = 0)
    @assert R > r "Ring tori, I told you."
    theta = range(0, 2π, length=lR)
    phi = range(0, 2π, length=lr)
    x_(Θ, ϕ) = (R + r*cos(Θ))*cos(ϕ)
    y_(Θ, ϕ) = (R + r*cos(Θ))*sin(ϕ)
    z_(Θ, ϕ) = z0 + r*sin(Θ)
    x = [x_(θ, ϕ) for θ in theta, ϕ in phi]
    y = [y_(θ, ϕ) for θ in theta, ϕ in phi]
    z = [z_(θ, ϕ) for θ in theta, ϕ in phi]
    return x, y, z
end

mx, my, mz = makietori(10,2, 500, 500);
surface(mx, my, mz, colormap=[:green], transparency = false, shading = true)

## MRD homework parameters
l10 = 0.8; #m
l11 = 0.6; #m
l20 = 0.3; #m

o = SVector(0,0,0.0);
v1 = SVector(0,0,l10);
v11 = SVector(0,l11,0);
v20 = SVector(0,0,l20);

s1 = o+v1;
s2 = s1+v11;
s3 = s2+v20;

ax = SVector(1,0,0.0);
ay = SVector(0,1,0.0);
az = SVector(0,0,1.0);

segments = [o => s1; s1 => s2; s2 => s3];

n = 100;
mx, my, mz = makietori(l11, l20, n, n, z0 = l10);
s = surface(mx, my, mz, color = fill(RGBA(1.,1.,1.,0.3), n, n), transparency = true, shading = true)
linesegments!(s, segments, linewidth = 5)
axis = s[Axis];
axis[:showgrid] = false

# x0: o, ax,ay,az
# x1: s1, ay,-ax,az
# x2: s2, ay,-az,-ax
# x3: s3, ax, ay, az
ll = 0.15;
lw = 4;
drawcf(s, o, ax, ay, az, llength = ll, lwidth = lw)
drawcf(s, s1, ay, -ax, az, llength = ll, lwidth = lw)
drawcf(s, s2, ay, -az, -ax, llength = ll, lwidth = lw)
drawcf(s, s3, ax, ay, az, llength = ll, lwidth = lw)

Makie.save("angle2.png", s)
