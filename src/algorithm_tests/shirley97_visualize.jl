# Visualization of shirley1997 - A low distortion map between disk and square

using LinearAlgebra
using StaticArrays: SVector
using Makie

include("../RANSAC.jl")
using .RANSAC

struct Circ{T<:Real}
    R::T
    db::Int
end

function concentriccircles(c)
    @assert 0 < c.R && c.R <= 1 "Radius must be: 0<R<=1 - unit disk."
    @assert c.db > 0 "How the *** should I sample $(c.db) samples?!"
    ϕ = 2*π/c.db
    [SVector(c.R*cos(k*ϕ), c.R*sin(k*ϕ)) for k in 0:c.db-1]
end

c1 = concentriccircles(Circ(0.4, 4));
c2 = concentriccircles(Circ(0.8, 12));
unitc = concentriccircles(Circ(1, 99));

s1 = Scene(resolution = (500,500));
scatter!(c1, color = :red);
scatter!(s1, c2, color = :blue);
lines!(s1, unitc, color = :black);

sq1 = unitdisk2square.(c1);
sq2 = unitdisk2square.(c2);
#unitsq = [SVector(i, j) for i in -1:2:1 for j in -1:2:1];
unitsq = [SVector(i, j) for i in 0:1 for j in 0:1];
unitsq[3:4] = unitsq[4:-1:3];
push!(unitsq, unitsq[1]);

s2 = Scene(resolution = (500,500));
scatter!(s2, sq1, color = :red);
scatter!(s2, sq2, color = :blue);
lines!(s2, unitsq, color = :black);

s = vbox(s1, s2)
