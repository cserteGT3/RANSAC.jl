using Pkg
Pkg.activate()

# every include
using LinearAlgebra
using StaticArrays
using Random
using Revise
using Colors
using Makie

includet("../RANSAC.jl")
using .RANSAC

function just_time(ss, itlength)
    # ss = [1,4,8,16]
    # itlength = (30000, 30000)

    vsp, nsp, _, _ = examplepc3()
    vsn, nsn, _, _ = examplepc3(true, all = true, mrotdeg = 5, vertscale = 0.3)

    # Pure
    p_aep = (ϵ = 0.3, α=deg2rad(5));
    cy_aep = (ϵ = 0.3, α=deg2rad(5));
    sp_aep = (ϵ = 0.3, α=deg2rad(5));
    p_ae = AlphSilon(sp_aep, p_aep, cy_aep);

    # noised
    p_aen = (ϵ = 0.3, α=deg2rad(5));
    cy_aen = (ϵ = 0.3, α=deg2rad(5));
    sp_aen = (ϵ = 0.3, α=deg2rad(5));
    n_ae = AlphSilon(sp_aen, p_aen, cy_aen);

    resu = zeros(length(ss),2)
    for i in 1:length(ss)
        tri = ss[i]
        # pure
        pc_p = PointCloud(vsp, nsp, tri)
        resu[i, 1] = @elapsed ransac(pc_p, p_ae, 15, 0.9, 900, itlength[1], 3, 500, true)

        # noised
        pc_n = PointCloud(vsn, nsn, tri)
        resu[i, 2] = @elapsed ransac(pc_n, n_ae, 15, 0.9, 900, itlength[2], 3, 500, true)
    end
    resu
end

function just_result(ss, itlength)
    # ss = [1,4,8,16]
    # itlength = (30000, 30000)

    vsp, nsp, _, _ = examplepc3()
    vsn, nsn, _, _ = examplepc3(true, all = true, mrotdeg = 5, vertscale = 0.3)

    # Pure
    p_aep = (ϵ = 0.3, α=deg2rad(5));
    cy_aep = (ϵ = 0.3, α=deg2rad(5));
    sp_aep = (ϵ = 0.3, α=deg2rad(5));
    p_ae = AlphSilon(sp_aep, p_aep, cy_aep);

    # noised
    p_aen = (ϵ = 0.3, α=deg2rad(5));
    cy_aen = (ϵ = 0.3, α=deg2rad(5));
    sp_aen = (ϵ = 0.3, α=deg2rad(5));
    n_ae = AlphSilon(sp_aen, p_aen, cy_aen);
    resu1 = [ransac(PointCloud(vsp, nsp, tri), p_ae, 15, 0.9, 900, itlength[1], 3, 500, true) for tri in ss]
    resu2 = [ransac(PointCloud(vsn, nsn, tri), p_ae, 15, 0.9, 900, itlength[2], 3, 500, true) for tri in ss]
    resu1, resu2
end

function just_time1(ss, itlength)
    # ss = [1,4,8,16]
    # itlength = (30000, 30000)

    vsp, nsp, _, _ = examplepc3()

    # Pure
    p_aep = (ϵ = 0.3, α=deg2rad(5));
    cy_aep = (ϵ = 0.3, α=deg2rad(5));
    sp_aep = (ϵ = 0.3, α=deg2rad(5));
    p_ae = AlphSilon(sp_aep, p_aep, cy_aep);

    resu = zeros(length(ss))
    for i in 1:length(ss)
        tri = ss[i]
        # pure
        pc_p = PointCloud(vsp, nsp, tri)
        Random.seed!(1234)
        resu[i, 1] = @elapsed ransac(pc_p, p_ae, 15, 0.9, 900, itlength, 3, 500, true)
    end
    resu
end

subst = range(1,step=4,length=2)
substt = just_time1(subst, 20)

lsc = lines(subst,substt)

ax = lsc[Axis]
ax[:names, :axisnames] = ("subsetek száma", "idő [s]")
Makie.save("timings.png", lsc )
