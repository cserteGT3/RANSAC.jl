# Generate code for the filters

using ControlSystems
using LinearAlgebra
using Logging

plotq = @isdefined plotyeah
if !plotq
	@info "Set `plotyeah=true` if you want to bodeplot()."
end

"""
    cSS(A, B)

Create a dumb state-space for `c2d()` use.
"""
function cSS(A, B)
    C = zeros(size(B))
    C[end] = 1
    ss(A,B,C',0)
end

function makeLyrics(sys)
    lyrics = Array{String}(undef, 0)
    it = 1:size(sys.A, 1)
    for i in it, j in it
        push!(lyrics, "float ad$i$j = $(sys.A[i,j]);")
    end
    for i in it
        push!(lyrics, "float bd$i = $(sys.B[i]);")
    end
    lyrics
end

function writeLyrics(lyricsname, lyrics)
    open(lyricsname, "w") do io
        for a in lyrics
            write(io, a*"\n")
        end
    end
end

function ploti(ss, n)
	plotq = @isdefined plotyeah
	if plotq
		if plotyeah
			pl = bodeplot(ss)
			Plots.savefig(pl, n)
		end
	end
end

"""
    rad2enc(x)

Compute radian into encoder ticks.
"""
rad2enc(x) = x*512/(2*π)

## Constants

Ts = 0.001;
ω_zaj = π/Ts;
Tc = 10/ω_zaj;

## Low pass filter

lpfA = [0 1 0; 0 0 1; -Tc^-3 -3/Tc^2 -3/Tc]
lpfB = [0, 0, Tc^-3];

Dlpf, _ = c2d(cSS(lpfA, lpfB), Ts);
writeLyrics("lpf.txt", makeLyrics(Dlpf));

ploti(Dlpf, "lpf.pdf")

## Harmad Bessel

b3A = [0 1 0; 0 0 1; -15/Tc^3 -15/Tc^2 -6/Tc]
b3B = [0, 0, 15/Tc^3];

Dbessel3, _ = c2d(cSS(b3A, b3B), Ts);
writeLyrics("bessel3.txt", makeLyrics(Dbessel3));

ploti(Dbessel3, "bessel3.pdf")

## Otod Bessel

b5A = vcat(hcat(zeros(4), I), [-945/Tc^5 -945/Tc^4 -420/Tc^3 -105/Tc^2 -15/Tc]);
b5B = [0, 0, 0, 0, 945/Tc^5];

Dbessel5, _ = c2d(cSS(b5A, b5B), Ts);
writeLyrics("bessel5.txt", makeLyrics(Dbessel5));

ploti(Dbessel5, "bessel5.pdf")
