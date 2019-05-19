using Pkg
pkg"activate"

using Makie
using Colors
using GeometryTypes

ps = Point2f0[[-pi, -pi], [-pi, pi], [pi, pi], [pi, -pi]]
sc = Scene(resolution = (500, 500))
sc = poly!(sc, ps, color = :gray23)
axis = sc[Axis]
axis[:names, :axisnames] = ("q1 [rad]", "q2 [rad]")

Makie.save("whataname.png", sc)