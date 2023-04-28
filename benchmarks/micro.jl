using GridHelpers
using BenchmarkTools

A = rand(512, 512)
position = (6.3, 421.7)
cell = Cell(position)

@btime interpolate_bilinear($A, $position)
@btime interpolate_bilinear($A, $position, $cell)
@btime estimate_gradient($A, $position)
@btime estimate_gradient($A, $position, $cell)
@btime bilinear_weights($cell, $position)
