module GridHelpers

struct GridPoint
  coords::Tuple{Int,Int}
end

Base.getindex(point::GridPoint) = point.coords
Base.iterate(point::GridPoint, args...) = iterate(point[], args...)
Base.length(point::GridPoint) = length(point[])
Base.convert(::Type{GridPoint}, coords::Tuple{Int,Int}) = GridPoint(coords)

function Base.getproperty(point::GridPoint, name::Symbol)
  name === :bottom && return point .+ (0, -1)
  name === :top && return point .+ (0, 1)
  name === :left && return point .+ (0, -1)
  name === :right && return point .+ (0, 1)
  getfield(point, name)
end

"""
Cell defined by four corner points.
"""
struct Cell
  bottom_left::GridPoint
  bottom_right::GridPoint
  top_right::GridPoint
  top_left::GridPoint
end

Base.length(::Cell) = 4
Base.iterate(cell::Cell, args...) = iterate((cell.bottom_left, cell.bottom_right, cell.top_left, cell.top_right), args...)
Base.getindex(cell::Cell, i) = iterate(cell, i)[1]

Base.getindex(A::AbstractArray, point::GridPoint) = A[point[]...]
Base.setindex!(A::AbstractArray, value, point::GridPoint) = A[point[]...] = value

nearest((x, y)) = Int.(round.((x, y)))

Cell(x, y) = Cell((x, y))
function Cell((x, y))
  (i, j) = Int.(floor.((x, y)))
  Cell((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
end

function bilinear_weights(cell::Cell, (x, y))
  (cx, cy) = cell.bottom_left
  w1 = (-(x - (1 + cx)) * -(y - (1 + cy)))
  w2 = ((x - cx) * -(y - (1 + cy)))
  w3 = (-(x - (1 + cx)) * (y - cy))
  w4 = ((x - cx) * (y - cy))
  (w1, w2, w3, w4)
end

function interpolate_bilinear(A, position)
  cell = Cell(position)
  weights = bilinear_weights(cell, position)
  sum(i -> A[cell[i]] * weights[i], 1:4; init = zero(eltype(A)))
end

export GridPoint, Cell, interpolate_bilinear, bilinear_weights

end
