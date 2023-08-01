module GridHelpers

using StaticArrays

"""
Point on a grid.

This is not a cell, which represents a whole face and not a vertex.
"""
struct GridPoint
  coords::Tuple{Int,Int}
end
GridPoint(i, j) = GridPoint((i, j))

Base.getindex(point::GridPoint) = point.coords
Base.getindex(point::GridPoint, i) = getindex(point[], i)
Base.iterate(point::GridPoint, args...) = iterate(point[], args...)
Base.length(point::GridPoint) = length(point[])
Base.convert(::Type{GridPoint}, coords::Tuple{Int,Int}) = GridPoint(coords)

function Base.getproperty(point::GridPoint, name::Symbol)
  name === :left && return GridPoint(point[] .+ (-1, 0))
  name === :right && return GridPoint(point[] .+ (1, 0))
  name === :bottom && return GridPoint(point[] .+ (0, -1))
  name === :top && return GridPoint(point[] .+ (0, 1))
  getfield(point, name)
end

neighbor(point::GridPoint, i) = i == 1 ? point.left : i == 2 ? point.right : i == 3 ? point.bottom : point.top
is_outside_grid(point, (ni, nj)) = !is_inside_grid(point, (ni, nj))
is_inside_grid((x, y), (ni, nj)) = 1 ≤ x ≤ ni && 1 ≤ y ≤ nj
materialize_grid((ni, nj)) = [GridPoint(i, j) for i in 1:ni, j in 1:nj]
materialize_grid(A::AbstractMatrix) = materialize_grid(size(A))

# Centered finite-difference method with a spatial step of 1.
function estimate_gradient(A, point::GridPoint, (ni, nj))
  sx = point[1] in (1, ni) ? 0.0 : lerp(A[point.left], A[point], 0.5) - lerp(A[point], A[point.right], 0.5)
  sy = point[2] in (1, nj) ? 0.0 : lerp(A[point.bottom], A[point], 0.5) - lerp(A[point], A[point.top], 0.5)
  (sx, sy)
end

abstract type Neighborhood end
"Four closest points, directly adjacent to the current grid point."
struct FourNeighbors <: Neighborhood end
"Eight closest points, including the four closest plus the slightly farthest grid point corners."
struct EightNeighbors <: Neighborhood end

neighbors(point::GridPoint) = neighbors(point, FourNeighbors())
neighbors(point::GridPoint, ::FourNeighbors) = @SVector [point.left, point.right, point.bottom, point.top]
neighbors(point::GridPoint, ::EightNeighbors) = @SVector [point.top.left, point.left, point.bottom.left, point.bottom, point.bottom.right, point.right, point.top.right, point.top]

"""
Cell defined by four corner points around coordinates `(x, y)`. `(x, y)` should be floating-point indices
associated to 1-based arrays, i.e. the rounded value of `x` and `y` should be an integer-valued grid index.

Iteration order is defined as follows: bottom left, bottom right, top left, top right.
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

Cell(x::Number, y::Number) = Cell((x, y))
Cell(position, bounds) = Cell(ifelse.(position .≥ bounds, position .- 1, position))
function Cell((x, y))
  (i, j) = floor.(Int, (x, y))
  Cell((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1))
end

"""
Return bilinear weights from all four corners, in cell iteration order.
Extracting bilinear weights with this function may be useful when they are to be combined with other weights.

!!! note
    If you want to interpolate a value directly, use `interpolate_bilinear` which will be much faster.
"""
function bilinear_weights(cell::Cell, (x, y))
  (cx, cy) = cell.bottom_left
  w1 = (-(x - (1 + cx)) * -(y - (1 + cy)))
  w2 = ((x - cx) * -(y - (1 + cy)))
  w3 = (-(x - (1 + cx)) * (y - cy))
  w4 = ((x - cx) * (y - cy))
  (w1, w2, w3, w4)
end

"""
Perform a bilinear interpolation of `A` at location `(x, y)`.

If the cell is provided, 
"""
function interpolate_bilinear(A, location, cell::Cell = Cell(location))
  (x, y) = location
  cx, cy = cell.bottom_left
  nx0 = lerp(A[cell.bottom_left], A[cell.bottom_right], x - cx)
  nx1 = lerp(A[cell.top_left], A[cell.top_right], x - cx)
  lerp(nx0, nx1, y - cy)
end

lerp(x, y, w) = x * (1 - w) + y * w

function estimate_gradient(A, location, cell::Cell = Cell(location))
  (x, y) = location
  cx, cy = cell.bottom_left
  gx = lerp(A[cell.bottom_right] - A[cell.bottom_left], A[cell.top_right] - A[cell.top_left], x - cx)
  gy = lerp(A[cell.top_right] - A[cell.bottom_right], A[cell.top_left] - A[cell.bottom_left], y - cy)
  (gx, gy)
end

export
  GridPoint,
  FourNeighbors, EightNeighbors, neighbors, neighbor,
  is_inside_grid, is_outside_grid,
  Cell,
  interpolate_bilinear,
  bilinear_weights,
  nearest,
  estimate_gradient,
  materialize_grid

end
