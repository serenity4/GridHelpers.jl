module GridHelpers

using StaticArrays

"""
Point on a grid.

This is not a cell, which represents a whole face and not a vertex.
"""
struct GridPoint{T<:Integer}
  coords::SVector{2,T}
end
GridPoint(coords::Tuple) = GridPoint(SVector(coords))
GridPoint(i, j) = GridPoint((i, j))

Base.convert(::Type{GridPoint{T}}, coords::Tuple{T, T}) where {T} = GridPoint(coords)

Base.getindex(point::GridPoint) = point.coords
Base.getindex(point::GridPoint, i) = getindex(point[], i)
Base.iterate(point::GridPoint, args...) = iterate(point[], args...)
Base.length(point::GridPoint) = length(point[])
Base.convert(::Type{GridPoint}, coords::Tuple{Int,Int}) = GridPoint(coords)

function Base.getproperty(point::GridPoint{T}, name::Symbol) where {T}
  name === :left && return GridPoint(point[] .- @SVector T[1, 0])
  name === :right && return GridPoint(point[] .+ @SVector T[1, 0])
  name === :bottom && return GridPoint(point[] .- @SVector T[0, 1])
  name === :top && return GridPoint(point[] .+ @SVector T[0, 1])
  getfield(point, name)
end

Base.getindex(A::AbstractArray, point::GridPoint) = A[point[]...]
Base.getindex(A, point::GridPoint) = A[point[]...]
Base.setindex!(A::AbstractArray, value, point::GridPoint) = A[point[]...] = value
Base.setindex!(A, value, point::GridPoint) = A[point[]...] = value

neighbor(point::GridPoint{T}, i) where {T} = i == T(1) ? point.left : i == T(2) ? point.right : i == T(3) ? point.bottom : point.top
is_outside_grid(point, (ni, nj)) = !is_inside_grid(point, (ni, nj))
is_inside_grid((x, y), (ni, nj)) = one(x) ≤ x ≤ ni && one(y) ≤ y ≤ nj
materialize_grid((ni, nj)) = [GridPoint(i, j) for i in one(ni):ni, j in one(nj):nj]
materialize_grid(A::AbstractMatrix) = materialize_grid(Int, size(A))

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
struct Cell{T}
  bottom_left::GridPoint{T}
  bottom_right::GridPoint{T}
  top_right::GridPoint{T}
  top_left::GridPoint{T}
end

Base.length(::Cell{T}) where {T} = T(4)
Base.iterate(cell::Cell, args...) = iterate((cell.bottom_left, cell.bottom_right, cell.top_left, cell.top_right), args...)
Base.getindex(cell::Cell{T}, i) where {T} = iterate(cell, i)[one(T)]

nearest((x, y)) = Int.(round.((x, y)))

Cell(x::Number, y::Number) = Cell((x, y))
Cell(position, bounds) = Cell{Int}(position, bounds)
Cell{T}(position, bounds) where {T} = Cell{T}(ifelse.(position .≥ bounds, position .- one(T), position))
Cell((x, y)::Tuple) = Cell{Int}((x, y))
function Cell{T}((x, y)) where {T}
  (i, j) = floor.(T, (x, y))
  Cell{T}((i, j), (i + one(T), j), (i + one(T), j + one(T)), (i, j + one(T)))
end
Cell((x, y)) = Cell{Int}((x, y))

"""
Return bilinear weights from all four corners, in cell iteration order.
Extracting bilinear weights with this function may be useful when they are to be combined with other weights.

!!! note
    If you want to interpolate a value directly, use `interpolate_bilinear` which will be much faster.
"""
function bilinear_weights(cell::Cell{T}, (x, y)) where {T}
  (cx, cy) = cell.bottom_left
  w1 = (-(x - (one(T) + cx)) * -(y - (one(T) + cy)))
  w2 = ((x - cx) * -(y - (one(T) + cy)))
  w3 = (-(x - (one(T) + cx)) * (y - cy))
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

lerp(x, y, w) = x * (one(w) - w) + y * w

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
