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

export GridPoint, Cell, interpolate_bilinear, bilinear_weights, nearest, estimate_gradient

end
