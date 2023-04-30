using GridHelpers
using Test

@testset "GridHelpers.jl" begin
    @testset "Grid points" begin
        point = GridPoint((5, 5))
        @test collect(point) == [5, 5]
        @test point[] == (5, 5)
        @test point.left === GridPoint((4, 5))
        @test point.right === GridPoint((6, 5))
        @test point.bottom === GridPoint((5, 4))
        @test point.top === GridPoint((5, 6))
        @test neighbor(point, 3) == point.bottom
        grid_size = (5, 6)
        @test !is_outside_grid(GridPoint((1, 1)), grid_size)
        @test is_outside_grid(GridPoint((1, 0)), grid_size)
        @test !is_outside_grid(GridPoint(grid_size), grid_size)
        @test is_outside_grid(GridPoint((6, 6)), grid_size)
    end

    @testset "Cell" begin
        position = (20.0, 30.0)
        cell = Cell(position)
        @test cell.bottom_left == GridPoint(position)
        @test cell.bottom_right == GridPoint(position .+ (1, 0))
        @test cell.top_left == GridPoint(position .+ (0, 1))
        @test cell.top_right == GridPoint(position .+ 1)
        @test collect(cell) == [cell[1], cell[2], cell[3], cell[4]]
        weights = bilinear_weights(cell, position)
        @test weights[1] == 1.0
        @test sum(weights) == 1
        @test nearest((4.3, 6.9)) == (4, 7)
        @test all(≥(0), weights)
        weights = bilinear_weights(cell, position .+ 0.5)
        @test all(==(0.25), weights)
        weights = bilinear_weights(cell, position .+ (0.13, 0.78))
        @test sum(weights) == 1
        @test all(≥(0), weights)
        @test all(bilinear_weights(cell, cell[i])[i] == 1.0 for i in 1:4)
        A = zeros(512, 512)
        @test interpolate_bilinear(A, (4.3, 6.1)) == 0.0
        A = ones(512, 512)
        @test interpolate_bilinear(A, (4.3, 6.1)) == 1.0
        A = rand(512, 512)
        @test 0 ≤ interpolate_bilinear(A, (4.3, 6.1)) ≤ 1.0
        @test estimate_gradient(A, (4.3, 6.1)) isa Tuple{Float64,Float64}
        @test estimate_gradient(A, GridPoint(4, 6), size(A)) isa Tuple{Float64,Float64}

        cell = Cell((30, 40))
        @test cell == Cell((30.1, 40.1))
        cell = Cell((30, 40), (30, 40))
        @test cell == Cell((29, 39))
    end
end;
