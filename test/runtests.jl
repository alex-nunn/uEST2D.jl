using uEST2D
using Test

@testset "uEST2D.jl" begin
    @testset "ion_trajectory" begin
        @testset "trajectory success" begin
            traj = ion_trajectory(1.4, 0.4, 0.00, 20.0)
            @test Int(traj.retcode) == 1
        end

        @testset "trajectory terminated" begin
            traj = ion_trajectory(1.4, 0.4, 0.1, 20.0)
            @test Int(traj.retcode) == 2
        end
    end
end