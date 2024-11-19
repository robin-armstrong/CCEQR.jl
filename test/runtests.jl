using LinearAlgebra
using Random
using CCEQR
using Test

function showInfo(msg, result)
    if typeof(result) != Test.Pass
        @info msg
    end
end

@testset "test shortwide gaussian" begin
    rng = MersenneTwister(1)
    A   = randn(rng, 20, 5000)
    p0  = qr(A, ColumnNorm()).p

    for k = 1:20
        for rho in exp10.(range(-10, -.1, 10))
            for full in [true, false]
                pstr = "parameters are: k = "*string(k)*", rho = "*string(rho)*", full = "*string(full)
                
                p, _, _, _ = cceqr!(deepcopy(A), k = k, rho = rho, full = full)
                showInfo(pstr, @test (p[1:k] == p0[1:k]))
            end
        end
    end
end

@testset "test tallthin gaussian" begin
    rng = MersenneTwister(1)
    A   = randn(rng, 5000, 20)
    p0  = qr(A, ColumnNorm()).p

    for k = 1:20
        for rho in exp10.(range(-10, -.1, 10))
            for full in [true, false]
                pstr = "parameters are: k = "*string(k)*", rho = "*string(rho)*", full = "*string(full)
                
                p, _, _, _ = cceqr!(deepcopy(A), k = k, rho = rho, full = full)
                showInfo(pstr, @test (p[1:k] == p0[1:k]))
            end
        end
    end
end
