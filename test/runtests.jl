using LinearAlgebra
using Random
using CCEQR
using Test

function showInfo(msg, result)
    if typeof(result) != Test.Pass
        @info msg
    end
end

@testset "test tallthin gaussian" begin
    rng = MersenneTwister(1)

    for k = 1:20
        for rho in exp10.(range(-10, -.01, 10))
            pstr = "parameters are k = "*string(k)*", rho = "*string(rho)
            
            for t = 1:20
                A          = randn(rng, 5000, 20)
                p0         = qr(A, ColumnNorm()).p
                p, _, _, _ = cceqr!(deepcopy(A), k = k, rho = rho)
                showInfo(pstr, @test (p[1:k] == p0[1:k]))
            end
        end
    end
end

@testset "test shortwide gaussian" begin
    rng = MersenneTwister(1)

    for k = 1:20
        for rho in exp10.(range(-10, -.01, 10))
            pstr = "parameters are k = "*string(k)*", rho = "*string(rho)
            
            for t = 1:20
                A          = randn(rng, 20, 5000)
                p0         = qr(A, ColumnNorm()).p
                p, _, _, _ = cceqr!(deepcopy(A), k = k, rho = rho)
                showInfo(pstr, @test (p[1:k] == p0[1:k]))
            end
        end
    end
end

@testset "test shortwide spiked" begin
    rng = MersenneTwister(1)

    for k = 1:20
        for rho in exp10.(range(-10, -.01, 10))
            pstr = "parameters are k = "*string(k)*", rho = "*string(rho)
            
            for t = 1:20
                A = randn(rng, 20, 5000)
                d = (.99).^(1:5000)
                d = d[randperm(rng, 5000)]
                rmul!(A, Diagonal(d))

                p0         = qr(A, ColumnNorm()).p
                p, _, _, _ = cceqr!(deepcopy(A), k = k, rho = rho)
                showInfo(pstr, @test (p[1:k] == p0[1:k]))
            end
        end
    end
end

@testset "test full qr" begin
    rng = MersenneTwister(1)
    
    for k = 1:20
        for t = 1:20
            A = randn(rng, 20, 5000)
            B = randn(rng, 5000, 20)

            cceqr!(A, k = k, full = true)
            cceqr!(B, k = k, full = true)
        end
    end
end
