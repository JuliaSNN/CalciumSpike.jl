using Test
using Random
using CalciumSpike

const CS = CalciumSpike

# noisy parameters so both the drift and the white noise are exercised
noisy = CaModel(η = 0.05, σ = 0.02)
spikes1 = Float32.(sort(rand(Xoshiro(1), 40) .* 9.9))
spikes_v = [Float32.(sort(rand(Xoshiro(10 + k), 30) .* 9.9)) for k in 1:5]
interval = (0.0, 10.0)
sr = 50.0

single(; kw...) = calcium_trace(spikes1, sr, interval; params = noisy, kw...)
multi(; kw...) = calcium_trace(spikes_v, sr, interval; params = noisy, kw...)

@testset "CalciumSpike seed" begin
    @testset "Vector{Float32}" begin
        Random.seed!(1); a = single(seed = 7)
        Random.seed!(999); rand(100); b = single(seed = 7)
        @test a.F == b.F && a.B == b.B && a.c == b.c       # bitwise
        c = single(seed = 8)
        @test c.F != a.F && c.B != a.B
        @test single(seed = Int32(7)).F == a.F
    end

    @testset "Vector{Vector{Float32}}" begin
        Random.seed!(1); a, ta = multi(seed = 7)
        Random.seed!(999); rand(100); b, tb = multi(seed = 7)
        @test a == b && ta == tb
        c, _ = multi(seed = 8)
        @test c != a
        @test all(a[k] != a[1] for k in 2:5)               # neurons draw different noise
        # common random numbers: perturbing neuron 3's spikes leaves other neurons' noise draws aligned
        pert = copy(spikes_v); pert[3] = Float32[]
        d, _ = calcium_trace(pert, sr, interval; params = noisy, seed = 7)
        @test d[1] == a[1] && d[5] == a[5] && d[3] != a[3]
    end

    @testset "default stays random (global RNG)" begin
        Random.seed!(1); a = single().F
        Random.seed!(1); b = single().F
        @test a == b
        @test single().F != single().F
        Random.seed!(1); a, _ = multi()
        Random.seed!(1); b, _ = multi()
        @test a == b
        @test multi()[1] != multi()[1]
        # seeded calls leave the global stream untouched
        Random.seed!(5); x = rand(); Random.seed!(5); single(seed = 3); @test rand() == x
        # explicit rng behaves like before
        @test single(rng = Xoshiro(2)).F == single(rng = Xoshiro(2)).F
    end

    @testset "rng and seed are exclusive" begin
        @test_throws ArgumentError single(rng = Xoshiro(1), seed = 1)
        @test_throws ArgumentError multi(rng = Xoshiro(1), seed = 1)
        @test_throws ArgumentError single(seed = 1.5)
    end

    @testset "noiseless model is deterministic" begin
        p = CaModel()
        @test calcium_trace(spikes1, sr, interval; params = p).F ==
              calcium_trace(spikes1, sr, interval; params = p).F
    end

    @testset "generate_synthetic_data" begin
        Random.seed!(1); a = CS.generate_synthetic_data(4, 6, 3; num_signals = 2, seed = 5)
        Random.seed!(2); b = CS.generate_synthetic_data(4, 6, 3; num_signals = 2, seed = 5)
        @test a == b
        Random.seed!(1); x = CS.generate_synthetic_data(4, 6, 3; num_signals = 2)
        Random.seed!(1); y = CS.generate_synthetic_data(4, 6, 3; num_signals = 2)
        @test x == y
        @test_throws ArgumentError CS.generate_synthetic_data(4, 6, 3; rng = Xoshiro(1), seed = 1)
    end
end
