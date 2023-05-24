using Fourier2
using HDF5
using PyCall
using Test

DATA_PATH = joinpath(@__DIR__, "data.h5")
data = h5open(DATA_PATH, "r")
np = pyimport("numpy")

function read_arr(name)
    arr = read(data[name])
    return permutedims(arr, reverse(1:ndims(arr)))
end

function test_fft1_small()
    for i in 1:4
        arr = read_arr("small_1d_$i")
        arr_copy = copy(arr)

        res_np = np.fft.fft(arr)
        res_my = fft_pow2(arr)

        @test isequal(arr_copy, arr)
        @test isapprox(res_np, res_my)
    end
end

function test_fft1_large()
    for i in 1:2
        arr = read_arr("large_1d_$i")
        arr_copy = copy(arr)

        res_np = np.fft.fft(arr)
        res_my = fft_pow2(arr)

        @test isequal(arr_copy, arr)
        @test isapprox(res_np, res_my)
    end
end

function test_fft1_axis_small()
    for i in 1:4
        arr = read_arr("small_5d_$i")
        arr_copy = copy(arr)

        for axis in 1:5
            res_np = np.fft.fft(arr, axis=axis-1)
            res_my = fft_pow2(arr, axis)

            @test isequal(arr_copy, arr)
            @test isapprox(res_np, res_my)
        end
    end
end

function test_fft1_axis_large()
    for i in 1:2
        arr = read_arr("large_5d_$i")
        arr_copy = copy(arr)

        for axis in 1:5
            res_np = np.fft.fft(arr, axis=axis-1)
            res_my = fft_pow2(arr, axis)

            @test isequal(arr_copy, arr)
            @test isapprox(res_np, res_my)
        end
    end
end

function test_fft2_default_small()
    for i in 1:4
        arr = read_arr("small_2d_$i")
        arr_copy = copy(arr)

        res_np = np.fft.fft2(arr_copy)
        res_my = fft2_pow2(arr, algorithm="default")

        @test isequal(arr_copy, arr)
        @test isapprox(res_np, res_my)
    end
end

function test_fft2_default_large()
    for i in 1:2
        arr = read_arr("large_2d_$i")
        arr_copy = copy(arr)

        res_np = np.fft.fft2(arr_copy)
        res_my = fft2_pow2(arr, algorithm="default")

        @test isequal(arr_copy, arr)
        @test isapprox(res_np, res_my)
    end
end

function test_fft2_cooley_small()
    for i in 1:4
        arr = read_arr("small_2d_$i")
        arr_copy = copy(arr)

        res_np = np.fft.fft2(arr_copy)
        res_my = fft2_pow2(arr, algorithm="cooley")

        @test isequal(arr_copy, arr)
        @test isapprox(res_np, res_my)
    end
end

function test_fft2_cooley_large()
    for i in 1:2
        arr = read_arr("large_2d_$i")
        arr_copy = copy(arr)

        res_np = np.fft.fft2(arr_copy)
        res_my = fft2_pow2(arr, algorithm="cooley")

        @test isequal(arr_copy, arr)
        @test isapprox(res_np, res_my)
    end
end

@testset "Fourier.jl" begin

    test_fft1_small()
    test_fft1_large()

    test_fft1_axis_small()
    test_fft1_axis_large()

    test_fft2_default_small()
    test_fft2_default_large()

    test_fft2_cooley_small()
    test_fft2_cooley_large()

end # testset
