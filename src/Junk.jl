"""
(-p)-th n-root of unity
"""
function W(n, p)
    return exp(-2 * 1im * π * p / n)
end


function _fft2_pow2_cooley_square_req!(arr)
    n = size(arr, 1)
    if n <= 2
        if n == 1
            return
        end
        tmp = [
            arr[1, 1] + arr[1, 2] + arr[2, 1] + arr[2, 2] arr[1, 1] - arr[1, 2] + arr[2, 1] - arr[2, 2];
            arr[1, 1] + arr[1, 2] - arr[2, 1] - arr[2, 2] arr[1, 1] - arr[1, 2] - arr[2, 1] + arr[2, 2]
        ]
        copy!(arr, tmp)
        return
    end

    new_n = n >> 1
    tmp = fill(Matrix{ComplexF64}(undef, new_n, new_n), 4)
    tmp[1] = arr[1:2:n, 1:2:n]
    tmp[2] = arr[2:2:n, 1:2:n]
    tmp[3] = arr[1:2:n, 2:2:n]
    tmp[4] = arr[2:2:n, 2:2:n]

    _fft2_pow2_cooley_square_req!.(tmp[1:4])

    println.(tmp)

    wnp = W.(n, 0:new_n-1)
    tmp[2] .*= wnp
    tmp[3] .*= transpose(wnp)
    tmp[4] .*= transpose(wnp)
    tmp[4] .*= wnp

    arr[1:new_n, 1:new_n] = tmp[1] + tmp[2] + tmp[3] + tmp[4]

    arr[new_n+1:n, 1:new_n] = tmp[1] - tmp[2] + tmp[3] - tmp[4]

    arr[1:new_n, new_n+1:n] = tmp[1] + tmp[2] - tmp[3] - tmp[4]

    arr[new_n+1:n, new_n+1:n] = tmp[1] - tmp[2] - tmp[3] + tmp[4]

    nothing
end

    # n = UInt64(size(arr, 1))
    # len = Int(round(log2int(n)))
    # k = (k == 0) ? len : k

    # for i in 2:n
    #     j = _rev_k_bits_and_move(i - 1, len, k) + 1
    #     if i < j
    #         arr[i], arr[j] = arr[j], arr[i]
    #     end
    # end

    # arr[block1 + 1, block2 + 1] += arr[block1 + 2, block2 + 1] + arr[block1 + 2, block2 + 2] + arr[block1 + 1, block2 + 2]
    # u12 = arr[block1 + 1, block2 + 2]
    # u21 = arr[block1 + 2, block2 + 1]
    # arr[block1 + 1, block2 + 2] = arr[block1 + 1, block2 + 1] - 2 * (arr[block1 + 1, block2 + 2] + arr[block1 + 2, block2 + 2])
    # arr[block1 + 2, block2 + 1] = arr[block1 + 1, block2 + 1] - 2 * (arr[block1 + 2, block2 + 1] + arr[block1 + 2, block2 + 2])
    # arr[block1 + 2, block2 + 2] = arr[block1 + 1, block2 + 1] - 2 * (u12 + u21)

        # m, n = size(arr)
    # j = 1
    # for i in 2:(axis==1 ? n : m)
    #     bit = n >> 1
    #     while j > bit
    #         j -= bit
    #         bit >>= 1
    #     end
    #     j += bit
    #     if i < j
    #         for k in 1:(axis==1 ? m : n)
    #             if axis == 1
    #                 arr[k, i], arr[k, j] = arr[k, j], arr[k, i]
    #             else
    #                 arr[i, k], arr[j, k] = arr[j, k], arr[i, k]
    #             end
    #         end
    #     end
    # end

    # n, m = size(arr)
    # if n == m
    #     _fft2_pow2_cooley_square!(arr)
    #     return
    # end

    # new_m = m >> 1
    # tmp = fill(Matrix{ComplexF64}(undef, n, new_m), 2)

    # tmp[1] = arr[:, 1:2:m]
    # tmp[2] = arr[:, 2:2:m]

    # _fft2_pow2_cooley_rect1!.(tmp[1:2])

    # wnp = W.(m, 0:new_m-1)
    # println(size(wnp))
    # println(size(tmp[2]))
    # tmp[2] .*= transpose(wnp)

    # arr[:, 1:new_m] = tmp[1] + tmp[2]
    # arr[:, new_m+1:m] = tmp[1] - tmp[2]

    # arr[block1 + 1, block2 + 1] += arr[block1 + 2, block2 + 1] + arr[block1 + 2, block2 + 2] + arr[block1 + 1, block2 + 2]
    # u12 = arr[block1 + 1, block2 + 2]
    # u21 = arr[block1 + 2, block2 + 1]
    # arr[block1 + 1, block2 + 2] = arr[block1 + 1, block2 + 1] - 2 * (arr[block1 + 1, block2 + 2] + arr[block1 + 2, block2 + 2])
    # arr[block1 + 2, block2 + 1] = arr[block1 + 1, block2 + 1] - 2 * (arr[block1 + 2, block2 + 1] + arr[block1 + 2, block2 + 2])
    # arr[block1 + 2, block2 + 2] = arr[block1 + 1, block2 + 1] - 2 * (u12 + u21)

    # m, n = size(arr)

    # levels = Int(intlen(m ÷ n))
    # _bit_reverse_permutation_2d!(arr, levels, 2)
    # for i in 1:n:m
    #     @views _fft2_pow2_cooley_square!(arr[i:i+n-1, :])
    # end

    # len = 2 * n
    # while len <= m
    #     half = len >> 1
    #     roots = transpose(cispi.(-2 * (0:half-1) / len))
    #     for block in 0:len:m-1
    #         for i in 1:half
    #             @. @views arr[block + i + half, :] *= roots[i]
    #             @. @views arr[block + i, :] += arr[block + i + half, :]
    #             @. @views arr[block + i + half, :] = arr[block + i, :] - 2 * arr[block + i + half, :]
    #         end
    #     end
    #     len <<= 1
    # end
