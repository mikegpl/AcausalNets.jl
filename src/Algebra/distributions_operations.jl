#=
distribution_operations:
- Julia version: 1.0
- Author: marcin
- Date: 2018-08-18
=#

using LinearAlgebra

# as defined in
# https://arxiv.org/pdf/0708.1337.pdf
# (equation 29)
function star_n(n::Float64)
    a_pow = 1 / (2 * n)
    b_pow = 1 / n
    (A::AbstractMatrix, B::AbstractMatrix) -> ((A ^ a_pow) * (B ^ b_pow) * (A ^ a_pow)) ^ n
end


const star = star_n(1.0)

# as defined in
# https://arxiv.org/pdf/0708.1337.pdf
# (equation 32)
unstar(C, B)::AbstractMatrix = star(pinv(B), C) # odwrotność star (za #6)

event(system::AbstractMatrix, e::AbstractMatrix) = (e * system * e) / tr(e * system)

