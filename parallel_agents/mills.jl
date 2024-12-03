using LinearAlgebra
using CUDA
using StaticArrays 
using GLMakie 
using OrdinaryDiffEq
import AcceleratedKernels as AK 
using OrdinaryDiffEq, TerminalLoggers

const use_gpu = true 
const dtype = use_gpu ? Float32 : Float64  
const SVec2 = SVector{2, dtype}
const SVec4 = SVector{4, dtype}

to_array(X) = reinterpret(reshape, dtype, X)
to_svecs(X) = reinterpret(reshape, SVec4, X)

@inline dist²(dx) = sum(z -> z^2, dx)
dist(x) = sqrt(dist(dx))

@inline wrap(dx, L, L_inv = inv(L)) = @. dx - L * round(dx * L_inv)
@inline dist²(dx, L, L_inv = inv(L)) = dist²(wrap(dx, L, L_inv))
dist(x, L, L_inv = inv(L)) = sqrt(dist(dx, L, L_inv))




N = 10_000 

mass = 1000.0
L = 10.0
Ca = 0.5
la = 2.0
Cr = 1.0
lr = 0.5
alpha = 1.6
beta = 0.5
v0 = sqrt(alpha/beta)

p = (;N , mass, L, Ca, la, Cr, lr, alpha, beta, v0) 
p = map(dtype, p)

pos = p.L * rand(SVec2, N)
vel = randn(SVec2, N)
@. vel = p.v0 * vel / norm.(vel)

dt = dtype(0.01) 

z0 = cu( [to_array(pos); to_array(vel)] )

function Morse_force(x, a, C, l, r)
    r_by_l = (r/l)^a 
    return (r_by_l / r^2 * a * C * exp(-r_by_l)) * x
 end
 
function mills!(dz, z, p, t)

    XV = to_svecs(z)
    VF = to_svecs(dz)

    pos_ind = SA[1,2]
    vel_ind = SA[3,4]
    L² = p.L^2
    L = p.L
    L_inv = inv(p.L)

    AK.foreachindex(XV, block_size = 256) do i
        XVi = XV[i]
        Xi = XVi[pos_ind]
        Vi = XVi[vel_ind]

        Fi = zero(Xi)

        for j in eachindex(XV) 
            Xj = XV[j][pos_ind]
            x = wrap(Xj .- Xi, L, L_inv) 
            r² = dist²(x)

            if zero(r²) < r² < L²  
                r = sqrt(r²)
                Fi += Morse_force(x, 1, p.Ca, p.la, r)
                Fi -= Morse_force(x, 1, p.Cr, p.lr, r)
            end
        end

        Fi *= p.mass / p.N 
        Fi += (p.alpha - p.beta * dot(Vi, Vi)) .* Vi

        VF[i] = [Vi; Fi] 
    end
    AK.synchronize(AK.get_backend(z))
end

dz = similar(z0)
@time mills!(dz, z0, p, nothing)

# z0_cpu = [to_array(pos); to_array(vel)]
# dz_cpu = similar(z0_cpu)
# @time mills!(dz_cpu, z0_cpu, p, nothing)

tspan = dtype.( (0.0, 50.0) )
prob = ODEProblem(mills!, z0, tspan, p)
@time sol = solve(prob, RK4(), dt = dt, adaptive = false, saveat = dtype[0,1,2,3,4,5,10,40,70,100], progress = true)
# 145.659184 seconds (9.15 M allocations: 379.076 MiB, 0.09% gc time, 0.78% compilation time: 93% of which was recompilation)


begin 
    fig = Figure(size = (1000, 400)) 

    offset = (1.0, 5.0) .- p.L / 2

    for (k, t) in enumerate(sol.t)
        i, j = divrem(k -1, 4) .+ 1

        ax = Axis(fig[i,j], aspect = DataAspect())
        xv = Array(sol.u[k])
        x = @. mod(xv[1,:] - offset[1], p.L)
        y = @. mod(xv[2,:] - offset[2], p.L)
        angle = atan.(xv[4,:], xv[3,:])
        scatter!(ax, x, y, rotation = angle, marker = '✈', markersize = 0.05, color = :black, markerspace = :data)

        if k > 1
            linkaxes!(content(fig[1,1]), ax)
        end
    end

    fig 
end