using GLMakie, StaticArrays
Makie.inline!(true)

const Dim = 2  # dimension (keep it as a global constant)

p = (
    N = 100, 
    R = 0.1, # repulsion radius
    γ = 0.2, # attraction to center
    μ = 1.0, # damping
    σ = 0.0, # noise
    #
    dt = 0.01,
    tend = 1.0
)

# init cell centers 
function initstate(p)
    (;N) = p

    X = rand(Dim, N)  # center of cells 
    V = zeros(Dim, N) # velocity of cells

    return (; X, V, t = Ref(0.0))
end

# simulate the model of 'N' overdamped spheres with:
# - non-overlap between spheres with radius R
# - attraction to center
# e.g. ̇Xᵢ = -∇U(Xᵢ) + γ(Xᵢ - X₀) + √(2μ)ξᵢ   where U is a non-smooth term for non-overlap

function simulate(s, p)
    (; N, R, γ, μ, σ, dt, tend) = p

    (;X, V, t) = s

    n_steps = ceil(Int, tend / dt)
    
    sol = [deepcopy(s)]

    for step in 1:n_steps

        # attraction to center
        V .= -γ .* (X .- 0.5)
        X .+= dt/μ * V
        X .+= sqrt(2σ*dt) .* randn(Dim, N)

        # non-overlap
        for i in 1:N
            for j in 1:i-1
                r = X[:,i] - X[:,j]
                d = norm(r)  # this is a bit slow, but let's keep it simple

                if d < 2*p.R 
                    @. X[:,i] += 0.5 * (2*p.R - d) * r / d
                    @. X[:,j] -= 0.5 * (2*p.R - d) * r / d
                end 
            end
        end

        t[] += dt

        push!(sol, deepcopy(s))
    end
    return sol
end


s = initstate(p)
sol = simulate(s, p)

# plot inital and terminal state
circle_properties = (marker = Makie.Circle, markersize = 2*p.R, markerspace = :data, color = :green)

fig = Figure()

ax1 = Axis(fig[1, 1], aspect = DataAspect())
scatter!(ax1, sol[2].X[1,:], sol[2].X[2,:]; circle_properties...)

ax2 = Axis(fig[1, 2], aspect = DataAspect())
scatter!(ax2, sol[end].X[1,:], sol[end].X[2,:]; circle_properties...)

linkaxes!(ax1, ax2)
fig

# make movie: 
state = Observable(sol[1])
X = @lift $state.X[1,:]
Y = @lift $state.X[2,:]

fig = Figure(resolution = (800, 800))
ax = Axis(fig[1, 1], aspect = DataAspect())
scatter!(ax, X, Y; circle_properties...)
limits!(ax, -1.0, 2.0, -1.0, 2.0)


record(fig, "minimal.mp4", 1:length(sol)) do i
    state[] = sol[i]
end