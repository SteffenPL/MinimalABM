using GLMakie, Graphs, BoundedDegreeGraphs, LinearAlgebra, StaticArrays
Makie.inline!(true)

# Note on Julia: The following Julia code uses a lot of 'NamedTuples'!
# This is a very convenient way to pass around parameters and state
# variables. For more information, see the Julia documentation:
# https://docs.julialang.org/en/v1/base/base/#Core.NamedTuple
#
# In short:
# - NamedTuples are created with nt = (; a=1, b=2)
# - NamedTuples are accessed with nt.a or nt.b
# - NamedTuples can be unpacked with (; a, b) = nt
# - NamedTuples can be merged with (; nt1..., nt2...)


# the following is a simulation with 'N' fibers which dynamically
# create and destroy connections between each other.
const Dim = 2  # dimension of the system

# Parameters   (in units μm, minutes, pN)
p = (
    N = 20, # number of fibers
    fiberlength = 1.0,  # unit: μm,
    μ = 1.0, # damping coefficient, unit: pN/μm
    I = 10.0, # moment of intertia
    κ = 1.0, # spring constant, unit: pN/μm
    k1 = 0.1, # rate of connection creation, unit: 1/min
    k0 = 0.1, # rate of connection destruction, unit: 1/min
    #
    f_ext = (xi, t) -> xi[1] < 0.5 ? (-1.0,0.0) : (1.0,0.0), # external force
    σ = 1e-4, # noise 
    # 
    dt = 0.01, 
    tend = 5.0,
)

# initial conditions
function initstate(p)
    (;N) = p  # unpack variables

    X = rand(Dim, N)  # center of fiber 
    P = randn(Dim, N)  # direction of fiber (normalized)
    
    for i in 1:N
        P[:,i] ./= norm(P[:,i])
    end

    maxdegree = 10
    bonds = BoundedDegreeMetaGraph(N, maxdegree, (;l1 = 0.0, l2 = 0.0)) 
    # graph of connections between fibers where each edge stores the information
    # about the points of attachment

    return (; X, P, bonds, t = Ref(0.0))
end

# this function is just norm(X[:,i] - X[:,j]), but it makes the code faster.
function dist(X, i, j) 
    d = 0.0 
    for k in 1:Dim
        d += (X[k,i] - X[k,j])^2
    end 
    return sqrt(d) 
end 

function orthcomplement(x)
    return (-x[2], x[1])
end

function cross2D(x,y)
    return x[1]*y[2] - x[2]*y[1]
end

function simulate(s, p) 

    s = deepcopy(s) # make a copy of the state
    (;N) = p  # unpack variables
    (;X, P, bonds) = s # unpack state

    sol = [deepcopy(s)]  # store the state at each time step
    obs = nothing 

    F = zeros(Dim, N)  # force on each fiber
    ω = zeros(Dim, N)  # torque on each fiber

    while s.t[] < p.tend 
        
        # update time
        s.t[] += p.dt
        F .= 0.0
        ω .= 0.0

        for j in 1:N
            for i in 1:j-1            
                d = dist(X, i, j)
                if d < p.fiberlength
                    if !has_edge(bonds, i, j)
                        Pi = P[:,i]
                        Pj = P[:,j]
                        r = X[:,j] - X[:,i]

                        # compute the intersection point
                        a1 = dot(r, Pi)
                        nij = a1*Pi - r 
                        l2 = norm(nij)^2 / dot(nij, Pj) 
                        l1 = dot(r + Pj*l2, Pi)

                        @assert X[:,j] + l2*Pj ≈ X[:,i] + l1*Pi

                        if abs(l1) < p.fiberlength / 2 && abs(l2) < p.fiberlength / 2
                            # the fibers are connected
                            add_edge!(bonds, i, j, (;l1, l2))  
                            # important: i < j since unordered graph; otherwise order l1/l2 is swapped.
                        end
                    end
                end
            end
        end

        for e in edges(bonds)
            i = e.src
            j = e.dst
            (;l1, l2) = bonds[e]

            # compute the force on the fiber
            r = X[:,j] - X[:,i]  # vector from Xi to Xj
            Pi = P[:,i]
            Pj = P[:,j]
            nij = r + Pj*l2 - Pi*l1  # vector from Xi + l1*Pi to Xj + l2*Pj
            F[:,i] .+= p.κ * nij
            F[:,j] .-= p.κ * nij

            if Dim == 2
                ω[:,i] .+= p.κ * cross2D(Pi*l1, nij) * orthcomplement(Pi)
                ω[:,j] .-= p.κ * cross2D(Pj*l2, nij) * orthcomplement(Pj)
            else 
                @error "not implemented"
            end
        end

        # update position
        for i in axes(X, 2)
            Xi = X[:,i]
            Xi .+= sqrt(p.dt * p.σ) * randn(Dim)
            Xi .+= p.dt / p.μ * (F[:,i] .+ p.f_ext(Xi, s.t[])) 
            X[:,i] .= Xi
        end

        # update rotation
        for i in axes(P, 2)
            Pi = P[:,i]
            Pi .+= p.dt / p.μ / p.I * ω[:,i] 
            P[:,i] .= Pi ./ norm(Pi)
        end

        # update solution 
        push!(sol, deepcopy(s))
    end

    return (;sol, obs) 
end

# This function is used to make a movie in the end.
# Essentially, we extract only the data needed for plotting using 
# the 'lift' function. This is a very powerful feature of Makie, as we can make 
# movies super fast.
# See: https://docs.makie.org/stable/explanations/nodes/ 
function plotdata(state, p) 
    
    fibers = lift(state) do s 
        [ 
            (   Point2f(s.X[:,i] - 0.5 * p.fiberlength * s.P[:,i]), 
                Point2f(s.X[:,i] + 0.5 * p.fiberlength * s.P[:,i])     ) 
            for i in axes(s.X, 2)
        ]
    end

    adh_bonds = lift(state) do s 
        [ 
            (   Point2f(s.X[:,e.src] + s.bonds[e.src,e.dst].l1 * s.P[:,e.src]), 
                Point2f(s.X[:,e.dst] + s.bonds[e.src,e.dst].l2 * s.P[:,e.dst]) ) 
            for e in edges(s.bonds)
        ]
    end

    return (;fibers, adh_bonds)
end



function plotstate!(ax::Axis, data)
    (;fibers, adh_bonds) = data
    linesegments!(ax, fibers, color=(:black, 0.5), linewidth = 2)
    linesegments!(ax, adh_bonds, color=(:red,0.5), linewidth = 2)    
end

function plotstate(s, p) 
    fig = Figure()
    ax = Axis(fig[1,1])
    plotstate!(ax, plotdata(s, p))
    fig
end



#####################
### Simulation ######
#####################


s = initstate(p)
#add_edge!(s.bonds, 1, 2)
#add_edge!(s.bonds, 2, 3)


sol, obs = simulate(s, p)

s_node = Observable(s)
plotstate(s_node, p)
s_node[] = sol[2]
plotstate(s_node, p)


# make a video 

fig = Figure()
ax = Axis(fig[1,1], title="ABM for fiber dynamics")
data = plotdata(s_node, p)
plotstate!(ax, data) 

record(fig, "demo.mp4", eachindex(sol)) do i 
    s_node[] = sol[i]  #update plots
end 

fig
