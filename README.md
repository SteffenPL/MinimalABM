# MinimalABM

Two minimal examples for agent-based models with pure Julia. 
One can consider using [Agents.jl]([ts.jl/stable/](https://juliadynamics.github.io/Agents.jl/stable/)) instead!

## fibers
The model `fibers.jl` generates a simulation with linesegments that create and destroy adhesive bonds (linear springs) between them.
If torn apart, this might model the behavior of biological fibers under stress:



https://github.com/SteffenPL/MinimalABM/assets/4634605/97e5a07c-4792-4fef-9970-4b3ce13ea831



## minimal
The model in `miminal.jl` just shows how to simulate a bunch of spheres with non-overlap and some attraction between them. This code 
shows a minimal structure typical for agent-based models.

https://github.com/SteffenPL/MinimalABM/assets/4634605/6aa18e9e-8d44-451b-90fb-9b0f25cf5d9c

