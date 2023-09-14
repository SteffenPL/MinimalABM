# MinimalABM

Two minimal examples for agent-based models with pure Julia. 
One can consider using [Agents.jl]([ts.jl/stable/](https://juliadynamics.github.io/Agents.jl/stable/)) instead!

## fibers
The model `fibers.jl` generates a simulation with linesegments that create and destroy adhesive bonds (linear springs) between them.
If torn apart, this might model the behavior of biological fibers under stress:



https://github.com/SteffenPL/MinimalABM/assets/4634605/808ee79f-33e7-4aeb-84b0-e1a28f54ca83



There are several bits which are missing for a relatistic biological model:
- the length of linesegments is too long (should be in the order relative to the typical fiber curvature)
- link detachment rate to current extension of the link
- fibers do not repulse each other
- there is no alignment between fibers (this is maybe a modelling question one needs to answer); this leads to linesegments rotating too fast
- fibers cannot have multiple links right now. This could improve the model.
- units are not on a realistic scale.

The implementation is still a bit simplified, for better performance one should:
- use static arrays to store positions
- use a better method to find fibers which are closeby (cell list method, spatial hash trees)

From a math perspective:
- taking the limit $y -> \infty$ might reduce to model to 
simple a bunch of linear springs which die with a certain rate. That
could lead to an analytical expression for the breaking point, simply in terms of the rates and number of springs.

## minimal
The model in `miminal.jl` just shows how to simulate a bunch of spheres with non-overlap and some attraction between them. This code 
shows a minimal structure typical for agent-based models.

https://github.com/SteffenPL/MinimalABM/assets/4634605/6aa18e9e-8d44-451b-90fb-9b0f25cf5d9c

