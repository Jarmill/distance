# distance
Estimate the distance of closest approach between points along trajectories of x'=f(t, x) starting from an initial set X0 and an unsafe set Xu over a time horizon of \[0, T\]. If possible, recover the trajectories that attain the closest approach.

An occupation-measure framework is used to find a convergeng sequence of upper bounds in rising degree to the true minimal distance value. The approximate-optimal trajectories may be recovered if moment matrices satisfy rank conditions (up to numerical accuracy).

## Dependencies

- Gloptipoly3: http://homepages.laas.fr/henrion/software/gloptipoly/. A customized version is included in as a submodule to this repository (gloptipoly_mod)
- YALMIP: https://yalmip.github.io/
- Mosek: https://www.mosek.com/ (or any solver compatible with YALMIP)

All code is written and tested on Matlab R2021a.


## Instructions
The `distance_manager` class has two input arguments: `lsupp` and `f`. The object `lsupp` of type `unsafe_support` is the options structure that defines properties of the system, such as its variables, the distance function, and the support set of (X, X0, Xu). The input `f` are the dynamics of the system.

Execute `distance_manager.run(order)` to formulate and solve the distance estimation program. Call `distance_manager.loc.recover()` to attempt rank-based recovery of near-optimal trajectories. `The `order` is the relaxation order involving moments of degree `2*order`. 

Examine and run `experiments/flow_distance_test.m` as an example. The README in the `experiments` folder describes all test scripts.

## Reference
https://arxiv.org/abs/2110.14047 "Bounding the Distance to Unsafe Sets with Convex Optimization"

https://homepages.laas.fr/vmagron/brainpop/DistanceJared.pdf Presentation given at BrainPOP seminar at LAAS-CNRS on June 28, 2021

## Contact
For comments and questions please email [Jared Miller](mailto:miller.jare@northeastern.edu?Subject=distance).

