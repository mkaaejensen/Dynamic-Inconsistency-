General notes: 

- For main demonstration run optimal_growth.m
- Default CMA-ES implementation is the recommended algorithm in growth/social planner models (except at very high precision where tsa may be preferable)
- A CMA-ES implementation which always generates some candidates with large step size will likely work very well in this class of problems (not implemented/included). This is suggested by the restart procedure described below as that effectively restarts the step size.


1) optimal_growth.m and optimal_growth10.m

Produce figure 4, panel A and B. Also plots the known analytical solution. 

Note that optimal_growth10.m initiates from a random starting point (hence will not reproduce exactly the same graph as in the paper except by coincidence).

2) optimal_growth0_90.m 

Parallel to 1) but with twice as many grid points (gridprec=40, meaning 40 grid points per unit interval). Recommended: Try reducing gridprec to for example 20 (the point here is that as present bias grows stronger, more precise computations become necessary)

3) optimal_growth_highprecision.m

Like 1) and 2) but with 500 grid points per unit interval (gridprec=500). Uses tsa. Used to compute figure in the paper (visually not distinguishable from 2)).

3) optimal_growthsigma2.m

By default loadoutput=1 so that precomputed figure is illustrated. This shows an optimal growth Markov equilibrium with jumps (a "Krussel and Smith 2003-type solution"). What is interesting is that the loss is at machine precision level meaning that this, without a doubt, well-approximates an exact equilibrium. 

To generate figure, the fastest way is to "restart": Change loadoutput to 0. Run algorithm and reduce loss to low e-6 region. Cancel, change loadoutput to 1 and run again. The resulting loss is at machine precision level (of course, even e-6 is an extremely low loss) 