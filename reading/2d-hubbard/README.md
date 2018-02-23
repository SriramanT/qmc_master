# Description

This is a fairly accessible article that shows how to convert the ground state projection operator with the Hubbard Hamiltonian into a huge sum which then allows you to do Monte Carlo sampling. It also touches on how to update the Green’s function efficiently; the checkerboard breakup; numerical stabilization of matrix multiplications in the algorithm; and how physical quantities can be calculated from the Green’s functions.

Note that this is written only for finite-temperature AFQMC with Metropolis-like sampling. It took me a while to understand the subtle difference between doing AFQMC with Metropolis-like and with branching random walk.