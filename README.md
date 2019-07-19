# HomologyLive
Source code of A. Rahm's software "HomologyLive"


What the software HomologyLive does: Statistical Physics simulations involve a large number of particles, whose geometrical distance to each other is calculated precisely and efficiently. This means that we have a natural distance function on the set of particles, and hence on its Vietoris-Rips complex, we can interpret naturally the 1- and 2-dimensional homology generators as loops and bubbles in the particle configuration. With HomologyLive, Alexander D. Rahm has made an implementation of the Vietoris-Rips complex adapted to this setting: Truncated at dimension 3, to stay within the physical dimensions of the ambient space, and with boundaries recorded as sparse matrices, in order to allow for huge numbers of particles in the simulations while keeping the linear algebra feasible on the machine. While a dense matrix computation with the full Vietoris-Rips complex - as implemented in SAGE - takes a big multiple of the machine resources of the simulations, HomologyLive in fact only takes a little fraction of them. The computation of the ranks of the sparse matrices is done using Linbox (see https://github.com/linbox-team ).


There are currently three branches, each one specialized to a specific application in natural sciences:

* Statistical Physics simulations

* Protein-Protein Interaction Networks in Genomics

* Gene occurrency in cells

Publications related to these applications can (respectively, will) be found on Alexander D. Rahm's homepage:

http://math.uni.lu/~rahm/
