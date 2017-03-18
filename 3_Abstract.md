

-------
## Abstract
-------

The Gauss-Seidel method is very efficient for solving problems such as tightly-coupled constraints with possible redundancies. 
However, the underlying algorithm is inherently sequential. A distributed memory parallel Gauss–Seidel 
algorithm for linear algebraic
systems is presented, in which a parameter is 
introduced to adapt the algorithm to different distributed memory parallel architectures. 
In this algorithm, the coefficient matrix and the right-hand side of the linear algebraic 
system are first divided into row-blocks in the natural rowwise-order according to the performance of 
the parallel architecture in use. And then these row-blocks are distributed among 
local memories of all processors through torus-wrap mapping techniques. The solution iteration 
vector is cyclically conveyed among processors at each iteration so as to decrease the communication.
The algorithm is a true Gauss–Seidel algorithm which maintains the convergence rate of the serial 
Gauss–Seidel algorithm and allows existing sequential codes to run in a parallel environment with a
little investment in recoding.

In this paper, we will also study several parallelization schemes for fully-coupled systems, 
unable to be parallelized by existing methods, taking advantage of recent 
many-cores architectures offering fast synchronization primitives.

