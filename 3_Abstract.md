

-------
## Abstract
-------

<p> In this work, we present the  parallel implementation of the Gauss-Seidel (GS) iterative algorithm for the solution of linear systems of equations  using MPI . 
<p> To approach this problem, we discuss and implement the distributed parallel Gauss-Seidel algorithm which exploits the sparsity of the system to extract parallelism. In this approach, the coefficient matrix and the right-hand side of the linear algebraic system are first divided into row-blocks in the natural rowwise-order. Thereafter the the row-blocks are distributed among local memories of all processors according to the torus-wrap mapping techniques. The solution iteration vector is cyclically conveyed among processors at each iteration so as to decrease the communication. The result of such approach is true parallel model with convergence rate as good as the sequential execution. This first approach  basically takes advantage  of the sparsity of the system to extract parallelism. </p>
<p>In addition to this, we   provide  parallelization scheme for fully-coupled systems (the dense Gauss-Siedel approach) which takes advantage of MPI processing paradigms. 
</p>

<p> 
Both strategies are evaluated through performance measurements on  high-performance computing architectures.
</p>
</p>
