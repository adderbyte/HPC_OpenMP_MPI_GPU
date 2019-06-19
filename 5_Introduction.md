------
## Introduction
------

<p> The solving of linear algebraic systems:
<p align = "center"> Ax = b</p>
<p>lies at the core of many scientific and engineering simulations.
Many practical problems can be translated into a large scale linear algebraic system. 
Therefore the parallel solving of large scale linear algebraic systems is of great importance in 
the scientific computation field. Methods for a linear algebraic system generally fall into two 
categories: <i> direct methods and iterative methods</i>.</p>

<p> A direct method is one in which  a fixed number of operations are carried out once, 
at the end of which the solution is produced. Gauss elimination and related strategies on a linear 
system is an example of such methods. From this, we know that the direct methods can arrive at the solution within
a limited number of steps. However, direct methods may be impractical if the coefficient matrix of the linear 
algebraic system to be solved is large and sparse because the sought-after factor can be dense. In addition to this,
Direct methods are often too expensive in either computation time or computer memory requirements, or possible both.</p>

<p>This is the 
 why the iterative methods, which are able to take advantage of sparse systems, are often preferable compared
to  direct methods in engineering fields. Through iterative methods, an approximate solution within the error 
tolerance could be obtained under the assumption that the algorithm is convergent. Among classical iterative methods, 
the Gauss–Seidel method has several interesting properties. In general, if the Jacobi method also converges, 
the Gauss–Seidel method will converge faster than the Jacobi method. Even though the SOR method with the optimal
relaxation parameter is faster than the Gauss–Seidel method, however, choosing an optimal SOR relaxation parameter 
can be difficult for many problems of practical interest [2]. Therefore, the Gauss–Seidel method is very attractive 
in practice and it is usually used as the smoother of the multigrid method for partial differential equations which 
typically yields good multigrid convergence properties.</p>

</p>
</p>
