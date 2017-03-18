
## HIGH PERFORMANCE COMPUTING (OpenMp, MPI and GPU).

---------------------------------------
###  Gauss-Seidel iterative method 
----------------------------------------
This project is about Parallelizing Gauss-Seidel Iterative method. To do this we make use of the inustrial standard libraries <i> OpenMP, MPI and GPU</i> in C++  and  then do some performance analysis.

--------
## Background: Direct Method Or Gaussian Elimination
--------
As a numerical technique the [Gaussian elimination technique](https://en.wikipedia.org/wiki/Gaussian_elimination), is rather unusual because it is <i> direct </i>. That is a result is obtained after a single application of Gaussian elimination. Once a "solution" has been obtained, Gaussian elimination offers no method of refinement. The lack of refinement can be a problem since Gaussian elimination is sensitive to rounding error.

-----
## Numerical Technique
----

Numerical techniques more commonly involve an iterative method. In this work we look at one such iterative method for approximating solution of system of linear equations : <b> The Gauss Siedel Method</b>. We will focus on how we can parallelize  this numerical method and report on its performance.
