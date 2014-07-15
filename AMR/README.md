# Interface for Adaptive Mesh Refinement with Concurrent Task-Based Models
---
We implemented an abstract interface for adaptive mesh refinement (AMR) with
visualizations.  We translated that interface into programming languages with
built-in concurrency primitives in order to test modern computing techniques on
scientific computation codes.


## Quad Tree
---
In this directory we have the code for an AMR interface implemented using a
pointer-based quad tree data structure written in C++.

## Concurrent Models
---
In this directory we have the code for two concurrent task-based AMR models: one
written in Go and one written in D.  We also include some of the basic
benchmarking tests to demonstrate speed up and the overhead cost of concurrency.

