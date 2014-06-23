MISH: Multiple Implementation Simple Hydrocode
======================

Description
-----------

MISH is a set of Miniapp implementation of a simple Godunov hydrocode using a number of standard HPC tools.

The codes are designed primarily to test the relative performance of the selected tools when running a standard simple compressible fluid dynamics model.

Each of the included subdirectories contains a single implementation of a Godunov Hydrocode using a some set of HPC tools such as MPI or OpenMP.

Usage
------

Command Line arguments
<dl>
  <dt>init</dt>
  <dd>Initial condition, one of:
    <dl>
    <dt>sod</dt>
    <dd>Sod Problem</dd>
    <dt>crn</dt>
    <dd>Quadrant Shock</dd>
    <dt>wsc *size*</dt>
    <dd>Weak Scaling, Sod Problem, # cells in active dim scaled by *size*</dd>
    <dt>scs *size*</dt>
    <dd>Adjustable Quadrant Shock, both dimensions scaled by *size*</dd>
    </dl>
  </dd>
</dl>

Output
----

In order to ease comparison and analysis using the timing results of the implementations, the data is printed using the following layout. Two of the entries are inserted with the intent that they be modified to record some additional data wwhich is not trivially available to the engine code such as a label for the precise type of machine run on or the initializaiton used.

````
TIME:cType,mType,init,nproc,nth,niters,ncells,runt
````

<dl>
<dt>cType</dt>
<dd>Type of implementation (C, OMP, MPI,...)</dd>
<dt>mType</dt>
<dd>"CPU:?" or "GPU:?" depending  on implementation</dd>
<dt>init</dt>
<dd>Always "Init"</dd>
<dt>nproc</dt>
</dd>Number of Processes used for computation</dd>
<dt>nth</dt>
<dd>Number of threads used for computation</dd>
<dt>niters</dt>
<dd>Number of timesteps (iterations of main loop) run</dd>
<dt>ncells</dt>
<dd>Number of cells in computation</dd>
<dt>wRunt</dt>
<dd>Wallclock runtime in seconds</dd>
</dl>
