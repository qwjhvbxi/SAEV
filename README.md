# SAEV
 
Shared Autonomous Electric Vehicles simulation and optimization model.


Setup
-----

The model tends to create large savefiles. For this reason, it is better to move the folder with savefiles in a location with enough space separated from the code.
The folder with savefiles is called `extdata`.
To change the location of `extdata`, simply move it somewhere else and update its positioning in the file `functions/setDataFolder.m`.
For example, change from:

`DataFolder='extdata/';`

to:

`DataFolder='C:\Users\example\Documents\extdata\';`


Using the model
---------------

How to launch a simulation:

`Result=generalC(P,extsave,info)`

The first input variable `P` is a struct-type with the parameters of the simulation, the second variable `extsave` indicates wheter to save/load the results, the third variable `info` indicates what level of information to show during a simulation run. See below for more information.
A typical call would be:

`Result=generalC(P,1,2)`

`P` can be constructed from scratch, but default values of `P` are available with the function `cpar`. For example, this function creates a `P` for simulations with the `NYC2016` dataset:

`P=cpar('NYC2016')`

This includes default values for fleet size, technological properties of vehicles, etc.

extsave and info variables
--------------------------

`extsave` can be:
* -1: run the simulation. Do not save results
* 0: try to load simulation if exists, otherwise run it. Do not save results
* 1: try to load simulation if exists, otherwise run it. Save results
* 2: run the simulation. Save results (may overwrite existing results)

Generally `extsave=1` is used.

`info` can be:
* 0: do not show progress information -- good for very short simulations
* 1: show simplified progress information (progress bar) -- good for short simulations
* 2: show detailed progress information (estimated time, etc.) -- good for long simulations
* negative numbers: periodically print a line with progress summary -- useful with parallel computing 

Generally `info=2` is used.