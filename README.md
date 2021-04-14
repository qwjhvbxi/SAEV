# SAEV

Shared Autonomous Electric Vehicles simulation and optimization model.


Setup
-----

To setup the model, first run:

`setup`

This will prompt you to choose a folder where to save output files and parameters.  The model tends to create large savefiles. For this reason, it is better to have the folder with savefiles in a location with enough space separated from the code.


Using the model
---------------

How to launch a simulation:

`Result=main(P,extsave,info)`

The first input variable `P` is a struct-type with the parameters of the simulation, the second variable `extsave` indicates whether to save/load the results, the third variable `info` indicates what level of information to show during a simulation run. See below for more information.
A typical call would be:

`Result=main(P,1,2)`

`P` can be constructed from scratch, but default values of `P` are available with the function `pdefault`. For example, this function creates a `P` for simulations with the `NYC2016` dataset:

`addpath functions`

`P=pdefault('NYC2016')`

The first line add the folder containing pdefault to the current path. `P` includes default values for fleet size, technological properties of vehicles, etc.

Previously saved default values of `P` can be called with the function `getp`, which reads the content of the subfolder `par` in the data folder.

## Parameter variable `P`

The parameter variable `P` need to have the following fields:

| Field         | format| description |
| ------------- | ------------- | ------------- |
| `scenario`    | char  | Name of scenario |
| `tripfolder`  | char  | Name of trip folder |
| `tripday`     | int   | Number of day | 
| `gridfile`    | char  | Name of file with electricity data |
| `gridday`     | int   | Number of day | 
| `m`           | int   | Fleet size | 
| `modechoice`  | logical   | Wether to use modechoice | 
| `carbonprice` | double   | Carbon price (per MWh) | 
| `Sim`         | struct | Simulation settings (see below) |
| `Tech`        | struct | Technical parameters of vehicles (see below) |
| `Operations`  | struct | Operational settings (see below) |
| `Charging`    | struct | Charging settings (see below) |
| `Relocation`  | struct | Relocation settings (see below) |

Additional optional fields:

`Pricing`


extsave and info variables
--------------------------

`extsave` can be:

* -1: run the simulation. Do not save results
* 0: try to load simulation if exists, otherwise run it. Do not save results
* 1: try to load simulation if exists, otherwise run it. Save results
* 2: run the simulation. Save results (may overwrite existing results)

Generally `extsave=1` is used (default value).

`info` can be:
* 0: do not show progress information -- good for very short simulations
* 1: show simplified progress information (progress bar) -- good for short simulations
* 2: show detailed progress information (estimated time, etc.) -- good for long simulations
* negative numbers: periodically print a line with progress summary -- useful with parallel computing 

Generally `info=2` is used (default value).



