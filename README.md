# SAEV

Shared Autonomous Electric Vehicles simulation and optimization model.

Code by Riccardo Iacobucci

https://riccardoiacobucci.com

## Cite

If you use this code in your research, please cite:

1. R.Iacobucci, R.Bruno, J.D.Schmöcker. *Integrated optimisation-simulation framework for scalable smart charging and relocation of shared autonomous electric vehicles*. Energies, 2021. doi: [10.3390/en14123633](https://doi.org/10.3390/en14123633)

If you use the FCR module, please additionally cite:

1. R.Iacobucci, J. Donhauser, J.D.Schmöcker, M. Pruckner. *Frequency Control Reserve Provision from a Fleet of Shared Autonomous Electric Vehicles*. IEEE MT-ITS 2021

If you use the pricing module, please additionally cite:

1. R.Iacobucci, J.D.Schmöcker. *Dynamic pricing for ride-hailing services considering relocation and mode choice*. IEEE MT-ITS 2021

## Quick start

To test the model, you can run the example data. 
1. Run `main`. on the command line. The function will list the available simulation scenarios. 
1. When requested, input the scenario name `NYC2018` or just `1`. 

The example data is based on the New York yellow taxicab data for 2018, available at https://www1.nyc.gov/site/tlc/about/tlc-trip-record-data.page .

## Using the model

How to launch a simulation:

`Result=main(P,extsave,info)`

The first input variable `P` is a struct-type with the parameters of the simulation. 
The second variable `extsave` indicates whether to save/load the results (to use this functionality first include `DataHash` in the path. You can find it at https://uk.mathworks.com/matlabcentral/fileexchange/31272-datahash).
The third variable `info` indicates what level of information to show during a simulation run. See below for more information.
A typical call would be:

`Result=main(P,1,2)`

`P` can be constructed from scratch, but default values of `P` are available with the function `pdefault`. For example, this function creates a `P` for simulations with the `NYC2016` dataset:

`addpath functions`

`P=pdefault('NYC2016')`

The first line add the folder containing pdefault to the current path. `P` includes default values for fleet size, technological properties of vehicles, etc.

Previously saved default values of `P` can be called with the function `getp`, which reads the content of the subfolder `par` in the data folder.

### Parameter variable `P`

The parameter variable `P` need to have the following fields:

| Field         | format| description |
| ------------- | ------------- | ------------- |
| `scenario`    | char  | Name of scenario |
| `tripfolder`  | char  | Name of trip folder |
| `tripday`     | int   | Number of day | 
| `gridfile`    | char  | Name of file with electricity data |
| `gridday`     | int   | Number of day | 
| `m`           | int   | Fleet size | 
| `modechoice`  | logical| Wether to use modechoice | 
| `carbonprice` | double | Carbon price (per MWh) | 
| `Sim`         | struct | Simulation settings (see below) |
| `Tech`        | struct | Technical parameters of vehicles (see below) |
| `Operations`  | struct | Operational settings (see below) |
| `Charging`    | struct | Charging settings (see below) |
| `Relocation`  | struct | Relocation settings (see below) |

Additional optional fields:

| Field         | format | description |
|`Pricing` 		| struct | Pricing info (see below) |
|`FCR` 			| struct | Settings for frequency control reserve (see below) |

#### `Sim`

| Field         | format| description |
| ------------- | ------------- | ------------- |
|e              | int       | Length of time step (minutes) |
|mpcpredict     | logical   | Perfect prediction? |

#### `Tech`

| Field         | format| description |
| ------------- | ------------- | ------------- |
|battery    | double | Vehicle battery capacity (kWh) |
|chargekw   | double | Vehicle max. charging power (kW) |
|consumption| double | Vehicle consumption (kWh/min) |
|cyclingcost| double | Battery cycling cost (per cycle) |
|efficiency | double | Battery round-trip efficiency [0,1] |

#### `Operations`

| Field         | format| description |
| ------------- | ------------- | ------------- |
|initialsoc| double| Initial state of charge [0,1] |
|minsoc| double| Minimum state of charge [0,1] |
|maxsoc| double| Maximum state of charge [0,1] |
|v2g| logical| Enable V2G discharge? |
|v2gminsoc| double| Minimum state of charge for V2G [0,1] |
|maxwait| double|  (minutes) |
|maxidle| double| (minutes) |

#### `Charging`

| Field         | format| description |
| ------------- | ------------- | ------------- |
| mthor     | int   | Time horizon for charging optimization (minutes) |
| extrasoc  | double| Extra charge level over min. SoC for charging optimization [0,1] |
| beta      | int   | Lenght of charging optimization interval (minutes) |

#### `Relocation`

| Field         | format| description |
| ------------- | ------------- | ------------- |
| alg | char | Name of algorithm |
| ... | ... | Parameters for specific algorithm | 

For 'Simplified' algorithm, which uses aggregate predictions:

| Field         | format| description |
| ------------- | ------------- | ------------- |
| tx    | double | Period of optimization call (minutes) | 
| ts    | double | Prediction delay | 
| tr    | double | Prediction horizon | 
| bmin  | double | Extra vehicles at nodes | 

#### `Pricing`

| Field         | format| description |
| ------------- | ------------- | ------------- |
| tp            | double | Period of pricing module call (minutes) | 
| movingcostkm  | double | Relocation cost (per km) | 
| basetariffkm  | double | Base tariff (per km) | 
| alternativecostfile   | double | File with alternative mode cost for each trip (see below for details). Leave empty if using `alternativecost` or `alternativecostkm` | 
| alternativecost   | double | Alternative mode cost (for each trip). Leave empty if using `alternativecostkm` | 
| alternativecostkm | double | Alternative mode cost (per km) | 
| traveltimecost   | logical | Consider perceived cost of travel time with SAEV in mode choice? | 
| mintariff     | double | Minimum tariff | 
| VOT           | double | Value of time (per hour) |
| dynamic       | logical | Optimize pricing dynamically?  | 

`alternativecostfile` is the name of a .mat file in the same folder as the trip files, with a variable of type cell named `alternativecost`. Each tripday correspond to an entry in the cell, which contains the alternative cost for each trip. Each array in the cell element must be the same length as the corresponding trip file for that day.

#### `FCR`

| Field         | format| description |
| ------------- | ------------- | ------------- |
| filename      | char | ... | 
| limits        | double [1x2] | ... | 
| contracted        | double | Service contract amount (MW) | 
| fastchargesoc     | double | SoC at which vehicles start charging slower | 
| slowchargeratio   | double | Ratio of slow to fast charging | 
| aggregatechargeratio  | double | Equivalent power exchange visible by charging module aggregator | 

### Data preparation

To run the simulations, there are at least 3 external files needed: 

1. a scenario file
1. a trip file
1. a file with electricity price data

#### Scenario files

Scenario files are .mat Matlab files with variables:

1. `T`, [n x n], travel time between each node in minutes;
1. (optional) `C`, [n x 2], list of coordinates of the nodes;
1. (optional) `D`, [n x n], travel distance between each node in km;
1. (optional) `Clusters`, [n x 1], cluster ID for each node;
1. (optional) `clusterIDs`, [nc x 1], node representing center of each cluster;
1. (optional) `chargingStations`, [ns x 1], list of nodes which have charging stations;

These files are stored in folder `data/scenarios/`.

For variable travel time during the day, `T` can also be given as a struct with fields `traveltime` and `hour` where each entry is a snapshot of the travel time matrix at the specified hour.
The model will calculate the travel time at each interval as linear interpolation between these points.

#### Trip files

Trip files are .mat Matlab files with variables `A` and `Atimes`, representing the origin and destination nodes of each trip and their request time (in minutes after midnight).
These files are stored in folder `data/trips/`.

#### Electricity price data

These are .csv files with the price of electricity at each time interval. 
The first row is reserved for statement about resolution: `Resolution in minutes, 30` means each value represent an interval of 30 minutes.
The second row are the headers, and the third row is the first data point.
If there are two columns, the second column represents the carbon intensity.
These files are stored in folder `data/grid/`.

### `extsave` and `info` variables 

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

## (optional) Advanced setup

The model tends to create large savefiles. 
For this reason, it is better to have the folder with savefiles in a location with enough space separated from the code.
To setup a data folder in another location, you can run:

`setupfolders`

This will prompt you to choose a folder where to save output files and parameters.  

## Requirements

The main required toolbox is the *Optimization Toolbox*. 
Other used toolboxes for minor functions (that can be replaced or ignored) is the *Statistics and Machine Learning Toolbox*.
Some functions run faster with *Parallel Computing Toolbox*.

The model is tested for Matlab 2017b or newer. 
The code can be run with older versions with some minor modifications (mostly explicit matrix expansions).


