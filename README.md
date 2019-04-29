# Exhaustive occupancy search

This exhaustive search algorithm explores the crystallographic occupancy & B factor of  ligands. 
The algorithm is designed to be used as a last stage step to confirm occupancy when 
the position of the ligand has been refined. 

We consider the crystal with a ligand bound to be in a superposition of bound and ground states, specifically 
we use the framework detailed in [Pearce, N. M., Krojer, T. & von Delft, F. (2017). Acta Cryst. D73, 256-266.] (https://doi.org/10.1107/S2059798317003412)

For a superposed model we vary the occupancy and B factor of the bound state; usually the bound ligand and any alternative conformations of the protein induced when the ligand fits, and the ground state; the conformations modelled for the protien without any ligand bound. i.e. Bound state occupancy 0.1, Ground state occupancy 0.9. At each value of occupancy we explore possible values of B facotr, currently fixing the B factor to be the same across the whole state. 
 
The project heavily relies on the [computational crystallographic toolbox](https://cci.lbl.gov/cctbx_docs/index.html), as distributed in [CCP4](http://www.ccp4.ac.uk/) as well as tools from [PanDDA](https://pandda.bitbucket.io).

## Installation
Install a recent version of CCP4.

Please install an up to date version of [panddas](https://pandda.bitbucket.io/index.html#install)

Download/clone this code from github. Unzip the download directory and cd into it:

```
cd "/path/to/download/directory" 
```

Then install using

```
ccp4-python setup.py install
```

## Usage

Scripts accessible are:
```
exhaustive
```

This can be run minimally by supplying a pdb and mtz file

```
exhaustive test/resources/FALZA-x0085/FALZA-x0085.free.mtz test/resources/FALZA-x0085/refine.pdb
```

Or with additional parameters:

```
exhaustive test/resources/FALZA-x0085/FALZA-x0085.free.mtz test/resources/FALZA-x0085/refine.pdb exhaustive.options.step=0.2'''
```
Options for this script are detailed in exhaustive.utils.phil, and can be supplied on command line as shown for the step size in the above example. The input is minimally a pdb file refined using giant.quick_refine, and the corresponding mtz file. 

```
exhaustive_multiple_sampling
```

As above but runs at step sizes 0.05 across occupancy [0,1], u_iso [0.2, 1.2] and then at the minima +/- 0.1 in step sizes of 0.01. This will increase the program runtime by 10 times, and is appropraite as in tested cases the minima of the difference maps are flat 


