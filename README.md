# Exhaustive Search

Exhaustive search algorithm explore occupancy & B factor of crystallographic ligands. 
The algorithm is designed to be used as a last stage step to confirm occupancy when 
the position of the ligand has been refined. 

Exhaustive search implies the exploration of a larger number of possible values 
of occupancy and B factor than in conventional refinement algorithms, 
such as phenix.refine and REFMAC5. 

We consider the crystal with a ligand bound to be in a 
superposition of bound and ground states, specifically 
we use the framework used in the PanDDA algorithm.
 
The project heavily relies on CCTBX libraries,as well as tools from the PanDDA repository.

## Installation
Install a recent version of CCP4.

Please install an up to date version of panddas, as per below instructions
https://pandda.bitbucket.io/index.html#install

Download from here. Unzip the download directory and cd into it:

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
Options for this script are detailed in exhaustive.utils.phil, and can be supplied on command line as shown for the step size in the above example. The input should be a pdb file refined using giant.quick_refine, and the corresponding mtz file. 

```
exhaustive_multiple_sampling
```

As above but runs at step sizes 0.05 across occupancy [0,1], u_iso [0.2, 1.2] and then at the minima +/- 0.1 in step sizes of 0.01. This will increase the program runtime by 10 times, and is appropraite as in tested cases the minima of the difference maps are flat 


