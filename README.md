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
  
  