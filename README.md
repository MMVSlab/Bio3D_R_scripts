# Bio3D_R_scripts

SCRIPT: PCA_analysis.R (by Davide Pirolli)

This script performs an automated PCA of your molecular dynamics simulation, employing the bio3d library, 
as described in the GrantLab website (http://thegrantlab.org/bio3d).
I found this implementation of bio3d very useful when working with Desmond-generated trajectories 
(that I convert in pdb/dcd by vmd, that is fast for managing with trajectories).

On Centos7 I had some prerequisites:
  1) libxml-devel
  2) openssl-devel
  3) motif-devel

Requisites: 
  1) The bio3d library installed in the R environment 
     [just start R and type: install.packages("bio3d", dependencies=TRUE) ]
  2) A pdb with the starting coordinates of the simulation 
  3) A dcd with the trajectory

Usage: Rscript PCA_analysis.R (and you will be asked the pdb, the dcd names and a prefix for output files)


The script will generate csv files importable into excel of:

1) The projection of protein motion onto the first three principal components (output_PCA{1,2,3}.csv)
2) The contribution of each residue to each of the first principal components (output_PC{1,2,3}_au_resposition.csv)
3) The trajectory of the first three principal components (pdb files)
4) The eigenvalues of each eigenvector in Angstrom^2 (output_eigenvalues.csv)
5) The TRACE value in Angstrom^2 (printed on screen, 'cause it's just the sum of eigenvalues)

SCRIPT: RMSIP.R

Briefly it calculates the overlap between the conformational space described by the first 10 eigenvectors.
It is commoly used to measure the similarity between conformational spaces explored by the first ten principal 
components of two different parts of a trajectory (es. PROTEINS: Structure, Function, and Genetics 36:419–424 (1999)).
RMSIP is also used to compare the conformational space of two different simulations 
of the same protein in different conditions (doi:10.1371/journal.pone.0121114) (but is not (yet) the case of this script).

All you have to do is provide a pdb, a dcd. Once you have given the number of the starting, intermediate and last frames,
the script divides the whole trajectory into two parts
part 1: starting frame -> intermediate frame
part 2: intermediate frame -> last frame
and compare the best 10 principal components of these two parts. The output is the Root Mean Square Inner Product, which goes from 0 to 1
and the more is near to 1 the more the two parts are identical.
