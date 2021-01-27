# Usefull script to measure the similarity of the conformational spaces captured between two parts of a simulations
# Commonly used to check the simulation convergence
# Script by Davide Pirolli
# usage: Rscript RMSIP.R 
# tested with R 3.6.0
# Read instruction below
# if you use it in your paper, I would appreciate if you could cite the repository
#
library(bio3d)

cat("coordinates pdb file: \n");
filecoord <- readLines("stdin",n=1);
cat("trajectory dcd file: \n");
filetraj <- readLines("stdin",n=1);
cat("")
cat("The simulation will be divided into two subsequent parts.\n");
cat("insert starting frame: \n");
startframe <- readLines("stdin",n=1);
cat("insert starting intermediate frame: \n");
interframe <- readLines("stdin",n=1);
cat("insert last frame: \n");
lastframe <- readLines("stdin",n=1);


pdb <- read.pdb(filecoord)
dcd <- read.dcd(filetraj)

trj.sub <- dcd[startframe:interframe, ] # Subpart 1 of the simulation
 
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=trj.sub, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
pc1 <- pca.xyz(xyz[,ca.inds$xyz])

trj2.sub <- dcd[interframe:lastframe, ] # subpart 2 of MD simulation
ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=trj2.sub, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)
pc2 <- pca.xyz(xyz[,ca.inds$xyz])
 
 
# Calculate the RMSIP between best 10 PCs of the two parts of the simulation
r <- rmsip(pc1, pc2, subset=10, row.name="sim_1-PC", col.name="sim_2-PC")

# Plot pairwise overlap values
plot(r, xlab="trajectory_1-PC", ylab="trajectory_2-PC")
 
print(r$rmsip)
