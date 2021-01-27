# Script by Davide Pirolli
# usage: Rscript PCA_analysis.R 
# tested with R 3.6.0
# Read instruction below
# if you use it in your paper, I would appreciate if you could cite the repository
#
library(bio3d)
#  non-interactive
cat("coordinates pdb file: \n");
filecoord <- readLines("stdin",n=1);
cat("trajectory dcd file: \n");
filetraj <- readLines("stdin",n=1); 
cat("give me an output prefix: \n");
outprefix <- readLines("stdin",n=1);

pdb <- read.pdb(filecoord)
dcd <- read.dcd(filetraj)

ca.inds <- atom.select(pdb, elety="CA")
xyz <- fit.xyz(fixed=pdb$xyz, mobile=dcd, fixed.inds=ca.inds$xyz, mobile.inds=ca.inds$xyz)

dim(xyz) == dim(dcd)

pc <- pca.xyz(xyz[,ca.inds$xyz])

PoV <- pc$sdev^2/sum(pc$sdev^2)
write.table(PoV, file=sprintf("%s_PoV.csv", outprefix))

pdfname1 <- paste(outprefix, "_PCA_glob.pdf", col="", sep="")
pdf(pdfname1) 
plot(pc, col=bwr.colors(nrow(xyz)) )
dev.off()

hc <- hclust(dist(pc$z[,1:2]))
grps <- cutree(hc, k=2)

pdfname2 <- paste(outprefix,"_PCA_glob_cluster.pdf", col="", sep="")
pdf(pdfname2)
plot(pc, col=grps)
dev.off()

write.table(pc$z[,1], file=sprintf("%s_PCA1.csv", outprefix))
write.table(pc$z[,2], file=sprintf("%s_PCA2.csv", outprefix))
write.table(pc$z[,3], file=sprintf("%s_PCA3.csv", outprefix))

pdfname3 <- paste0(outprefix,"_PCA_Res.pdf")
pdf(pdfname3)
plot.bio3d(pc$au[,1], ylab="PC1 (A)", xlab="Residue Position", typ="l")
points(pc$au[,2], typ="l", col="blue")
points(pc$au[,3], typ="l", col="red")
dev.off()


write.table(pc$au[,1], file=sprintf("%s_PC1_au_resposition.csv", outprefix))
write.table(pc$au[,2], file=sprintf("%s_PC2_au_resposition.csv", outprefix))
write.table(pc$au[,3], file=sprintf("%s_PC3_au_resposition.csv", outprefix))


p1 <- mktrj.pca(pc, pc=1,b=pc$au[,1], file=sprintf("%s_pc1.pdb", outprefix))
p2 <- mktrj.pca(pc, pc=2,b=pc$au[,2], file=sprintf("%s_pc2.pdb", outprefix))
p3 <- mktrj.pca(pc, pc=2,b=pc$au[,3], file=sprintf("%s_pc3.pdb", outprefix))

eigv <- pc$L
trace <- sum(eigv)
write.table(pc$L, file=sprintf("%s_eigenvalues.csv", outprefix))
print()
print(The trace value is:  )
trace
print(Angstrom^2)
print ()
proc.time()
