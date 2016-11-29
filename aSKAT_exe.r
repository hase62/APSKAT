#install.packages("SKAT")
library(SKAT)
source("aSKAT_func.r")

###Sample Run For All Genes Using bfile###
result<-aSKAT_bfile(prefix = "hapmap3_r3_b36_fwd.consensus.qc.poly.remove.ch10",
                    SNPSetID = "snpset.txt", significance_level = 0.05/20000, 
                    ID = "job" , n.maxres = 5 / (0.05/20000))

###Write Results###
write.table(result, paste(commandArgs(TRUE)[1], "-", ID, ".txt", sep=""), 
    row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")

###Sample###
if(FALSE){
  data(SKATBinary.example)
  attach(SKATBinary.example)
  result<-aSKAT_matrix(Z, y, significance_level = 0.05/20000, 
                       ID = "job" , n.maxres = 5 / (0.05/20000))
}