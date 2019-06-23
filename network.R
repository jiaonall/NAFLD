library(Hmisc)

### reform of correlation matrix
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}
NAFLD_sparcc_cor <- read.table("sparcc_NAFLD.txt",sep = '\t',fill = T,header = T,row.names = 1)
NAFLD_sparcc_p <- read.table("pvals_NAFLD.txt",sep = '\t',fill = T,header = T,row.names = 1)
NAFLD_sparcc <- flattenCorrMatrix(NAFLD_sparcc_cor,NAFLD_sparcc_p)
NAFLD_sparcc_sig <-NAFLD_sparcc[which(NAFLD_sparcc$p<0.05),]
NAFLD_sparcc_sig_high <-NAFLD_sparcc_sig[which(abs(NAFLD_sparcc_sig$cor)>0.4),]
write.table(NAFLD_sparcc_sig,'NAFLD_sparcc_sig.txt', sep = '\t',col.names= T, row.names = T)
write.table(NAFLD_sparcc_sig_high,'NAFLD_sparcc_sig_strong.txt', sep = '\t',col.names= T, row.names = T)

CRC_sparcc_cor <- read.table("sparcc_CRC.txt",sep = '\t',fill = T,header = T,row.names = 1)
CRC_sparcc_p <- read.table("pvals_CRC.txt",sep = '\t',fill = T,header = T,row.names = 1)
CRC_sparcc <- flattenCorrMatrix(CRC_sparcc_cor,CRC_sparcc_p)
CRC_sparcc_sig <-CRC_sparcc[which(CRC_sparcc$p<0.05),]
CRC_sparcc_sig_high <-CRC_sparcc_sig[which(abs(CRC_sparcc_sig$cor)>0.2),]
write.table(CRC_sparcc_sig,'CRC_sparcc_sig.txt', sep = '\t',col.names= T, row.names = T)
write.table(CRC_sparcc_sig_high,'CRC_sparcc_sig_strong.txt', sep = '\t',col.names= T, row.names = T)
