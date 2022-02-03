library(data.table)
library(stringr)

path<-"[path to the 1Kg_30x_SNP_info.txt]"

setDTthreads(8)

and <- fread(paste0(path, "/1Kg_30x_SNP_info.txt"), header = FALSE)

and <- and[and$V5 > 0.005 & and$V5 < 0.995, ]

and$ID <- paste0(and$V1, ":", and$V2, "_", and$V3, "_", and$V4)

fwrite(as.data.table(and$ID), 
paste0(path, "/SNPlist_1kG_30X_MAF0005.txt"), 
sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
