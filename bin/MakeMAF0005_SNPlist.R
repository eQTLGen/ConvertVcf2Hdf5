library(data.table)
library(stringr)
library(digest)
library(hash)

path<-"[path to the 1Kg_30x_SNP_info.txt]"

setDTthreads(8)

and <- fread(paste0(path, "/1Kg_30x_SNP_info.txt"), header = FALSE)

and <- and[and$V5 > 0.005 & and$V5 < 0.995, ]

and$ID <- paste0(and$V1, ":", and$V2, "_", and$V3, "_", and$V4)

fwrite(as.data.table(and$ID), 
paste0(path, "/SNPlist_1kG_30X_MAF0005.txt"), 
sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# Reference for HASE:
# col names: ID, bp, str_allele1, str_allele2, CHR, allele1, allele2
chr <- str_replace(and$V1, "chr", "")

ref <- data.table(ID = and$ID, 
bp = and$V2, 
str_allele1 = and$V3, 
str_allele2 = and$V4,
CHR = chr)

fwrite(ref, paste0(path, "/wo_hash.txt"), sep = "\t", quote = FALSE, row.names = FALSE)

# calculate hashes for each allele
#[python]
ref <- fread(paste0(path, "wo_hash.txt"))
hashes <- fread(paste0(path, "Alleles_and_hashes.txt"))
hashes <- hashes[, -1, with = FALSE]

ref <- merge(ref, hashes, by.x = "str_allele1", by.y = "Allele")
colnames(ref)[6] <- "allele1"
ref <- merge(ref, hashes, by.x = "str_allele2", by.y = "Allele")
colnames(ref)[7] <- "allele2"

ref <- ref[, c(3, 4, 2, 1, 5, 6, 7)]
ref <- ref[order(CHR, bp), ]
fwrite(ref, paste0(path, "1000Gp1v3.ref.gz"), sep = "\t", quote = FALSE, row.names = FALSE)
