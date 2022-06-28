#reading data
locus_to_gene=read.csv2("locus-to-geneMOD.csv", sep=",")

#adding column names
colnames(locus_to_gene)=c("locus", "organism", "geneSymbol")

#editing dataframe column type
locus_to_gene$locus = as.character(locus_to_gene$locus) #the type was integer
locus_to_gene$organism = as.character(locus_to_gene$organism) #the type was integer
locus_to_gene$geneSymbol = as.character(locus_to_gene$geneSymbol) #the type was integer

#finding all homo sapiens inlcuding the typos
#some of the typos have a 0 instead of the 1st o
homo=grep(pattern="(.*(H|h)(O|o|0)(M|m)(O|o|0).*)", perl=TRUE,x = locus_to_gene$organism)
homo_locus_gene=locus_to_gene[homo,]
rm(homo, locus_to_gene)

#finding/fixing the homo sapiens typos 
typo=which(homo_locus_gene$organism!="Homo sapiens")
for (i in typo)
{
  homo_locus_gene$organism[i]="Homo sapiens"
}
rm(i,typo)

#checking for leading or trailing whitespaces
grep(pattern=" ", perl=TRUE,x = homo_locus_gene$locus)
grep(pattern=" ", perl=TRUE,x = homo_locus_gene$geneSymbol)

#using which(is.na()) for each column there is no missing data
which(is.na(homo_locus_gene$locus))
which(is.na(homo_locus_gene$organism))
which(is.na(homo_locus_gene$geneSymbol))

#using sum(duplicated(dataframe)) -> there are 0 duplicated lines
sum(duplicated(homo_locus_gene))

#comparing locuses from locus-to-geneMOD.csv to locuses from miRNA-locusMOD.txt
#homo_locus=unique(homo_locus_gene$locus)
#miRNA_locus=unique(miRNA_to_locus$locus)
#intersect(homo_locus,miRNA_locus) #returns everything in common in both 
#setdiff(homo_locus,miRNA_locus) #returns the locuses in homo_locus that are not present in miRNA_locus
#rm(homo_locus, miRNA_locus)

#outputting data frame
write.csv(homo_locus_gene, "../locus-to-gene-filtered.csv",row.names = T)
