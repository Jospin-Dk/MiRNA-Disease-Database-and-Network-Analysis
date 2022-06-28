#reading data
gene_to_disease= read.delim2("curated_gene_disease_associationsMOD.tsv")

#filtering
#1- keeping columns of interest
gene_to_disease=subset(gene_to_disease, select = c(geneSymbol,diseaseName, score))

#editing dataframe column type
gene_to_disease$geneSymbol = as.character(gene_to_disease$geneSymbol) #the type was integer
gene_to_disease$diseaseName = as.character(gene_to_disease$diseaseName) #the type was integer
gene_to_disease$score= as.numeric(as.character(gene_to_disease$score)) #the type was integer

#2- keeping gene/disease confidence score of >0.8 
temp=which(gene_to_disease$score<=0.8)
gene_to_disease=gene_to_disease[-temp,]
rm(temp)

#checking for leading or trailing spaces
grep(pattern=" ", perl=TRUE,x = gene_to_disease$geneSymbol)

#making the score column /100
gene_to_disease$score=gene_to_disease$score*100

#lower case everything in the diseaseName column
gene_to_disease$diseaseName=tolower(gene_to_disease$diseaseName)

#using which(is.na()) for each column there is no missing data
which(is.na(gene_to_disease$geneSymbol))
which(is.na(gene_to_disease$diseaseName))
which(is.na(gene_to_disease$score))

#using sum(duplicated(dataframe)) -> there are 0 duplicated lines
sum(duplicated(gene_to_disease))
test=unique(gene_to_disease$diseaseName)


#outputting data frame
write.csv(gene_to_disease, "../curated_gene_disease_associations_filtered.csv",row.names = T)
