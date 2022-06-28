#reading data
miRNA_to_locus=read.delim2("miRNA-locusMOD.txt", header = F) #file does not have a header

#adding column names
colnames(miRNA_to_locus)=c("number", "miRNA", "locus", "miRNAPredScore")

#editing dataframe column type
miRNA_to_locus$miRNA = as.character(miRNA_to_locus$miRNA) #the type was integer
miRNA_to_locus$locus = as.character(miRNA_to_locus$locus) #the type was integer
miRNA_to_locus$miRNAPredScore= as.numeric(as.character(miRNA_to_locus$miRNAPredScore)) #the type was integer

#filtering
#1- removing the "number" column (no need for it)
miRNA_to_locus=subset(miRNA_to_locus, select = -c(number))
#2- keeping miRNA/locus confidence score of >80 
temp=which(miRNA_to_locus$miRNAPredScore<=80.0)
miRNA_to_locus=miRNA_to_locus[-temp,]
rm(temp)
#3- checking/fixing miss-spelling (miR notation) in miRNA column
typo= grep(pattern="(.*miR.*)|(.*let.*)", invert = TRUE, perl=TRUE,x = miRNA_to_locus$miRNA)
for(x in typo)
{
  temp=unlist(strsplit(miRNA_to_locus$miRNA[x], split="-"))
  final=paste(temp[1],"miR",sep="-")
  for (i in 3:length(temp))
  {
    final=paste(final,temp[i],sep="-")
  }
  miRNA_to_locus$miRNA[x]=final
}
rm(x,i,temp,typo,final)
#4- checking miss-spelling in microRNA column pre-miR
typo= grep(pattern="(cfa.*)|(mmu.*)|(gga.*)|(hsa.*)|(rno.*)", invert = TRUE, perl=TRUE,x = miRNA_to_locus$miRNA)
for (x in typo) 
{
  temp= unlist(strsplit(miRNA_to_locus$miRNA[x], split="-"))
  if(tolower(temp[1])=="cfa")
  {
    final="cfa"
    for (i in 2:length(temp)) 
    {
      final=paste(final,temp[i],sep = "-")
    }
    miRNA_to_locus$miRNA[x]=final
  }
  else if(tolower(temp[1])=="mmu")
  {
    final="mmu"
    for (i in 2:length(temp)) 
    {
      final=paste(final,temp[i],sep = "-")
    } 
    miRNA_to_locus$miRNA[x]=final
  }
  else if (tolower(temp[1])=="gga")
  {
    final="gga"
    for (i in 2:length(temp)) 
    {
      final=paste(final,temp[i],sep = "-")
    }
    miRNA_to_locus$miRNA[x]=final
  }
  else if (tolower(temp[1])=="hsa")
  {
    final="hsa"
    for (i in 2:length(temp)) 
    {
      final=paste(final,temp[i],sep = "-")
    }
    miRNA_to_locus$miRNA[x]=final
  }
  else
  {
    final="rno"
    for (i in 2:length(temp)) 
    {
      final=paste(final,temp[i],sep = "-")
    }
    miRNA_to_locus$miRNA[x]=final
  }
}
rm(final,i,x,temp,typo)
#5- trimming leading and trailing white space from the miRNA and locus columns 
miRNA_to_locus$miRNA= trimws(miRNA_to_locus$miRNA)
miRNA_to_locus$locus= trimws(miRNA_to_locus$locus)

#using which(is.na()) for each column there is no missing data
which(is.na(miRNA_to_locus$miRNA))
which(is.na(miRNA_to_locus$locus))
which(is.na(miRNA_to_locus$miRNAPredScore))
#using sum(duplicated(dataframe)) -> there are 0 duplicated lines
sum(duplicated(miRNA_to_locus[, c("miRNA","locus")]))

#outputting data frame
write.csv(miRNA_to_locus, "./miRNA-locus-filtered.csv",row.names = T)
