##################################################################################################
#PACKAGES
##################################################################################################
import pandas as pd
import argparse

##################################################################################################
#ARGPARSE
##################################################################################################
parser= argparse.ArgumentParser()
group= parser.add_mutually_exclusive_group(required= True) #making sure that the user can not use the --mirna, --disease, --sm and --sd in tandem
group.add_argument("--mirna", help= "search the database for all genes and diseases associated with a miRNA") #option 1
group.add_argument("--disease", help= "search the database for all miRNA and genes associated with a disease" ) #option 2
group.add_argument("--sm", help= "given part of a miRNA, search the database for all miRNA that match") #option 3
group.add_argument("--sd", help= "given part of a disease, search the database for all diseases that match") #option 4
parser.add_argument("-geneConf", help="filter gene to disease results based on the confidence score", type= int, choices=range(80,101)) #option 5. used with option 1
parser.add_argument("-locusConf", help="filter mirna to locus results based on the predicted score", type= int, choices=range(80,101)) #option 5'. used with option 2
parser.add_argument("-score", help= "choose whether confidence or predicted scores are added to output", choices=["true","false"], default= "true") #option 6. used with options 1 and 2
args=parser.parse_args()

##################################################################################################
#DATA
##################################################################################################
print("####################")
print("LOADING DATABASE")
print("####################")

#reading the data
data=pd.read_csv("./miRNA_gene_disease.csv")
#removing 1st column (Unnamed: 0)
data.drop('Unnamed: 0', inplace=True, axis=1)

print("####################")
print("DATABASE LOADED")
print("####################")

##################################################################################################
#METHODS
##################################################################################################

#--miRNA option: return the genes and diseases associated with a certain miRNA
#-geneConf: filter for a certain gene confidence score
#if -score is false, do not add score column to the output, else return score column with output
def mirnaFiltering(miRNA, score, geneConf):
    #filtering for a specific mirna
    ans=data[data['miRNA']==miRNA]
    ans=ans[['geneSymbol', 'diseaseName','score']]
    unique_ans=ans.drop_duplicates() #keeping only unique combinations of gene-disease-score (miRNA can bind to more than 1 locus that bind to the same gene)
    #check whether the dataframe is empty. report error message if so
    if len(unique_ans)==0:
        print("There is no miRNA entry in the database named: %s" % miRNA)
        exit()
    #checking geneConf argument for extra filtering
    if geneConf != None:
        filt= unique_ans['score']> geneConf
        filt_ans= unique_ans[filt]
        #check whether the dataframe is empty. report error message if so
        if len(filt_ans)==0:
            print("There are no genes and diseases in the database associated with %s with a confidence score of %s" % (miRNA, geneConf))
            exit()
        #checking score argument in tandem with geneConf
        if score =="false":
            print(filt_ans[["geneSymbol","diseaseName"]].to_markdown(tablefmt="grid", index=False))
        else:
            print(filt_ans[["geneSymbol","diseaseName","score"]].to_markdown(tablefmt="grid", index=False))
    else:
        #checking score argument (geneConf argument not used)
        if score =="false":
            print(unique_ans[["geneSymbol","diseaseName"]].to_markdown(tablefmt="grid", index=False))
        else:
            print(unique_ans[["geneSymbol","diseaseName","score"]].to_markdown(tablefmt="grid", index=False))

#--disease option: returning the mirna and genes associated with a certain disease
#-locusConf: filter for a certain mirna to locus confidence score
#if -score is false, do not add score column to the output, else return score column with output
def diseaseFiltering(disease, score, locusConf):
    #filtering for a specific disease
    ans=data[data['diseaseName']==disease.lower()]
    ans=ans[['geneSymbol', 'miRNA','miRNAPredScore']]
    unique_ans= ans.drop_duplicates() #duplicates may arise from same miRNA, different locus, same miRNAPredScore
    #checking locusConf argument for extra filtering
    #check whether the dataframe is empty. report error message if so
    if len(unique_ans)==0:
        print("There is no disease entry in the database named: %s" % disease)
        exit()
    #checking locusConf argument for extra filtering
    if locusConf != None:
        filt=unique_ans['miRNAPredScore']>locusConf
        filt_ans=unique_ans[filt]
        #check whether the dataframe is empty. report error message if so
        if len(filt_ans)==0:
            print("There are no miRNA and genes in the database associated with %s with a predicted score of %s" % (disease, locusConf))
            exit()
        #checking score argument in tandem with locusConf
        if score =="false":
            print(filt_ans[["miRNA","geneSymbol"]].to_markdown(tablefmt="grid", index=False))
        else:
            print(filt_ans[["miRNA","geneSymbol","miRNAPredScore"]].to_markdown(tablefmt="grid", index=False))
    else:
        #checking score argument
        if score =="false":
            print(unique_ans[["miRNA", "geneSymbol"]].to_markdown(tablefmt="grid", index=False))
        else:
            print(unique_ans[["miRNA","geneSymbol","miRNAPredScore"]].to_markdown(tablefmt="grid", index=False))

#searches the database for miRNA containing part of a miRNA(input)
def mirnaSimilar(mir):
    ans=data[data['miRNA'].str.count(mir)>0]
    unique_ans=ans['miRNA'].unique()
    #check whether the dataframe is empty. report error message if so
    if len(unique_ans)!=0:
        print(unique_ans.tolist())
    else:
        print("There are no miRNA entries in the database that contain: %s" % mir)

#searches the database for disease names containing part of the name(input)
def diseaseSimilar(dis):
    ans=data[data['diseaseName'].str.count(dis)>0]
    unique_ans=ans['diseaseName'].unique()
    #check whether the dataframe is empty. report error message if so
    if len(unique_ans)!=0:
        print(unique_ans.tolist())
    else:
        print("There are no disease entries in the database that contain: %s" % dis)

#making sure that no incorrect argument combination is used by the user
def argumentChecker():
    if args.mirna!=None and args.locusConf!=None:
        print("These arguments cannot be used in tandem")
        exit()
    if args.disease!=None and args.geneConf!=None:
        print("These arguments cannot be used in tandem")
        exit()
    if args.sm!=None and (args.locusConf!=None or args.geneConf!=None or args.score=="false"):
        print("These arguments cannot be used in tandem")
        exit()
    if args.sd!=None and (args.locusConf!=None or args.geneConf!=None or args.score=="false"):
        print("These arguments cannot be used in tandem")
        exit()

#checks arguments and call appropriate functions
def progExecution():
    if args.mirna!=None:
        mirnaFiltering(args.mirna, args.score, args.geneConf)
    elif args.disease!=None: #can use "" on terminal to capture disease name with more than 1 word
        diseaseFiltering(args.disease, args.score, args.locusConf)
    elif args.sm!=None:
        mirnaSimilar(args.sm)
    elif args.sd!=None:
        diseaseSimilar(args.sd)

##################################################################################################
#MAIN
##################################################################################################
argumentChecker()
progExecution()