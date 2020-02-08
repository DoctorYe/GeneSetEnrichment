# Path Enrichment Run
setwd('/Users/ningrongye/Desktop/python/Gene/JN/')



TCGA_Inform=read.csv('/Users/ningrongye/Desktop/python/Gene/GBM_clinicalMatrix',sep='\t')
TCGA_Inform_GBM=TCGA_Inform[TCGA_Inform$histological_type%in%
                                levels(TCGA_Inform$histological_type)[c(2,4)],]


TCGA_SeqDfPath='./gdac.broadinstitute.org_GBM.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.2016012800.0.0/GBM.rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.data.txt'
TCGA_ArrayDfPath='/Users/ningrongye/Desktop/python/Gene/HT_HG-U133A'
CGGA_SeqDfPath='./CGGA.mRNAseq_693.RSEM-genes.20190531.txt'
FAO_Genes=levels(read.csv('./fatty acid oxidation gene.csv')$gene.name)
FAS_Genes=levels(read.csv('./fatty acid synthesis gene.csv')$gene.name)

TargetGenes=c('CD274', 'CD47') # CD274 not exist

TCGA_SeqDf=read.csv(TCGA_SeqDfPath,  sep = '\t')
TCGA_ArrayDf=read.csv(TCGA_ArrayDfPath,  sep = '\t')
CGGA_SeqDf=read.csv(CGGA_SeqDfPath,  sep = '\t')



#### seq data TCGA
TCGA_SeqDf1=TCGA_SeqDf[-1,]
for( i in c(2:length(colnames(TCGA_SeqDf)))){
  TCGA_SeqDf1[,i]=as.numeric(as.character(TCGA_SeqDf1[,i]))
}

TCGA_SeqDf1[,1]=as.character(TCGA_SeqDf1[,1])
for(i in c(1:length(TCGA_SeqDf1[,1]))){
  temp=strsplit(TCGA_SeqDf1[i,1],'\\|')[[1]]
  TCGA_SeqDf1$GeneName[i]=temp[1]
  TCGA_SeqDf1$GeneID[i]=temp[2]
}

DfCharacter_Col=c('GeneName','GeneID','Hybridization.REF')

row.names(TCGA_SeqDf1)=TCGA_SeqDf1$GeneID
TCGA_SeqDf2=Df_Clean(Df=TCGA_SeqDf1,CharacterCol = DfCharacter_Col)
Df_CUse=gsub("\\.",'-',rownames(TCGA_SeqDf2))
TCGA_SeqDf2=as.data.frame(TCGA_SeqDf2)
for(i in c(1:length(rownames(TCGA_SeqDf2)))){
  TCGA_SeqDf2$PatientID[i]=paste0(c(strsplit(rownames(TCGA_SeqDf2)[i],
                                             '\\.')[[1]][1:3],
                                    '01'),
                                  collapse = '-')
  
}
rownames(TCGA_SeqDf2)=TCGA_SeqDf2$PatientID
TCGA_SeqDf3=TCGA_SeqDf2[,!colnames(TCGA_SeqDf2)%in%c('PatientID')]
dim(TCGA_SeqDf3)
ColumnNameNeeded=TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName%in%c(FAO_Genes,
                                                              FAS_Genes,
                                                              TargetGenes)]

FAO_GenesID=TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName%in%FAO_Genes]
FAS_GenesID=TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName%in%FAS_Genes]
TargetGenesID=TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName%in%TargetGenes]

TCGA_SeqDfUse1=TCGA_SeqDf3[,colnames(TCGA_SeqDf3)%in%ColumnNameNeeded]

# FAO vs FAS
DiffExpression(Df = TCGA_SeqDfUse1,GeneSetNeeded = FAO_GenesID,TargetGene = TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName=='CD47'])
# 0.7082784 0.6686831
DiffExpression(Df = TCGA_SeqDfUse1,GeneSetNeeded = FAO_GenesID,TargetGene = TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName=='CD274')

# 0.3673433 0.8802031

# FAO vs all
DiffExpression(Df = TCGA_SeqDf3,GeneSetNeeded = FAO_GenesID,
               TargetGene = TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName=='CD47'])
#0.0008751896 0.9998650971
DiffExpression(Df = TCGA_SeqDf3,GeneSetNeeded = FAO_GenesID,
               TargetGene = TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName=='CD47'],Plot='./test.png')

#
#   GeneSetNeeded
#Dif     0     1
#0 17443   146
#1  1410     2
#pathway      pval      padj        ES      NES nMoreExtreme size
#1:    Path 0.1262411 0.1262411 0.2491839 1.200955          889  148
#leadingEdge
#1: 1892,549,117145,36,23659,27163,...


DiffExpression(Df = TCGA_SeqDf3,GeneSetNeeded = FAO_GenesID,
               TargetGene = TCGA_SeqDf1$GeneID[TCGA_SeqDf1$GeneName=='CD274'])

#0.09472672 0.96542632



#### array data TCGA
# in 20190623.R





#### seq data CGGA

rownames(CGGA_SeqDf)=as.character(CGGA_SeqDf$Gene_Name)
CGGA_SeqDf1=Df_Clean(Df = CGGA_SeqDf,CharacterCol = 'Gene_Name')
any(colnames(CGGA_SeqDf1)=='CD47')
#TRUE

any(colnames(CGGA_SeqDf1)=='CD274')
#TRUE

# FAO vs FAS
# CD47
CGGA_SeqDf1=as.data.frame(CGGA_SeqDf1)
CGGA_SeqDf1Use=CGGA_SeqDf1[,colnames(CGGA_SeqDf1)%in%c(FAO_Genes,
                                        FAS_Genes,
                                        TargetGenes)]
DiffExpression(Df = CGGA_SeqDf1Use,
               GeneSetNeeded =FAO_Genes,
               TargetGene ="CD47")
#0.4995138
# CD274
DiffExpression(Df = CGGA_SeqDf1Use,
               GeneSetNeeded =FAO_Genes,
               TargetGene ="CD274")
# 0.2169626


# FAO vs all
# CD47
DiffExpression(Df = CGGA_SeqDf1,
               GeneSetNeeded =FAO_Genes,
               TargetGene ="CD47")
# fisher p=0.01146978
DiffExpression(Df = CGGA_SeqDf1,
               GeneSetNeeded =FAO_Genes,
               TargetGene ="CD47",
               Plot = './CGGA_FAO_CD47.png')
# fgsea p: 0.1073107  

# CD274

DiffExpression(Df = CGGA_SeqDf1,
               GeneSetNeeded =FAO_Genes,
               TargetGene ="CD274")
#  0.1035727




