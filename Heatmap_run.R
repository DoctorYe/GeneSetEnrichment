# Heatmap run
dim(TCGA_SeqDfUse1)
#[1] 171 303
dim(CGGA_SeqDf1Use)
# [1] 693 293

# TCGA
colnames(TCGA_SeqDfUse1)=TCGA_SeqDf1$GeneName[TCGA_SeqDf1$GeneID%in%colnames(TCGA_SeqDfUse1)]
TCGAHeatmap_prepare=HeatmapPrepare(Df=TCGA_SeqDfUse1,
                TargetGene='CD47',
                GeneSet_List=list(FAO=FAO_Genes,
                                  FAS=FAS_Genes),
                CutMethod=mean)
DrawHeatmap(Df = t(TCGAHeatmap_prepare$Df),
            ColumnGroup = TCGAHeatmap_prepare$RowGroup,
            ColumnColorList=list(Groups = c("CD47(+)" = "red", "CD47(-)" = "green")),
            RowGroup = TCGAHeatmap_prepare$ColumnGroup,
            SavePNGPath ='./TCGAHeatmapSeq_CD47.png' )  
tempfun=function(Df,TargetGene,GeneSet_List,SavePNGPath){
  Heatmap_prepare=HeatmapPrepare(Df=Df,
                                 TargetGene=TargetGene,
                                 GeneSet_List=GeneSet_List,
                                 CutMethod=mean)
  temp_g=c("red","green")
  names(temp_g)=c(1,2)
  DrawHeatmap(Df = t(Heatmap_prepare$Df),
              ColumnGroup = Heatmap_prepare$RowGroup,
              ColumnColorList=list(Groups=temp_g),
              RowGroup = Heatmap_prepare$ColumnGroup,
              SavePNGPath =SavePNGPath)  
  
}

tempfun(Df = TCGA_SeqDfUse1,
        TargetGene = 'CD274',
        GeneSet_List = list(FAO=FAO_Genes,
                            FAS=FAS_Genes),
        SavePNGPath='./TCGAHeatmapSeq_CD274.png')


# CGGA

tempfun(Df = CGGA_SeqDf1Use,
        TargetGene = 'CD47',
        GeneSet_List = list(FAO=FAO_Genes,
                            FAS=FAS_Genes),
        SavePNGPath='./CGGAHeatmapSeq_CD47.png')

tempfun(Df = CGGA_SeqDf1Use,
        TargetGene = 'CD274',
        GeneSet_List = list(FAO=FAO_Genes,
                            FAS=FAS_Genes),
        SavePNGPath='./CGGAHeatmapSeq_CD274.png')



