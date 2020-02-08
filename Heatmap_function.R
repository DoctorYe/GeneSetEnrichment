# Plot heatmap 
# the script was run after PathEnrichment_run.R

library(ComplexHeatmap)
HeatmapPrepare=function(Df,TargetGene,GeneSet_List,CutMethod=mean){
  Df$Group=paste0(TargetGene,'(-)')
  Df$Group[Df[,TargetGene]>=CutMethod(Df[,TargetGene])]=paste0(TargetGene,'(+)')
  RowGroup=Df$Group
  GeneNeeded=c()
  ColumnGroups_temp=colnames(Df)
  for(i in names(GeneSet_List)){
    temp=as.vector(GeneSet_List[[i]])
    ColumnGroups_temp[ColumnGroups_temp%in%temp]=i
  }
  DfUse=Df[,ColumnGroups_temp%in%names(GeneSet_List)]
  ColumnGroups=ColumnGroups_temp[ColumnGroups_temp%in%names(GeneSet_List)]
  Result=list(ColumnGroup=ColumnGroups,RowGroup=RowGroup,Df=DfUse)
  return(Result)
}
#library(circlize)
#col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
#col_fun(seq(-3, 3))
DrawHeatmap=function(Df,
         DfValueName='Zscore',
         ColumnGroup=c('a','a','b','b'),
         ColumnColorList=list(Groups = c("a" = "red", "b" = "green")),
         RowGroup=c('1','2'),
         SavePNGPath='./Heatmap.png'){
  
  ha = HeatmapAnnotation(bar = ColumnGroup,name = 'Groups',
                         col = ColumnColorList)
  #col_fun = colorRamp2(c(-2, 0, 2), c("green", "white", "red"))
  #col_fun(seq(-3, 3))
  ht=Heatmap(Df, name = DfValueName, 
          show_column_dend = FALSE,show_row_dend = FALSE,
          row_split =RowGroup,column_split = ColumnGroup,
          top_annotation = ha,
          show_row_names = FALSE,
          show_column_names = FALSE,#col = col_fun,
          width = unit(20, "cm"), height = unit(20, "cm"))
  #ColumnGroup_L = Legend(col = ColumnColorList, 
  #                       at = c(0:length(ColumnColorList$bar)-1), 
  #                    labels = names(ColumnColorList$bar))
  png(file=SavePNGPath,width =30,height = 30,units = "cm", res  = 1200)
  draw(ht)
  dev.off()
}

test$RowGroup


