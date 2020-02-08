# correlation plot
library("ggpubr")

CorDf=function(Df,Target_x_Col,Target_y_Cols)
  {
  ResultDf=as.data.frame(Target_y_Cols)
  for(i in c(1:length(Target_y_Cols)))
      {
    temp=cor.test(Df[,Target_x_Col],Df[,Target_y_Cols[i]],method = 'pearson')
    ResultDf$r[i]=temp$estimate
    ResultDf$p[i]=temp$p.value
    }
  ResultDf$adjp=p.adjust(ResultDf$p,method = 'BH')
  return(ResultDf)
}
for( i in TargetGenes){
  CorDf_test=CorDf(Df=TCGA_SeqDfUse1,
                   Target_x_Col=i,
                   Target_y_Cols=colnames(TCGA_SeqDfUse1)[
                     colnames(TCGA_SeqDfUse1)%in%c(FAO_Genes,FAS_Genes)])
  CorDf_test$Group='FAO'
  CorDf_test$Group[CorDf_test$Target_y_Cols%in%FAS_Genes]='FAS'
  write.csv(CorDf_test,paste0('./Plot20190626/TCGA_Correlation_',i,'.csv'))
}

for( i in TargetGenes){
  CorDf_test=CorDf(Df=CGGA_SeqDf1Use,
                   Target_x_Col=i,
                   Target_y_Cols=colnames(CGGA_SeqDf1Use)[
                     colnames(CGGA_SeqDf1Use)%in%c(FAO_Genes,FAS_Genes)])
  CorDf_test$Group='FAO'
  CorDf_test$Group[CorDf_test$Target_y_Cols%in%FAS_Genes]='FAS'
  write.csv(CorDf_test,paste0('./Plot20190626/CGGA_Correlation_',i,'.csv'))
}



