RemoveNa=function(Df){
  temp=c()
  for(i in c(1:dim(Df)[2])){
    temp[i]=all(is.na(Df[,i]))
  }
  #print(temp)
  return(Df[,!temp]) 
}

# Path enrichment Function

# remove low variance
library(caret)
library(data.table)
library(ggplot2)
library(ggpubr)
RemoveLowVar=function(Df){
  # col
  NeedRemove=nearZeroVar(Df)
  if(length(NeedRemove)!=0){
    Df1=Df[,-NeedRemove]
  }else{
    Df1=Df
  }
  # row
  Df_t=t(Df1)
  NeedRemoveR=nearZeroVar(Df_t)
  if(length(NeedRemoveR)!=0){
    Df2=Df_t[-NeedRemove,]
  }else{
    Df2=Df_t
  }
  return(t(Df2))
}

Df_Ztrans=function(Df){
  Df1=scale(Df, center = TRUE, scale = TRUE)
  return(Df1)
}

Df_Clean=function(Df,CharacterCol){
  Df0=RemoveNa(Df)
  Df1=t(Df0[,!colnames(Df0)%in%CharacterCol])
  Df2=RemoveLowVar(Df=Df1)
  Df3=as.data.frame(Df_Ztrans(Df2))
  return(Df3)
}

library(limma)
library(fgsea)
DiffExpression=function(Df,GeneSetNeeded, TargetGene,
                        CutMethod=mean,AdjMethod="BH",
                        Plot=FALSE){
  Df$Group=0
  Df$Group[Df[,TargetGene]>=CutMethod(Df[,TargetGene])]=1
  design <- model.matrix(~0+Group, data=Df)
  fit=lmFit(t(Df), design)
  fit_eBayes=eBayes(fit)
  ResultDf=data.table(topTable(fit_eBayes, 
                        adjust.method="BH", 
                        number=length(fit_eBayes$coefficients), 
                        sort.by = "B"), 
               keep.rownames = T)
  ResultDf$Dif=0
  ResultDf$Dif[ResultDf$adj.P.Val<0.05]=1
  ResultDf$GeneSetNeeded=0
  ResultDf$GeneSetNeeded[ResultDf$rn%in%GeneSetNeeded]=1
  
  print(table(ResultDf[,c('Dif','GeneSetNeeded')]))
  #L=fisher.test(table(ResultDf[,c('Dif','GeneSetNeeded')]), alternative = "less")
  G=fisher.test(table(ResultDf[,c('Dif','GeneSetNeeded')]), alternative = "greater")
  #L_p=L$p.value
  G_p=G$p.value
  #p=c(L_p,G_p)
  p=G_p
  if (class(Plot)=="character") {
    PathNeeded=list()
    PathNeeded[['Path']]=GeneSetNeeded
    RankedDf=ResultDf[order(t), list(rn, t)]
    
    RankedDf$Order=c(1:dim(RankedDf)[1])
    RankedDf$t_z=scale(RankedDf$t)
    
    Ranks=ResultDf[order(t), list(rn, t)]$t
    names(Ranks)=ResultDf[order(t), list(rn, t)]$rn
    GSEA=fgsea(pathways = PathNeeded, 
          stats = Ranks,
          minSize=15,
          maxSize=500,
          nperm=10000)
    print(GSEA)
    P1=plotEnrichment(GeneSetNeeded,
                   Ranks) + xlab("") 
    P2=ggplot(RankedDf,aes(Order,t_z)) +
      geom_bar(stat = "identity",fill='grey') + 
      theme_bw()+ xlab("") +ylab("Ranked list Metric")+ylim(c(-2.5,3))
    
    ggarrange(P1, P2, heights = c(2, 0.7),
              ncol = 1, nrow = 2)
    ggsave(Plot)
  
  }
  return(p)
}
#DiffExpression(Df=DfUse_Z,GeneSetNeeded=FAO_Genes,TargetGene = 'CD47')

