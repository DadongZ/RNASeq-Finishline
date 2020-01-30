rm(list=ls())
library(tidyverse)
library(readxl)
library(shiny)
library(shinyWidgets)
library(ggplot2)
library(DESeq2)
library(reshape2)
library(calibrate)
library(DEFormats)
library(gage)
library(HGNChelper)
library(org.Hs.eg.db)
data("hgnc.table")
options(stringsAsFactors=FALSE)

##Local func
getdgeres<-function(count, manifest, comparison, level1, level0){
  gcol<-colnames(count)[1]
  count<-count%>%data.frame%>%
    column_to_rownames(var=gcol)
  pcol<-colnames(manifest)[1]
  manifest<-manifest%>%data.frame%>%
    column_to_rownames(var=pcol)
  
  if (!identical(colnames(count), rownames(manifest))) {
    message("Sample names in count matrix must be identical to the sample names in pheno data")
    break
  } 
  
  #dds object
  dds <- DESeqDataSetFromMatrix(countData=count, 
                                colData=manifest, 
                                design=as.formula(paste0("~", comparison)))
  dds<-dds[rowSums(counts(dds))>0, ]
  
  normalized_dat<-assay(rlog(dds, blind=TRUE))%>%data.frame%>%rownames_to_column(var="Gene")
  dds<-DESeq(dds)
  #results
  res<-results(dds, contrast=c(comparison, level1, level0))%>%data.frame%>%
    rownames_to_column(var="Gene")%>%
    filter(complete.cases(.))%>%
    mutate(log10padj=-log10(padj))%>%
    arrange(pvalue)
  return(list(results=res, 
              normal=normalized_dat))
}

getvolcano<-function(dgeres, AdjustedCutoff=0.05, FCCutoff=1){
  brushdat<-dgeres%>%
    mutate(Significance=factor(ifelse(log2FoldChange < (-FCCutoff) & padj < AdjustedCutoff, "FC_FDR_Down",
                               ifelse(log2FoldChange > FCCutoff & padj < AdjustedCutoff, "FC_FDR_Up", "NS"))))
  
  #plot
  p <-ggplot(brushdat, aes(x=log2FoldChange, y=log10padj)) +
    geom_point(aes(color=Significance), alpha=1/2, size=2) +
    theme(legend.position = "bottom", legend.title = element_blank())+
    theme_bw(base_size=16) +
    xlab(bquote(~Log[2]~ "fold change")) +
    ylab(bquote(~-Log[10]~adjusted~italic(P))) +
    geom_vline(xintercept=c(-FCCutoff, FCCutoff), linetype="longdash", colour="black", size=0.4) +
    geom_hline(yintercept=-log10(AdjustedCutoff), linetype="longdash", colour="black", size=0.4)
  
  return(list(plot=p,
              brush=brushdat))
}

getgores<-function(dgeres, species="Human"){
  GO<-go.gsets(species = species)
  brushdat<-dgeres%>%
    mutate(entrez=mapIds(org.Hs.eg.db, dgeres$Gene, 'ENTREZID', 'SYMBOL', multiVals = "first"))
  fcs<-brushdat$log2FoldChange
  names(fcs)<-brushdat$entrez
  bpres<-gage(fcs, gsets=GO$go.sets[GO$go.subs$BP], ref=NULL, samp = NULL)
  ccres<-gage(fcs, gsets=GO$go.sets[GO$go.subs$CC], ref=NULL, samp = NULL)
  mfres<-gage(fcs, gsets=GO$go.sets[GO$go.subs$MF], ref=NULL, samp = NULL)
  return(list(bpupper=bpres[['greater']],
              bpless=bpres[['less']],
              ccupper=ccres[['greater']],
              ccless=ccres[['less']],
              mfupper=mfres[['greater']],
              mfless=mfres[['less']]))
}

