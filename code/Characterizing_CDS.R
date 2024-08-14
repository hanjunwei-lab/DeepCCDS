library(openxlsx)
library(dplyr)
library(tinyarray)
library(GSVA)
library(rcdk)
library(rcellminer)
library(rcellminerData)
library(clusterProfiler)
library(org.Hs.eg.db)
library(tinyarray)
library(PharmacoGx)
library(ggstatsplot)
library(ggplot2)
library(ggpubr)
library(pROC)
library(e1071)
library(glmnet)
################
######GDSC######
################

###Data processing###
#Import model annotations Select model with information on tumor and cell line (model contains SIDM number and corresponding cell line).
GDSC_model<- read.csv("E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Model Annotation\\model_list_20230923.csv",header=T)
GDSC_model <- GDSC_model[,c(1,5,7,10,13)]
GDSC_model <- GDSC_model[which(GDSC_model[,5]=="Cell Line"),]
GDSC_model <- GDSC_model[which(GDSC_model[,4]%in%c("Tumour","Metastasis")),]

##Mutation data##
#Import gene annotation information and retain information with ensembl for subsequent mRNA retention
GDSC_gene_annotated <- read.csv("E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Gene Annotation\\gene_identifiers_20191101.csv",header = T)
GDSC_gene_annotated <- GDSC_gene_annotated[which(GDSC_gene_annotated[,3]!=""),]
GDSC_gene_annotated <- GDSC_gene_annotated[,c(1,3)]

#Importing driver gene mutation information, again only tumor and cell line models were selected to obtain driver gene
GDSC_mutation_info <- read.csv("E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Mutation Data\\mutations_all_20230202.csv",header = T)
GDSC_mutation_info <- GDSC_mutation_info[,c(1,3)]
GDSC_mutation_info <- GDSC_mutation_info[which(GDSC_mutation_info[,2]%in%GDSC_model[,1]),]

#Matching ensembl by SIDG number of the mutation spectrum
GDSC_mutation_info <- inner_join(GDSC_gene_annotated,GDSC_mutation_info,by="gene_id")

GDSC_mutation_matrix <- matrix(0,length(unique(GDSC_mutation_info[,2])),length(unique(GDSC_mutation_info[,3])),dimnames=list(unique(GDSC_mutation_info[,2]),unique(GDSC_mutation_info[,3])))
for(i in 1:nrow(GDSC_mutation_info)){
  GDSC_mutation_matrix[GDSC_mutation_info[i,2],GDSC_mutation_info[i,3]] <- 1
  print(i)
}

#Converting ensembl to symbol 
GDSC_mutation_matrix <- as.data.frame(GDSC_mutation_matrix)
GDSC_mutation_matrix <- trans_exp(GDSC_mutation_matrix,mrna_only = T) 

##Expression data##
#Import tpm matrix
GDSC_exp <- read.csv("E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Expression Data\\rnaseq_tpm_20220624.csv",header=F)
GDSC_exp <- GDSC_exp[-c(2,3,4,5),-2]
GDSC_exp <- GDSC_exp[,c(1,which(GDSC_exp[1,]%in%GDSC_model[,1]))]
colnames(GDSC_exp) <- GDSC_exp[1,]
colnames(GDSC_exp)[1] <- "gene_id"
GDSC_exp <- GDSC_exp[-1,]
#Transformation to ensemble based on gene annotation
GDSC_exp <- inner_join(GDSC_gene_annotated,GDSC_exp,by="gene_id")

write.table(GDSC_exp,file="E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Expression Data\\GDSC_exp.txt",sep="\t",quote=F,row.names=F)
GDSC_exp <- read.table("E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Expression Data\\GDSC_exp.txt",sep="\t",header=T,row.names=1)

GDSC_exp <- na.omit(GDSC_exp)
GDSC_exp <- GDSC_exp[which(rowSums(GDSC_exp[,-1])>0),]
GDSC_exp <- aggregate(.~ensembl_gene_id,mean,data=GDSC_exp)
GDSC_exp <- aggregate(GDSC_exp,by=list(ensembl_gene_id=GDSC_exp$ensembl_gene_id),mean)
rownames(GDSC_exp) <- GDSC_exp[,1]
GDSC_exp <- GDSC_exp[,-1]
GDSC_exp <- log2(GDSC_exp+1)

GDSC_exp <- trans_exp(GDSC_exp,mrna_only = T)
#Cell Lineage Taking Intersection
GDSC_cell_line <- intersect(GDSC_res_cell_line,intersect(colnames(GDSC_exp),colnames(GDSC_mutation_matrix)))
GDSC_cell_line <- sort(GDSC_cell_line)

GDSC_mutation_matrix <- GDSC_mutation_matrix[,GDSC_cell_line]
GDSC_exp <- GDSC_exp[,GDSC_cell_line]
#z-score standardization
GDSC_exp <- scale(t(GDSC_exp))
GDSC_exp <- t(GDSC_exp)
GDSC_exp <- as.data.frame(GDSC_exp)

#######################
######Random Walk######
#######################
#Driver genes
cancer_gene <- read.csv("E:\\Task\\Task1_deep_learning_IC50\\Drug Response\\GDSC\\Data\\Cancer Gene\\Census_allMon.csv")
cancer_gene <- cancer_gene[which(cancer_gene[,5]==1),1]

cancer_mutaition_matrix <- GDSC_mutation_matrix[which(rownames(GDSC_mutation_matrix)%in%cancer_gene),]


cancer_mutaition_gene <- c()
for(i in 1:nrow(cancer_mutaition_matrix)){
  if(as.numeric(table(as.numeric(cancer_mutaition_matrix[i,]))[2])/ncol(cancer_mutaition_matrix) >= 0.01){
    cancer_mutaition_gene <- c(cancer_mutaition_gene,rownames(cancer_mutaition_matrix)[i])
  }
}

#Random Walk Algorithm
rw <- function(W,p0,gamma) {
  
  p0 <- t(p0)
  p0 <- p0/sum(p0)
  PT <- p0
  
  k <- 0
  delta <- 1
  
  Ng <- dim(W)[2]
  for (i in 1:Ng) {
    sumr<-sum(W[i,])
    if(sumr==0)
    {
      W[i,] <- as.numeric(length=length(W[i,]))
    }
    if(sumr>0)
    {
      W[i,] <- W[i,]/sum(W[i,])
    }
  }
  W <- t(W)
  
  while(delta>1e-10) {
    PT1 <- (1-gamma)*W
    PT2 <- PT1 %*% t(PT)
    PT3 <- (gamma*p0)
    PT4 <- t(PT2) + PT3
    delta <- sum(abs(PT4 - PT))
    PT <- PT4
    k <- k + 1
  }
  PT <- t(PT)
  rownames(PT)<-NULL
  return(PT)
}
#Import PPI network
load("C:/Users/Nek/Downloads/adjM_PPI.rdata")
p0 <- c(rep(0,nrow(adjM)))
names(p0) <- rownames(adjM)
gene_start <- cancer_mutaition_gene[which(cancer_mutaition_gene%in%rownames(adjM))]
p0[gene_start] <- 1

r <- 0.9
PT <- rw(adjM,p0,r)
PT <- cbind(rownames(adjM),PT)
PT <- as.data.frame(PT[which(PT[,2]!=0),])
PT[,2] <- as.numeric(PT[,2])
#PT[,2] <- 1/-log2(as.numeric(PT[,2]))
PT <- PT[order(PT[,2],decreasing = T),]
#PT[,2] <- scale(PT[,2])
colnames(PT) <- c("gene","score")

################
######GSEA######
################

FastSEAscore<-function(labels.list,correl.vector = NULL){
  tag.indicator <- labels.list
  no.tag.indicator <- 1 - tag.indicator
  N <- length(labels.list)
  Nh <- length(labels.list[labels.list==1])
  Nm <-  N - Nh
  correl.vector <- abs(correl.vector)
  sum.correl.tag    <- sum(correl.vector[tag.indicator == 1])
  norm.tag    <- 1.0/sum.correl.tag
  norm.no.tag <- 1.0/Nm
  RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
  
  max.ES <- max(RES)
  min.ES <- min(RES)
  ES <- signif(ifelse(max.ES > - min.ES, max.ES, min.ES), digits=5)
  return(ES)
}
spw <- function(df_rank,spwlist){
  r_names<-names(df_rank)
  
  spwlist1 <- lapply(spwlist,function(x ,y){ 
    x1<-x[which(is.element(x,y)==TRUE)]
    return(x1)},r_names)
  
  
  spw_l<-NULL
  qc<-NULL
  for(i in 1:length(spwlist1)){
    if(length(spwlist1[[i]])>1){
      spw_l<-c(spw_l,length(spwlist1[[i]]))
    }else{
      qc<-c(qc,i)
    }
  }
  spwlist1<-spwlist1[-qc]
  
  if(length(spwlist1)>0){
    
    spw_ES<-sapply(spwlist1, function(s){
      tag.indicator <- sign(match(r_names, s, nomatch = 0))
      ES <- FastSEAscore(tag.indicator,correl.vector = df_rank)
      return(ES)
    })
    
    
    rd_geneid<-lapply(c(1:1000),function(r){
      rd_geneid1<-sample(r_names,replace = F)
      return(rd_geneid1)
    })
    rd_geneid <- do.call(cbind, rd_geneid)
    
    
    spw_rd_ES<-lapply(c(1:1000),function(r){
      rdES<-sapply(spwlist1, function(s){
        tag.indicator <- sign(match(rd_geneid[,r], s, nomatch = 0))
        ES <- FastSEAscore(tag.indicator,correl.vector = df_rank)
      })
      return(rdES)
    })
    spw_rd_ES <- do.call(cbind, spw_rd_ES)
    
   
    p.vals <- NULL
    for(i in 1:length(spw_ES)){
      if(is.na(spw_ES[i])==T){
        p.vals[i]<-1
      }else{
        if(spw_ES[i]>=0){
          
          p.vals[i]<-sum(spw_rd_ES[i,] >= spw_ES[i])/1000
        }else{
          
          p.vals[i]<-sum(spw_rd_ES[i,] <= spw_ES[i])/1000
        }
      }
      
    }
    spw_ES[which(is.na(spw_ES)==T)]<-0
    
    drug_data<-data.frame(SubpathwayId=names(spwlist1),Length=spw_l,ES=spw_ES,P_value=p.vals,stringsAsFactors=F)
  }else{
    drug_data<-NULL
  }
  
  return(drug_data)
}
#Import KEGG pathway. Disease pathways have been excluded.
load("D:/2work/listp243.rdata")
genelist <- as.numeric(PT[,2])
names(genelist) <- as.character(PT[,1])
pathway <- spw(genelist,listp243)
fdr <- p.adjust(
  pathway$P_value,  
  method ="fdr"    
)

pathway <- cbind(pathway,fdr)

pathway_0.2 <- pathway[which(pathway[,5]<0.2&pathway[,3]>0),]

geneset <- split(kegg_323_gmt$V2,kegg_323_gmt$V1)
geneset <- geneset[rownames(pathway_0.2)]



######################
######ssGSEA######
######################
GDSC_ssgsea <- gsva(as.matrix(GDSC_exp),geneset,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE,ssgsea.norm=T)
GDSC_ssgsea <- as.data.frame(t(GDSC_ssgsea))
GDSC_ssgsea <- GDSC_ssgsea[sort(rownames(GDSC_ssgsea)),sort(colnames(GDSC_ssgsea))]
GDSC_exp <- GDSC_exp[,sort(colnames(GDSC_exp))]
GDSC_cell_line_map <- as.data.frame(GDSC_cell_line)
colnames(GDSC_cell_line_map) <- "model_id"
GDSC_cell_line_map <- inner_join(GDSC_cell_line_map,GDSC_model,by="model_id")
GDSC_cell_line_map <- GDSC_cell_line_map[order(GDSC_cell_line_map$model_id),]
colnames(GDSC_exp) <- GDSC_cell_line_map[,2]
colnames(GDSC_mutation_matrix) <- GDSC_cell_line_map[,2]
GDSC_cell_line_index <- data.frame(cell_idx=0:(nrow(GDSC_cell_line_map)-1),cell_line=GDSC_cell_line_map[,2])
GDSC_ssgsea <- GDSC_ssgsea[,sort(colnames(GDSC_ssgsea))]
GDSC_ssgsea_input <- cbind(GDSC_cell_line_index,GDSC_ssgsea)
write.csv(GDSC_ssgsea_input,"GDSC_ssgsea_input.csv",row.names = F)
