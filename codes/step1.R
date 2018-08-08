##### inputs
options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
dataset_file<-args[2]
num_sub<-as.numeric(args[1])
pvalue_th<-as.numeric(args[3])
fc_th<-as.numeric(args[4])
column_seq<-as.matrix(as.numeric(args[5:length(args)]))



library('Matching')
fc_soglia<-fc_th
pval_soglia<-pvalue_th
dim_sign<-num_sub


###computation of fold change and KS on the miRNA matrix
library('Matching')
X<-read.delim(dataset_file,header=TRUE,row.names=1,sep="\t")
X<-as.matrix(X)
cluster<-as.matrix(read.delim("../data/mirbase_cluster_g2.txt",header=FALSE,sep="\t"))
conversion_mi_miamt<-as.matrix(read.delim("../data/conversion_IDpremirna_namemimat.txt",header=TRUE,sep="\t")) 
summary<-numeric(0)  

for(i in 1:dim(cluster)[1]) {
  ###selection miRNA in cluster
  m<-match(cluster[i,],"")
  nnvuoti<-which(is.na(m))
  cluster_i<-as.matrix(cluster[i,nnvuoti])
  list_mirna<-numeric(0)
  for(j in 2:dim(cluster_i)[1]){
    m<-match(conversion_mi_miamt[,2],cluster_i[j,1])
    w<-which(!is.na(m))
    list_mirna<-c(list_mirna,conversion_mi_miamt[w,1])
  }
  list_mirna<-gsub("-5p", "", list_mirna)
  list_mirna<-gsub("-3p", "", list_mirna)
  list_mirna2<-paste(list_mirna,'-1',sep='')
  list_mirna3<-paste(list_mirna,'-2',sep='')
  list_mirna<-as.matrix(rbind(as.matrix(list_mirna),as.matrix(list_mirna2),as.matrix(list_mirna3)))
  m<-match(row.names(X),list_mirna[,1])
  w<-which(!is.na(m))
  X_mir<-as.matrix(X[w,])
  if(dim(X_mir)[1]!=0){
      if(dim(X_mir)[2]==1){
        X_mir<-t(X_mir)
      }
    class<-1
    for(j in seq(1,2*num_sub,2)){
      p1<-numeric(0)
      fc1<-numeric(0)
      for(l in 1:dim(X_mir)[1]){
          N1<-X_mir[i,column_seq[j,1]:column_seq[j+1,1]]
          if(column_seq[j,1]==1){
            S1<-X_mir[i,column_seq[j+2,1]:dim(X)[2]]
          }else if(j+1 != 2*num_sub){
            S1<-rbind(as.matrix(X_mir[i,1:column_seq[j-1,1]]),as.matrix(X_mir[i,column_seq[j+2,1]:dim(X_mir)[2]]))
          }else{
            S1<-X_mir[i,1:column_seq[j-1,1]]
          }
          s1<-ks.boot(N1,S1,alternative="two.sided")$ks.boot.pvalue
          p1<-cbind(p1,s1)
          fc1<-cbind(fc1,mean(N1)-mean(S1))
      }
      p1<-as.matrix(p1)
      fc1<-as.matrix(fc1)
      #####control on significance
      at_least_one_up<-0
      at_least_one_dn<-0
      for(jj in 1:dim(X_mir)[1]){
        if( as.numeric(p1[jj,1]) < pval_soglia && abs(fc1[jj,1]) > fc_soglia){
          if(fc1[jj,1]<0){
            at_least_one_dn<-at_least_one_dn+1
          } 
          if(fc1[jj,1]>=0){
            at_least_one_up<-at_least_one_up+1
          } 
        }
      }
      if(at_least_one_up >=2){
        summary <- rbind(summary,cbind(cluster[i,1],class,'up'))
      }
      if(at_least_one_dn >=2){
        summary <- rbind(summary,cbind(cluster[i,1],class,'down'))
      }
      ##########
      class<-class+1
    }
  }
}    
summary<-as.matrix(summary)
colnames(summary)<-c("miRNA_cluster","subype/class","sign(up/down)")
write.table(summary,file="../results/step1_results.txt",sep="\t",row.names=FALSE,col.names=TRUE)



