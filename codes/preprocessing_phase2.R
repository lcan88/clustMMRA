options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
print(args)
dataset_mirna_file<-args[1]
dataset_mrna_file<-args[2]

library('preprocessCore')
TCGA_mirna<-as.matrix(read.table(dataset_mirna_file,header=TRUE,row.names=1,sep='\t'))
TCGA_mrna<-as.matrix(read.table(dataset_mrna_file,header=TRUE,row.names=1,sep='\t'))
if(all(colnames(TCGA_mirna)==colnames(TCGA_mrna))){
  print("columns of the datasets correctly ordered")
}else{
  print("ATTENTION! the columns of the mRNA and microRNA dataset are not in the same order. Reorder them and rerun") 
}
cluster<-as.matrix(read.delim("../results/step2_results.txt",header=T,sep="\t"))
conversion_mi_miamt<-as.matrix(read.delim("../data/conversion_IDpremirna_namemimat.txt",header=TRUE,sep="\t")) 
conv_clust_singlemir<-as.matrix(read.delim("../data/mirbase_cluster_g2.txt",header=FALSE,sep="\t"))


for(i in 1:dim(cluster)[1]){
  ########### conv cluster mir
  print(i)
  m<-match(conv_clust_singlemir[,1],cluster[i,1])
  w<-which(!is.na(m))
  m<-match(conv_clust_singlemir[w,],"")
  nnvuoti<-which(is.na(m))
  cluster_i<-as.matrix(conv_clust_singlemir[w,nnvuoti])
  list_mirna<-numeric(0)
  for(j in 2:dim(cluster_i)[1]){
    m<-match(conversion_mi_miamt[,2],cluster_i[j,1])
    w<-which(!is.na(m))
    if(length(w)==0){
      print('error')
    }
    list_mirna<-c(list_mirna,conversion_mi_miamt[w,1])
  }
  list_mirna<-gsub("-5p", "", list_mirna)
  list_mirna<-gsub("-3p", "", list_mirna)
  list_mirna2<-paste(list_mirna,'-1',sep='')
  list_mirna3<-paste(list_mirna,'-2',sep='')
  list_mirna<-as.matrix(rbind(as.matrix(list_mirna),as.matrix(list_mirna2),as.matrix(list_mirna3)))
  ##########
  mat<-numeric(0)
  names<-numeric(0)
    for(l in 1:dim(list_mirna)[1]){
      m<-match(row.names(TCGA_mirna),list_mirna[l,1])
      w<-which(!is.na(m))
      if(length(w)!=0){
        mat<-rbind(mat,TCGA_mirna[w,])
        names<-rbind(names,row.names(TCGA_mirna)[w])
      }
    }
  mat<-as.matrix(mat)
  dim_cluster<-dim(mat)[1]
  names<-as.matrix(names)
	matrice<-rbind(TCGA_mrna,mat)
	row.names(matrice)[(dim(matrice)[1]-dim(names)[1]+1):dim(matrice)[1]]<-names

	##################### matrix normalization
	matrice<-exp(matrice)
	matrice2<-normalize.quantiles(matrice)
	rownames(matrice2)<-rownames(matrice)
	colnames(matrice2)<-colnames(matrice)
	matrice2<-log2(matrice2)
	standev<-numeric(0)
	for(j in 1:(dim(matrice2)[1]-dim_cluster)){
	  standev<-rbind(standev,sd(matrice2[j,]))
	}
	w<-which(standev >1.2)
	matrice_new<-rbind(matrice2[w,],matrice2[(dim(matrice2)[1]-dim_cluster):dim(matrice2)[1],])
	x<-cbind(rownames(matrice_new),rownames(matrice_new),matrice_new)
	vet<-cbind('genes','desc',t(as.matrix(colnames(matrice_new))))
	mat<-rbind(vet,x)
	write.table(mat,paste('../results/boot/',mirna[i],'.txt',sep=''),sep='\t',col.names=FALSE,row.names=FALSE)
}