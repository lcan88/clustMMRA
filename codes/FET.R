options(echo=TRUE) 
args <- commandArgs(trailingOnly = TRUE)
print(args)
dim_sign<-as.numeric(args[1])
th<-as.numeric(args[2])


step2<-as.matrix(read.delim("../results/step2_results.txt",header=TRUE,sep="\t"))
res<-numeric(0)
for(j in 1:dim(step2)[1]){
    mirna_DB<-step2[j,1]
    step2[j,2]<-gsub(" ", "", step2[j,2]) 
    my.data <- read.delim(paste('../results/net_',mirna_DB,'.txt',sep=''), header=FALSE,sep = '\t')
    mir_link<-grep('hsa-',my.data[,2])
    if(length(mir_link)!=0){
      my.data<-my.data[-mir_link,]
    }
    TCGA<-read.delim(paste('../results/boot/',mirna_DB,'.txt',sep=''), row.names=1,header=TRUE,sep = '\t')
    class<-read.delim(paste('../data/signatures/',step2[j,4],step2[j,2],'.txt',sep=''),header=FALSE,sep = '\t')
    net<-as.matrix(unique(my.data[,2]))
    m<-match(row.names(TCGA),net[,1])
    w<-which(is.na(m))
    non_net<-as.matrix(row.names(TCGA[w,]))
    m1<-match(row.names(TCGA),class[,1])
    w1<-which(is.na(m1))
    non_class<-as.matrix(row.names(TCGA[w1,]))
    m2<-match(class[,1],net[,1])
    w2<-which(!is.na(m2))
    class_net<-dim(as.matrix(class[w2,]))[1]
    m3<-match(class[,1],non_net[,1])
    w3<-which(!is.na(m3))
    nonnet_class<-dim(as.matrix(class[w3,]))[1]
    m4<-match(non_class[,1],net[,1])
    w4<-which(!is.na(m4))
    noclass_net<-dim(as.matrix(non_class[w4,]))[1]
    m5<-match(non_class[,1],non_net[,1])
    w5<-which(!is.na(m5))
    nonet_noclass<-dim(as.matrix(non_class[w5,]))[1]
    mat<-rbind(c(class_net,noclass_net),c(nonnet_class,nonet_noclass))
    FET<-fisher.test(mat)
    if(FET$p.value <= th_STR){
        res<-rbind(res,c(mirna_DB,step2[j,2],step2[j,3],step2[j,4],FET$p.value))
    }    
    
}
colnames(res)<-c("mirna","class","mirna_sig","signature_sign","p-value")
write.table(res,'../results/clustMMRA_output.txt',sep='\t',row.names=FALSE,col.names=TRUE)
