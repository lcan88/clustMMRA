ttest<-read.delim("./results/step1_results.txt",header=T,sep="\t")
ttest<-as.matrix(ttest)
cluster<-as.matrix(read.delim("../data/mirbase_cluster_g2.txt",header=FALSE,sep="\t"))
conversion_mi_miamt<-as.matrix(read.delim("../data/conversion_IDpremirna_namemimat.txt",header=TRUE,sep="\t")) 
X<-read.delim("../data/all_no_mTB.txt",header=FALSE,sep="\t")
X<-as.matrix(X)
Y<-read.delim("../data/sperimentale.txt",header=FALSE,sep="\t")
Y<-as.matrix(Y)
genes<-as.matrix(unique(X[,2]))
riassunto<-numeric(0)
####### I classify miRNA-target interaction according to number of DBs
DB_2<-numeric(0)
DB_3<-numeric(0)
DB_4<-numeric(0)
for(supp in 4:2){
    m<-match(X[,3],supp)
    w<-which(!is.na(m))
    if(supp==4){
        DB_4<-X[w,]
        DB_3<-rbind(DB_3,X[w,])
        DB_2<-rbind(DB_2,X[w,])
    }
    if(supp==3){
        DB_3<-rbind(DB_3,X[w,])
        DB_2<-rbind(DB_2,X[w,])
    }
    if(supp==2){
        DB_2<-rbind(DB_2,X[w,])
    }
}
###############################################################################
########## target enrichment test
#### count microR's targets
for(j in 1:dim(ttest)[1]){
  print(j)
  class<-as.numeric(ttest[j,2])
  for(segno in c('up','down')){
    sign<-as.matrix(read.delim(paste("../data/signatures/",segno,class,".txt",sep=""),header=FALSE,sep="\t"))
    cl<-ttest[j,1]
    m<-match(cluster[,1],cl)
    w<-which(!is.na(m))
    if(length(w)==0){
      print('error cluster name not found')
    }
    m<-match(cluster[w,],"")
    nnvuoti<-which(is.na(m))
    cluster_i<-as.matrix(cluster[w,nnvuoti])
    list_mirna<-numeric(0)
    for(ss in 2:dim(cluster_i)[1]){
      m<-match(conversion_mi_miamt[,2],cluster_i[ss,1])
      w<-which(!is.na(m))
      list_mirna<-c(list_mirna,conversion_mi_miamt[w,1])
    }
    list_mirna<-as.matrix(list_mirna) #####contiens mirna in the cluster 
    ###DB restricted to the signature
    #DB2 in sign
    m2<-match(DB_2[,2],sign[,1])
    w2<-which(!is.na(m2))
    DB_2_sign<-as.matrix(DB_2[w2,])
    if(dim(DB_2_sign)[2]==1){
      DB_2_sign<-t(DB_2_sign)
    }
    #DB3 in sign
    m3<-match(DB_3[,2],sign[,1])
    w3<-which(!is.na(m3))
    DB_3_sign<-as.matrix(DB_3[w3,])
    if(dim(DB_3_sign)[2]==1){
      DB_3_sign<-t(DB_3_sign)
    }
    #DB4 in sign
    m4<-match(DB_4[,2],sign[,1])
    w4<-which(!is.na(m4))
    DB_4_sign<-as.matrix(DB_4[w4,])
    if(dim(DB_4_sign)[2]==1){
      DB_4_sign<-t(DB_4_sign)
    }
    #Db sper sign
    ms<-match(Y[,2],sign[,1])
    ws<-which(!is.na(ms))
    Y_sign<-as.matrix(Y[ws,])
    if(dim(Y_sign)[2]==1){
      Y_sign<-t(Y_sign)
    }
    #### compute targets of cluster in the signature
    t_2<-numeric(0)
    t_3<-numeric(0)
    t_4<-numeric(0)
    t_s<-numeric(0)
    for(mir in 1:dim(list_mirna)[1]){
      ####target in sign 2db
      m2<-match(DB_2_sign[,1],list_mirna[mir,1])
      w2<-which(!is.na(m2))
      if(length(w2)!=0){
        t_2<-union(t_2,as.matrix(DB_2_sign[w2,2]))
      }
      ####target in sign 3db
      m3<-match(DB_3_sign[,1],list_mirna[mir,1])
      w3<-which(!is.na(m3))
      if(length(w3)!=0){
        t_3<-union(t_3,as.matrix(DB_3_sign[w3,2]))
      }
      ####target in sign 4db
      m4<-match(DB_4_sign[,1],list_mirna[mir,1])
      w4<-which(!is.na(m4))
      if(length(w4)!=0){
        t_4<-union(t_4,as.matrix(DB_4_sign[w4,2]))
      }
      ####target in sign sperimentally derived
      ms<-match(Y_sign[,1],list_mirna[mir,1])
      ws<-which(!is.na(ms))
      if(length(ws)!=0){
        t_s<-union(t_s,as.matrix(Y_sign[ws,2]))
      }
    }
    t_num_2<-length(t_2)
    t_num_3<-length(t_3)
    t_num_4<-length(t_4)
    t_num_s<-length(t_s)
###########################null model for random genes of same dimension of targets of each cluster    
    target_2r<-numeric(0)
    target_3r<-numeric(0)
    target_4r<-numeric(0)
    target_sr<-numeric(0)
    for(run in 1:1000){
        print(run)
      s<-as.matrix(sample(1:dim(genes)[1], dim(sign)[1], replace = FALSE, prob = NULL))
      sign_r<-as.matrix(genes[s,1])
      ###DB ristrected to sign_r
      #DB2 in sign
      m2<-match(DB_2[,2],sign_r[,1])
      w2<-which(!is.na(m2))
      DB_2_sign_r<-as.matrix(DB_2[w2,])
      if(dim(DB_2_sign_r)[2]==1){
        DB_2_sign_r<-t(DB_2_sign_r)
      }
      #DB3 in sign
      m3<-match(DB_3[,2],sign_r[,1])
      w3<-which(!is.na(m3))
      DB_3_sign_r<-as.matrix(DB_3[w3,])
      if(dim(DB_3_sign_r)[2]==1){
        DB_3_sign_r<-t(DB_3_sign_r)
      }
      #DB4 in sign
      m4<-match(DB_4[,2],sign_r[,1])
      w4<-which(!is.na(m4))
      DB_4_sign_r<-as.matrix(DB_4[w4,])
      if(dim(DB_4_sign_r)[2]==1){
        DB_4_sign_r<-t(DB_4_sign_r)
      }
      #Db sper sign
      ms<-match(Y[,2],sign_r[,1])
      ws<-which(!is.na(ms))
      Y_sign_r<-as.matrix(Y[ws,])
      if(dim(Y_sign_r)[2]==1){
        Y_sign_r<-t(Y_sign_r)
      }
      #### compute targets of the cluster in sign_r
      t_2_r<-numeric(0)
      t_3_r<-numeric(0)
      t_4_r<-numeric(0)
      t_s_r<-numeric(0)
      for(mir in 1:dim(list_mirna)[1]){
        ####target in sign 2db
        m2<-match(DB_2_sign_r[,1],list_mirna[mir,1])
        w2<-which(!is.na(m2))
        if(length(w2)!=0){
          t_2_r<-union(t_2_r,as.matrix(DB_2_sign_r[w2,2]))
        }
        ####target in sign 3db
        m3<-match(DB_3_sign_r[,1],list_mirna[mir,1])
        w3<-which(!is.na(m3))
        if(length(w3)!=0){
          t_3_r<-union(t_3_r,as.matrix(DB_3_sign_r[w3,2]))
        }
        ####target in sign 4db
        m4<-match(DB_4_sign_r[,1],list_mirna[mir,1])
        w4<-which(!is.na(m4))
        if(length(w4)!=0){
          t_4_r<-union(t_4_r,as.matrix(DB_4_sign_r[w4,2]))
        }
        ####target in sign sperimentally oderived interactions
        ms<-match(Y_sign_r[,1],list_mirna[mir,1])
        ws<-which(!is.na(ms))
        if(length(ws)!=0){
          t_s_r<-union(t_s_r,as.matrix(Y_sign_r[ws,2]))
        }
      }
      t_num_2_r<-length(t_2_r)
      t_num_3_r<-length(t_3_r)
      t_num_4_r<-length(t_4_r)
      t_num_s_r<-length(t_s_r)
      target_2r<-c(target_2r,t_num_2_r)
      target_3r<-c(target_3r,t_num_3_r)
      target_4r<-c(target_4r,t_num_4_r)
      target_sr<-c(target_sr,t_num_s_r)
    }
    soglia2<-quantile(target_2r,.95)
    soglia3<-quantile(target_3r,.95)
    soglia4<-quantile(target_4r,.95)
    soglias<-quantile(target_sr,.95)
    if(t_num_2>soglia2){
      riassunto<-rbind(riassunto,cbind(cl,class,ttest[j,3],segno,2))
    }
    if(t_num_3>soglia3){
      riassunto<-rbind(riassunto,cbind(cl,class,ttest[j,3],segno,3))
    }
    if(t_num_4>soglia4){
      riassunto<-rbind(riassunto,cbind(cl,class,ttest[j,3],segno,4))
    }
    if(t_num_s>soglias){
      riassunto<-rbind(riassunto,cbind(cl,class,ttest[j,3],segno,'s'))
    }
  }
}

colnames(riassunto)<-c("cluster","Expression_class","miRNAcluster_sign","signature","minimal_number_DB")

write.table(riassunto,"../results/step2_results.txt",col.names=TRUE,row.names=FALSE,sep="\t")





