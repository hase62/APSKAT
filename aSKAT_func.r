#install.packages("SKAT")
library(SKAT)

#Efficient Reading
read.scan<-function(x, cut){
  junk<-scan(x, "character", sep="\n")
  p <- length(junk)
  ID <- rep("---", p)
  for(i in 1:p){
    junk2<-unlist(strsplit(junk[i], cut))
    score<-junk2[-c(1)]
    ID[i]<-junk2[1]
    if(i==1)book_n<-matrix(0, p, length(score))
    book_n[i,]<-score
  }
  rownames(book_n)<-ID
  colnames(book_n)<-scan(x, "character", sep=cut, nline=1)[-1]
  return(book_n)
}

#aSKAT Using bfile (PLINK Format)
aSKAT_bfile<-function(prefix, SNPSetID, significance_level = 0.05 / 20000, ID="job", 
                      n.maxres = 5.0 / (0.05 / 20000), conf_1=significance_level, conf_2=conf_1, 
                      covariate = NA){
  ###Generate SSD###
  bfile<-paste(prefix, c("bed", "bim", "fam"), sep=".")
  SSD<-paste(ID, ".SSD", sep="")
  INFO<-paste(ID, ".Info", sep="")
  Generate_SSD_SetID(bfile[1], bfile[2], bfile[3], SNPSetID, SSD, INFO)
  
  ###Preparation of Generating obj###
  phenotype<-Read_Plink_FAM(bfile[3])[,6]
  flag<-FALSE
  if(!is.na(covariate)){
     flag<-TRUE
     fam<-Read_Plink_FAM(bfile[3])
     covariate<-t(sapply(scan(covariate,"character",sep="\n"), function(x) strsplit(x, " ")[[1]]))
     covariate<-covariate[match(fam[,1],covariate[,1]),][,-c(1)]
     covariate<-apply(covariate, 2, function(x) as.numeric(x))
  }
  print(phenotype)
  sInfo<-Open_SSD(SSD, INFO)
  setIndex<-sInfo$SetInfo[,1]
    
  ###Continue / Binary###
  if(length(unique(phenotype))>2){
    type<-"C"
  }else{
    type<-"D"
  }
  
  ###Main###
  if(!flag){
    obj<-SKAT_Null_Model(phenotype ~ 1, out_type=type, n.Resampling=0, 
      type.Resampling="bootstrap", Adjustment=TRUE)
  }else{
    obj<-SKAT_Null_Model(phenotype ~ covariate, out_type=type, n.Resampling=0, 
      type.Resampling="bootstrap", Adjustment=TRUE)
  }
  result<-SKAT.SSD.All(sInfo, obj)
  to_be_evaluated<-ifelse(result$results$P.value < 0.5, 0, 1)
  
  ###Adaptive Procedure###
  n.res=n.maxres / 10
  if(n.res > 1.0e4) n.res <- 1.0e4
  current_count<-0
  pvalues<-sapply(1:length(to_be_evaluated), function(x) as.list(1))
  while(length(which(to_be_evaluated==0)) > 0 & current_count <= n.maxres){
    if(!flag){
      obj<-SKAT_Null_Model(phenotype ~ 1, out_type=type, n.Resampling=n.res, 
        type.Resampling="bootstrap", Adjustment=TRUE)
    }else{
      obj<-SKAT_Null_Model(phenotype ~ covariate, out_type=type, n.Resampling=n.res, 
        type.Resampling="bootstrap", Adjustment=TRUE)
    }
    th_2<-which(dnbinom(0, seq(1:100), .05) < conf_2)[1]
    for(i in 1:length(to_be_evaluated)){
      if(to_be_evaluated[i]==0){
        temp<-SKAT.SSD.OneSet_SetIndex(sInfo, setIndex[i], obj)
  	    pvalues[[i]]<-c(pvalues[[i]], temp$p.value.resampling)
  	    if(current_count==0) pvalues[[i]]<-pvalues[[i]][-c(1)]
  	    result$results[i,2]<-Get_Resampling_Pvalue_1(temp$p.value, pvalues[[i]])$p.value
        if(temp$p.value==1) result$results[i,2]<-temp$p.value
        
  	    #Stop Criterion: 1
  	    if(length(which(temp$p.value > pvalues[[i]][1:th_2])) == th_2){
  	      to_be_evaluated[i]<-1
  	      #Stop Criterion: 2
  	    }else if(length(which(temp$p.value < pvalues[[i]])) == length(pvalues[[i]])){
  	      if(dnbinom(0, length(pvalues[[i]]), (1 - significance_level)) < conf_2){
  	        to_be_evaluated[i]<-1
  	        result$results[i,2]<-1/(length(pvalues[[i]]) + 1)
  	      }
  	      #Stop Criterion: 3
  	    }else{
  	      dif<-qnorm(p = conf_1, mean = 0, 
  	                 sd = sqrt(result$results[i,2] * (1 - result$results[i,2]) / length(pvalues[[i]])))
  	      if(result$results[i,2] + dif > significance_level | result$results[i,2] - dif < significance_level){
  	        to_be_evaluated[i]<-1
  	      }
  	    }
  	    current_count<-length(pvalues[[i]])
      }
    }
    gc();gc();
  }
  
  result<-result$results
  rownames(result)<-sInfo$SetInfo[,2]
  Close_SSD()
  warnings()
  return(result)
}

#aSKAT Using Matrix
aSKAT_matrix<-function(geno_matrix, phenotype, significance_level, ID, n.maxres, 
                       conf_1 = significance_level, conf_2=conf_1, covariate=NA, is_dose=FALSE){
  ###Continue / Binary###
  if(length(unique(phenotype)) > 2){
    type<-"C"
  }else{
    type<-"D"
  }
  
  ###Main###
  if(is.na(covariate)){
    obj<-SKAT_Null_Model(phenotype ~ 1, out_type=type, n.Resampling=0, 
                         type.Resampling="bootstrap", Adjustment=TRUE)
  }else{
    obj<-SKAT_Null_Model(phenotype ~ covariate, out_type=type, n.Resampling=0, 
                         type.Resampling="bootstrap", Adjustment=TRUE)
  }
  if(type=="D"){
    out<-SKATBinary(geno_matrix, obj, is_dosage = is_dose)
  }else{
    out<-SKAT(geno_matrix, obj, is_dosage = is_dose)
  }
  
  ###Adaptive Procedure###
  n.res=n.maxres / 10
  if(n.res > 1.0e4) n.res <- 1.0e4
  current_count<-0
  th_2<-which(dnbinom(0, seq(1:100), .05) < conf_2)[1]
  pvalues<-NULL
  while(current_count <= n.maxres){
    if(is.na(covariate)){
      obj<-SKAT_Null_Model(phenotype ~ 1, out_type=type, n.Resampling=n.res, 
                           type.Resampling="bootstrap", Adjustment=TRUE)
    }else{
      obj<-SKAT_Null_Model(phenotype ~ cov, out_type=type, n.Resampling=n.res, 
                           type.Resampling="bootstrap", Adjustment=TRUE)
    }
    if(type=="D"){
      perm<-SKATBinary(geno_matrix, obj, is_dosage = is_dose)
    }else{
      perm<-SKAT(geno_matrix, obj, is_dosage = is_dose)
    }
    pvalues<-c(pvalues, perm$p.value.resampling)
    result<-Get_Resampling_Pvalue_1(out$p.value, pvalues)$p.value

    #Stop Criterion: 1
    if(length(which(out$p.value > pvalues[1:th_2])) == th_2){
      break
      #Stop Criterion: 2
    }else if(length(which(out$p.value < pvalues)) == length(pvalues)){
      if(dnbinom(0, length(pvalues), (1 - significance_level)) < conf_2){
        break
      }
      #Stop Criterion: 3
    }else{
      dif<-qnorm(p = conf_1, mean = 0, 
                 sd = sqrt(result * (1 - result) / length(pvalues)))
      if(result + dif > significance_level | result - dif < significance_level){
        break
      }
    }
    current_count<-length(pvalues)
  }
  out$p.value<-result
  return(out)
}
