#' @title Create Univariable Cox PH Regression Table 
#'
#' @description The function creates a result table for univariable Cox PH regression. \code{veteran} dataset in the \code{survival} package is employed as an example data.
#' @param time a variable for time to event
#' @param status a variable for event (or censoring) indicator (1 indicates event, 0 is censored)
#' @param covariate a vector of variables included in the univariable analysis
#' @param label a vector of labels shown in the table
#' @param type a vector of types either categorical or numerical (\code{"cat"} or \code{"num"})
#' @param reflevel a vector of reference levels for categorical variables. NA is used for numerical variable.
#' @param data data frame dataset
#' @export
#' @import dplyr 
#' @import survival
#' @return matrix with estimates of coefficient, HR, and P-value
#' @examples 
#' library(survival)
#' data(veteran)
#' df <- veteran
#' uni.cox.reg(time="time", status="status", 
#'             covariate=c("age", "trt","celltype"), 
#'             label=c("Age", "Treatment", "Celltype"), 
#'             type=c("num","cat","cat"), 
#'             reflevel=c(NA, "1","squamous"), 
#'             data=df)     
#'                         
#' 


uni.cox.reg <- function(time,status,covariate,label,type,reflevel=NULL,data){
  
  data <- easytab::na.handle.2(data, variable=covariate, type, useNA=FALSE)  

  covariate <- as.character(covariate)
  label     <- as.character(label)
  
  numvar=length(covariate)
  
  for (i in 1:numvar){
    if (type[i]=="cat"){  
      #Assign reference levels for categorical variables
      data[,which(colnames(data)==covariate[i])]<-as.factor( data[,which(colnames(data)==covariate[i])])
      data[,which(colnames(data)==covariate[i])]<- relevel(data[,which(colnames(data)==covariate[i])], ref = as.character(reflevel[i])) 
    }
  }
  
  #Create Survival formula
  a=paste("Surv(",time,sep='')
  b=paste(status,")~",sep='')
  func=paste(paste(a,",",sep=""),b,sep="")
  
  univ_formulas <- sapply(covariate,
                          function(x) as.formula(paste(func, x))) #create univariate formulas
  
  univ_models <- lapply( univ_formulas, function(x){coxph(x, data = data)})
  # Extract data 
  allres<-NULL
  for (k in 1:length(univ_models)){
    x=univ_models[[k]]
    
    x <- summary(x)
    p.value<- round(x$coef[,5], 5)
    beta<-round(x$coef[,1], 2);#coeficient beta
    HR <-round(x$coef[,2], 2); #exp(beta)
    HR.confint.lower <- round(x$conf.int[,"lower .95"], 2)
    HR.confint.upper <- round(x$conf.int[,"upper .95"],2)
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-cbind(beta, HR, p.value)
    colnames(res)<-c("coefficient", "HR (95% CI for HR)",  
                     "P-value")
    allres<-c(allres,list(res))
  }
  
  Table=do.call("rbind", allres)
  
  ET <- list()
  for (i in 1:numvar){
    df <- data %>% dplyr::filter(!is.na(eval(parse(text=covariate[i])))) %>%
      dplyr::filter(eval(parse(text=covariate[i]))!="")
    
    comm <- paste("(Surv(",time,",",status," == 1) ~ 1", ")", sep="")
    fit  <- survfit(eval(parse(text=comm)), data=df) 
    tmp.e <- summary(fit)$table["events"]
    tmp.n <- summary(fit)$table["records"]
    ET[i] <- paste(tmp.e,"/",tmp.n,sep="") 
  }
  
  #Create labels for the first/second columns of the table
  Col1<-NULL
  for (i in 1:numvar){
    if (type[i]=="cat"){  
      catvar <- data[,which(colnames(data)==covariate[i])]
      catvar <- catvar[!is.na(as.character(catvar)) & as.character(catvar)!=""]   
      catvar<-factor(catvar)
      new.level=levels(catvar)[levels(catvar)!=reflevel[i]]
      numCat<-length(new.level)  
      
      cat.tab=NULL
      for(j in 1:numCat){
        cat.tab<-rbind(cat.tab,c(NA,paste(new.level[j]," vs. ",reflevel[i],sep=""),NA))
      }
      cat.tab[1,1] <- c(label[i])
      cat.tab[1,3] <- ET[[i]]    
      Col1<-rbind(Col1,cat.tab)
    }
    if (type[i]=="num"){  
      num.var<-c(label[i],NA, ET[[i]]) 
      Col1<-rbind(Col1,num.var)
    }
    colnames(Col1) <- c("","","Event/Total")
  }
  Table<-cbind(Col1,Table)
  
  rownames(Table) <- NULL
  Table[is.na(Table)]<-rep('',sum(is.na(Table)))
  
  # P values NEJM format
  tmp <- tmpc <- as.numeric(Table[,'P-value'])
  tmp[tmpc >= 0.005 & !is.na(tmpc)] <- formatC(tmpc[tmpc >= 0.005 & !is.na(tmpc)], digit=2, format="f")
  tmp[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)] <- formatC(tmpc[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)], digit=3, format="f")
  tmp[tmpc < 0.001] <- "<.001"
  tmp[tmpc < 0.05 & !is.na(tmpc)]  <- as.character(paste("**",tmp[tmpc < 0.05 & !is.na(tmpc)],"**",sep=""))  
  Table[,'P-value'] <- tmp
  
  # return Table
  return(Table)
}




#' @title Create Multivariable Cox PH Regression Table 
#'
#' @description The function creates a result table for multivariable Cox PH regression. \code{veteran} dataset in the \code{survival} package is employed as an example data.
#' @param time a variable for time to event
#' @param status a variable for event (or censoring) indicator (1 indicates event, 0 is censored)
#' @param covariate a vector of variables included in the multivariable analysis
#' @param label a vector of labels shown in the table
#' @param type a vector of types either categorical or numerical (\code{"cat"} or \code{"num"})
#' @param reflevel a vector of reference levels for categorical variables. NA is used for numerical variable.
#' @param data data frame dataset
#' @export
#' @import dplyr 
#' @import survival
#' @return matrix with estimates of coefficient, HR, and P-value
#' @examples 
#' library(survival)
#' data(veteran)
#' df <- veteran
#' mul.cox.reg(time="time", status="status", 
#'             covariate=c("age", "trt","celltype"), 
#'             label=c("Age", "Treatment", "Celltype"), 
#'             type=c("num","cat","cat"), 
#'             reflevel=c(NA, "1","squamous"), 
#'             data=df)
#'                             



mul.cox.reg <- function(time,status,covariate,label,type,reflevel=NULL, data){
  
  data <- na.handle.2(data, variable=covariate, type, useNA=FALSE)  

  covariate <- as.character(covariate)
  label     <- as.character(label)
  
  numvar=length(covariate)
  
  for (i in 1:numvar){
    if (type[i]=="cat"){  
      #Assign reference levels for categorical variables
      data[,which(colnames(data)==covariate[i])]<-as.factor( data[,which(colnames(data)==covariate[i])])
      data[,which(colnames(data)==covariate[i])]<- relevel(data[,which(colnames(data)==covariate[i])], ref = as.character(reflevel[i])) 
    }
  }
  
  
  #Create Survival formula
  a=paste("Surv(",time,sep='')
  b=paste(status,")",sep='')
  func=paste(paste(a,",",sep=""),b,sep="")
  
  #fmla=as.formula(paste(func, covariate))
  
  fmla <- as.formula(paste(paste(func,"~ "), paste(covariate, collapse= "+")))
  
  model <-  coxph(fmla, data = data)
  
  x <- summary(model)
  p.value<-round(x$coef[,5], digits=4)
  beta<-round(x$coef[,1], digits=2);#coeficient beta
  HR <-round(x$coef[,2], digits=2);#exp(beta)
  HR.confint.lower <- round(x$conf.int[,"lower .95"], 2)
  HR.confint.upper <- round(x$conf.int[,"upper .95"],2)
  HR <- paste0(HR, " (", 
               HR.confint.lower, "-", HR.confint.upper, ")")
  res<-as.matrix(cbind(beta, HR, p.value))
  colnames(res)<-c("coefficient", "HR (95% CI for HR)",  
                   "P-value")
  
  #Create labels for the first/second columns of the table
  Col1<-NULL
  for (i in 1:numvar){
    if (type[i]=="cat"){  
      catvar <- data[,which(colnames(data)==covariate[i])]
      catvar<-factor(catvar)
      new.level=levels(catvar)[levels(catvar)!=reflevel[i]]
      numCat<-length(new.level)
      
      cat.tab=NULL
      for(j in 1:numCat){
        cat.tab<-rbind(cat.tab,c(NA,paste(new.level[j]," vs. ",reflevel[i],sep="")))
      }
      cat.tab[1,1]<-c(label[i])
      Col1<-rbind(Col1,cat.tab)
    }
    if (type[i]=="num"){  
      num.var<-c(label[i],NA)
      Col1<-rbind(Col1,num.var)
    }
  }
  res<-cbind(Col1,res)
  
  rownames(res) <- NULL
  res[is.na(res)]<-rep('',sum(is.na(res)))
  
  # P values NEJM format
  tmp <- tmpc <- as.numeric(res[,'P-value'])
  tmp[tmpc > 0.005 & !is.na(tmpc)] <- formatC(tmpc[tmpc > 0.005 & !is.na(tmpc)], digit=2, format="f")
  tmp[tmpc <= 0.005 & tmpc >= 0.001 & !is.na(tmpc)] <- formatC(tmpc[tmpc <= 0.005 & tmpc >= 0.001 & !is.na(tmpc)], digit=3, format="f")
  tmp[tmpc < 0.001] <- "<.001"
  tmp[tmpc < 0.05 & !is.na(tmpc)]  <- as.character(paste("**",tmp[tmpc < 0.05 & !is.na(tmpc)],"**",sep=""))  ## added 01/04/2019
  res[,'P-value'] <- tmp
  
  return(res)
}

