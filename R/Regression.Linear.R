#' @title Create Univariable Linear Regression Table 
#'
#' @description The function creates a result table of univariable linear regression.
#' @param outcome a continous outcome variable
#' @param covariate a vector of variables included in the univariable linear analysis
#' @param label a vector of labels shown in the table for covariates
#' @param type a vector of types either categorical or numerical ("\code{cat}" or "\code{num}")
#' @param reflevel a vector of reference levels for categorical variables. \code{NA} is used for continuous variable.
#' @param data data frame dataset
#' @export
#' @importFrom stats glm
#' @return matrix with estimates of coefficients and P-value
#' @examples 
#' # Use mtcars dataset for example
#' df <- mtcars
#' df$mpg.cat <- ifelse(df$mpg > 18, 1, 0)
#' 
#' # Create variable, label, type, and reference information
#' r1 <- c("hp", "Horse power", "num", NA)
#' r2 <- c("wt", "Weight", "num", NA)
#' r3 <- c("gear", "Gear", "cat", "3")
#' input  <- rb(3)
#' 
#' uni.linear.reg(outcome="mpg", 
#'                covariate=input[,"variable"],
#'                label=input[,"label"],
#'                type=input[,"type"],
#'                reflevel=input[,"ref"], 
#'                data=df
#'                )  
#'                              



uni.linear.reg <- function(outcome,covariate,label,type,reflevel=NULL,data){
  
  data <- na.handle.2(data, variable=covariate, type, useNA=FALSE)  
  
  digit=2
  label <- as.character(label)
  
  numvar=length(covariate)

  for (i in 1:numvar){
    if (type[i]=="cat"){  
      #Assign reference levels for categorical variables
      data[,which(colnames(data)==covariate[i])]<-as.factor( data[,which(colnames(data)==covariate[i])])
      data[,which(colnames(data)==covariate[i])]<- relevel(data[,which(colnames(data)==covariate[i])], ref = as.character(reflevel[i])) # HH: Added as.character
    }
  }
  
  #Create outcome formula
  fmla <- paste(outcome,"~ ")
  
  univ_formulas <- sapply(covariate,
                          function(x) as.formula(paste(fmla, x))) #create univariate formulas
  
  univ_models <- lapply( univ_formulas, function(x){glm(x, family=gaussian(link="identity"), data=data, na.action = na.omit)}) 
  #univ_models <- lapply( univ_formulas, function(x){glm(x, family=binomial(link='logit'), data=na.omit(data))})
  # HH: Fixed data=na.omit(data) to avoid whole na data frame.
  
  # Extract data 
  allres<-NULL
  for (k in 1:length(univ_models)){
    x=univ_models[[k]]
    
    coef<-round(x$coefficients[-1],digit)
    OR<-round(exp(x$coefficients[-1]),digit)
    se<-round(coef(summary(x))[,2],digit)[-1]
    pvalues<-round(coef(summary(x))[,4],4)[-1]
    
    #Confidence interval for the coefficients
    CI.LL.beta<-x$coefficients[-1]-1.96*se
    CI.UL.beta<-x$coefficients[-1]+1.96*se
    interval.beta<-paste('(',round(CI.LL.beta,digit),';',round(CI.UL.beta,digit),')',sep='')
    
    
    res<-cbind(coef, interval.beta, pvalues)
    colnames(res)<-c("Coefficient", "95% CI ", "P-value")
    
    allres<-c(allres,list(res))
  }
  
  Table=do.call("rbind", allres)
  
  #Create labels for the first/second columns of the table
  Col1<-NULL
  for (i in 1:numvar){
    if (type[i]=="cat"){  
      catvar <- data[,which(colnames(data)==covariate[i])]
      catvar <- catvar[!is.na(as.character(catvar)) & as.character(catvar)!=""]  #HH: added additional condition to avoid "". 
      catvar<-factor(catvar)
      new.level=levels(catvar)[levels(catvar)!=reflevel[i]]
      numCat<-length(new.level)    # H: numcat: number of category other than reference
      
      cat.tab=NULL
      for(j in 1:numCat){
        cat.tab<-rbind(cat.tab,c(NA,paste(new.level[j]," vs.",reflevel[i],sep="")))
      }
      cat.tab[1,1]<-c(as.character(label[i])) # H: Insert variable label in the top left corner. 
      Col1<-rbind(Col1,cat.tab)
    }
    if (type[i]=="num"){  
      num.var<-c(as.character(label[i]),NA)
      Col1<-rbind(Col1,num.var)
    }
  }
  Table<-cbind(Col1,Table)
  
  rownames(Table) <- NULL
  Table[is.na(Table)]<-rep('',sum(is.na(Table)))

  res<-Table
  # P values NEJM format -----------------------------------------------
  tmp <- tmpc <- as.numeric(res[,"P-value"])
  tmp[tmpc >= 0.005 & !is.na(tmpc)] <- formatC(tmpc[tmpc >= 0.005 & !is.na(tmpc)], digit=2, format="f")
  tmp[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)] <- formatC(tmpc[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)], digit=3, format="f")
  tmp[tmpc >=0.995] <- ">0.99"
  tmp[tmpc < 0.001] <- "<.001"
  tmp[tmpc < 0.05 & !is.na(tmpc)]  <- as.character(paste("**",tmp[tmpc < 0.05 & !is.na(tmpc)],"**",sep=""))  ## added 01/04/2019
  res[,"P-value"] <- tmp
  
  return(res)
}



#' @title Create Multivariable Linear Regression Table 
#'
#' @description The function creates a result table of multivariable linear regression.
#' @param outcome a continous outcome variable
#' @param covariate a vector of variables included in the multivariable linear analysis
#' @param label a vector of labels shown in the table for covariates
#' @param type a vector of types either categorical or numerical ("\code{cat}" or "\code{num}")
#' @param reflevel a vector of reference levels for categorical variables. \code{NA} is used for continuous variable.
#' @param data data frame dataset
#' @export
#' @importFrom stats glm
#' @return matrix with estimates of coefficients and P-value
#' @examples 
#' # Use mtcars dataset for example
#' df <- mtcars
#' df$mpg.cat <- ifelse(df$mpg > 18, 1, 0)
#' 
#' # Create variable, label, type, and reference information
#' r1 <- c("hp", "Horse power", "num", NA)
#' r2 <- c("wt", "Weight", "num", NA)
#' r3 <- c("gear", "Gear", "cat", "3")
#' input  <- rb(3)
#' 
#' mul.linear.reg(outcome="mpg", 
#'                covariate=input[,"variable"],
#'                label=input[,"label"],
#'                type=input[,"type"],
#'                reflevel=input[,"ref"], 
#'                data=df
#'                )   
#'                               


mul.linear.reg <- function(outcome,covariate,label,type,reflevel=NULL,data){
  
  data <- na.handle.2(data, variable=covariate, type, useNA=FALSE)  
  
  digit=2
  label <- as.character(label)
  
  numvar=length(covariate)
  
  #reflevel: vector same length as covariate, label, and type. For continuous variables, assign NA value.
  
  for (i in 1:numvar){
    if (type[i]=="cat"){  
      #Assign reference levels for categorical variables
      data[,which(colnames(data)==covariate[i])]<-as.factor( data[,which(colnames(data)==covariate[i])])
      data[,which(colnames(data)==covariate[i])]<- relevel(data[,which(colnames(data)==covariate[i])], ref = as.character(reflevel[i])) # HH: Added as.character
    }
  }
  
  #Create outcome formula
  func <- paste(outcome,"~ ")
  
  fmla <- as.formula(paste(func, paste(covariate, collapse = "+")))
  
  model <- glm(fmla, family=gaussian(link="identity"), data=data, na.action = na.omit)
  # HH: Fixed data=na.omit(data) to avoid whole na data frame.
  
  x <- summary(model)
  
  coef<-round(coef(x)[-1,1],digit)
  #OR<-round(exp(coef(x)[-1,1]),digit)
  se<-round(coef(x)[,2],digit)[-1]
  pvalues<-round(coef(x)[,4],4)[-1]
  
  #Confidence interval for the coefficients
  CI.LL.beta<-coef-1.96*se
  CI.UL.beta<-coef+1.96*se
  interval.beta<-paste('(',round(CI.LL.beta,digit),';',round(CI.UL.beta,digit),')',sep='')
  
  res<-cbind(coef, interval.beta, pvalues)
  colnames(res)<-c("Coefficient", "95% CI ", "P-value")
  
  #Create labels for the first/second columns of the table
  Col1<-NULL
  for (i in 1:numvar){
    
    if (type[i]=="cat"){  
      catvar <- data[,which(colnames(data)==covariate[i])]
      catvar <- catvar[!is.na(as.character(catvar)) & as.character(catvar)!=""]  #HH: added additional condition to avoid "". 
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
    colnames(Col1) <- c("Variable", "")
  }
  
  res<-cbind(Col1,res)
  
  rownames(res) <- NULL
  # res[is.na(res)]<-rep('',sum(is.na(res)))  # Just change NA to "" Not necessary
  
  # P values NEJM format -----------------------------------------------
  tmp <- tmpc <- as.numeric(res[,'P-value'])
  tmp[tmpc >= 0.005 & !is.na(tmpc)] <- formatC(tmpc[tmpc >= 0.005 & !is.na(tmpc)], digit=2, format="f")
  tmp[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)] <- formatC(tmpc[tmpc < 0.005 & tmpc >= 0.001 & !is.na(tmpc)], digit=3, format="f")
  tmp[tmpc >=0.995] <- ">0.99"
  tmp[tmpc < 0.001] <- "<.001"
  tmp[tmpc < 0.05 & !is.na(tmpc)]  <- as.character(paste("**",tmp[tmpc < 0.05 & !is.na(tmpc)],"**",sep=""))  ## added 01/04/2019
  res[,'P-value'] <- tmp
  
  res[is.na(res)] <- ""
  return(res)
}


