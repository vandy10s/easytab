#' @title Rowbind Vectors of Input
#'
#' @description The function combines vectors of input by row-wise, where vector is composed of variable, label, type, and reference level.
#' @param x Number of vectors to be combined
#' @export
#' @return A matrix with following columns c(variable, label, type, reference) or c(variable, label, type, reference, order)
#' @examples
#' r1 <- c("mpg","MPG","num",NA,1)
#' r2 <- c("gear","Gear","cat","3",2)
#' rb(2)

rb <- function(x){
  len    <- seq(1,x)
  rs     <- rep("r", length(len))
  rlists <- paste(rs,len,sep="")
  com    <- paste("rbind(", paste(rlists, collapse=","),")")
  mat    <- eval(parse(text=com))

  if (ncol(mat)==5){
    colnames(mat) <- c("variable", "label", "type", "ref", "sort")
    mat <- mat[order(as.numeric(mat[,5])),]
    mat[,5] <- seq(1,nrow(mat))
  } else if (ncol(mat)==4) {
    colnames(mat) <- c("variable", "label", "type", "ref")
  }

  mat = (mat[complete.cases(mat[,1]),])
  return(mat)
}








#' @title Median Calculation with Min and Max
#'
#' @description Calculate median, min and max for continuous variables
#' @param var is a vector of continuous variable
#' @param digits defines a decimal point
#' @return A character vector with median, min, and max
#' @keywords internal

getMedian <- function(var,digits=2){
  output <- rep(NA,2)
  output[1] <- paste('N = ',sum(!is.na(var)),sep="") # Total N for continuous variables
  output[2] <- paste(round(median(var,na.rm=TRUE),digits)," (",
                     round(min(var,na.rm=TRUE),digits),", ",
                     round(max(var,na.rm=TRUE),digits),")",sep="")
  output
}



#' @title Median Calculation with Q25 and Q75
#'
#' @description Calculate median, q25, and q75 for continuous variables
#' @param var is a vector of continuous variable
#' @param digits defines a decimal point
#' @return A character vector with median, q25, and q75
#' @keywords internal

getMedian.mid <- function(var,digits=2){
  output <- rep(NA,2)
  output[1] <- paste('N = ',sum(!is.na(var)),sep="") # Total N for continuous variables
  output[2] <- paste(round(quantile(var,na.rm=TRUE)[3],digits)," (",
                     round(quantile(var,na.rm=TRUE)[2],digits),", ",
                     round(quantile(var,na.rm=TRUE)[4],digits),")",sep="")
  output
}



#' @title Mean calculation with Standard Deviation
#'
#' @description Calculate mean and sd for continuous variables
#' @param var is a vector of continuous variable
#' @param digits defines a decimal point
#' @return A character vector with sd
#' @keywords internal

getMean <- function(var,digits=2){
  output <- rep(NA,2)
  output[1] <- paste('N = ',sum(!is.na(var)),sep="") #Total N for continuous variables
  output[2] <- paste(round(mean(var,na.rm=TRUE),digits)," (",
                     round(sqrt(var(var,na.rm=TRUE)),digits),")",sep="")
  output
}


#' @title Percentages Calculation with Frequencies
#'
#' @description Calculate percentages for categorical variables with frequencies
#' @param var is a categorical variable
#' @param category is a list of levels in categorical variable
#' @param digits defines a decimal point
#' @return A character vector with sd
#' @keywords internal
#'
getNP <- function(var,category,digits=2){
  freqTab <- table(var, useNA="no")
  propTab <- prop.table(freqTab)
  outFreq <- as.numeric(freqTab[names(freqTab)==category])
  outPct <- round(as.numeric(propTab[names(propTab)==category])*100,digits)
  output <- paste(outFreq," (",outPct,"%)",sep="")
  return( output)
}




#' @title Transform "" levels as NA in factor
#'
#' @description Calculate percentages for categorical variables with frequencies
#' @param data is a dataset
#' @param variable is a variable
#' @param type is a type of that variable
#' @param useNA is whether or not to include NULL
#' @param na is how to express NULL values
#' @export
#' @return A updated vector with either NULL value or whatever values for NULL values such as "Unknown"
#' @keywords internal
#'

na.handle.2 <- function(data, variable, type, useNA, na){

  numvar <- length(variable)
  ## Exclusion list from the levels
  exclusion <- c(""," ","  ","N/A","n/a","NA","na")

  for (i in 1:numvar){

    var <- data[,which(colnames(data)==variable[i])]

    if (type[i]=="cat"){

      ## Save the order of levels if class==FACTOR
      if (class(var)=="factor"){
        lev <- levels(var)
      }else{
        lev <- levels(as.factor(var))
      }

      lev <- setdiff(lev, exclusion)      # lev <- lev[!(lev %in% exclusion)]

      ## Transform factor into character to remove unnecessary levels
      tmp <- as.character(var)

      if (useNA==T){
        tmp[tmp %in% exclusion|is.na(tmp)] <- na
        tmp <- as.factor(tmp)
        last <- na
        new_levels = c(setdiff(lev, last), last)  # make NA the last level!
        var <- factor(tmp, levels = new_levels)
      } else {
        tmp[tmp %in% exclusion] <- NA
        var <- factor(tmp, levels = lev)
      }
    }
    data[,which(colnames(data)==variable[i])] <- var
  }
  return(data)
}



