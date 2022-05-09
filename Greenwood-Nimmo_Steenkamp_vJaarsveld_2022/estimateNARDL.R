#======================================================================================
# Generalized NARDL code
#======================================================================================
# NB NB when converting to zoo you need to do it like this!!!
# !!! dep <- as.zoo(data_est[,1, drop=F])  !!!
# This function is dependent on dynlm and stringr

# Check what hap[pens with first obs in data]

estimate_nardl <- function(intercept=TRUE, dep  , indep  , decomp  , p , qp, qn , r_mat = as.matrix(c(1), row=T)){
  # IN function testing 
  #p=c(1:p_lag); qp=c(0:q_lag);qn=c(0:q_lag); r_mat=rmat
  #  Build the formula
  #---------------------------------------------------------------------------   
  # 1) Fist get the decompositions
  
  xp <- list()
  xn <-  list()
  l_n <- 1
  for (i in 1:ncol(indep)) {
    if(colnames(indep)[i] %in% decomp){
      diff_x <- diff(indep[,i])
      
      diff_x_p <- ((diff_x > 0)*diff_x) # Conceptually should perhaps be an NA
      
      xp_series <- cumsum(na.omit(diff_x_p))
      xp_series <- c(0,xp_series)
      
      xp[[l_n]] <- xp_series
      
      diff_x_n <- ((diff_x < 0)*diff_x)
       
      xn_series <- cumsum(na.omit(diff_x_n))
      xn_series <- c(0,xn_series)
      xn[[l_n]] <- xn_series
      
      
      l_n <- l_n+1
    } # end if (colnames(indep)[i] %in% decomp)
  }# end for i
  
  # place all the xp and xn in a single zoo object
  xp <- do.call(cbind, xp); colnames(xp) <- paste(decomp, "p", sep="") 
  xn <- do.call(cbind, xn); colnames(xn) <- paste(decomp, "n", sep="")
  
  # Convert xp and xn to zoo objects
  xp <- zoo(xp,order.by = index(indep))
  xn <- zoo(xn,order.by = index(indep))
  
  #---------------------------------------------------------------------------   
  # 2) Build the decomposition part of the formula 
  f_decomp <- c(colnames(xp),colnames(xn))
  f_decomp <- paste("L(", f_decomp,",1)", sep = "", collapse = "+")

  
  
  #---------------------------------------------------------------------------    
  # 3) Create lags for the dependent variable of the formula (p)
  
    temp <- list()
    l_n <- 1 
    for(i in p){ 
      temp[[l_n]] <- paste("L(d(",colnames(dep),"),",i, ")",sep = "")
      l_n <- l_n+1
      }
    f_dep_L <- paste(unlist(temp),collapse="+")
  
  #---------------------------------------------------------------------------    
  # 4) Create lags for the decomposition terms of the formula (q)
  decomp_terms <- c(colnames(xp),colnames(xn))
    
    temp1 <- list()
    l_n <- 1
    for (i in qp) {
      temp1[[l_n]] <- paste("L(d(",decomp_terms[1],"),",i,")", sep="")
      l_n <- l_n+1
    }
    
    temp2 <- list()
    l_n <- 1
    for (i in qn) {
      temp2[[l_n]] <- paste("L(d(",decomp_terms[2],"),",i,")", sep="")
      l_n <- l_n+1
    }
    
    
    f_decomp_L <- paste(paste(unlist(temp1),collapse = "+"), paste(unlist(temp2), collapse = "+"), sep = "+")
    
  #---------------------------------------------------------------------------    
  # 5) Create lags for the independent variables that are not decomposed (r)
    
    # 5.1) First get a single lag of the non-decomposed variables
    l_n <- 1
    temp <- list()
    for (i in 1:ncol(indep)) {
      if(!(colnames(indep)[i] %in% decomp)){
         temp[[l_n]] <- paste("L(",colnames(indep)[i], ",1)", sep = "")
         l_n <- l_n+1
      } # end if (colnames(indep)[i] %in% decomp)
    }# end for i   
    f_indep_sL <- paste(unlist(temp), collapse = "+")
    
    # 5.2) Create lag structure for the non-decomposed differenced independent variables
    # First check if there are independent variables that should not be decomposed
    indep_temp <- indep[,!colnames(indep)%in%decomp]
    
    l_n <- 1
    temp <- list()
    for (i in 1:(ncol(indep_temp))) {
      
        for (j in 1:ncol(r_mat)) { # starts at first row in matrix. Columns have the lags in while rows are the variables
         temp[[l_n]] <- paste("L(d(",colnames(indep_temp)[i], "),", r_mat[i,j], ")", sep = "") 
         l_n <- l_n+1
        } # end for (j in 1:(r-1)) 
        
      
    }# end for i   
    f_indep_L <- paste(unlist(temp), collapse = "+")

    # Do an operation with r so that if there are independent variables that are not 
    # getting decomposed the dunction will not return an error
    decomp_chk <- str_detect(colnames(indep), decomp)
    if (TRUE %in% decomp_chk){
      r <- 0
    }

    #---------------------------------------------------------------------------    
    # 5) Build the formula
      # First check if there are independent variables that are not decomposed
    if (intercept==TRUE){
    decomp_chk <- str_detect(colnames(indep), decomp)
    if(FALSE %in% decomp_chk){
    f_combined <- paste("d(",colnames(dep),")", " ~ 1 +", "L(",colnames(dep),",1)", sep="")
    f_combined <- c(f_combined,
                       f_decomp,
                       f_dep_L,
                       f_decomp_L,
                       f_indep_sL,
                       f_indep_L)
    f_combined <- as.formula(paste(f_combined, collapse = "+")) 
    } else{
      f_combined <- paste("d(",colnames(dep),")", " ~ 1 + ", "L(",colnames(dep),",1)", sep="")
      f_combined <- c(f_combined,
                      f_decomp,
                      f_dep_L,
                      f_decomp_L) 
      f_combined <- as.formula(paste(f_combined, collapse = "+")) 
    }
    } # end intercept true
    
    if (intercept==FALSE){
      decomp_chk <- str_detect(colnames(indep), decomp)
      if(FALSE %in% decomp_chk){
        f_combined <- paste("d(",colnames(dep),")", " ~ 0 +", "L(",colnames(dep),",1)", sep="")
        f_combined <- c(f_combined,
                        f_decomp,
                        f_dep_L,
                        f_decomp_L,
                        f_indep_sL,
                        f_indep_L)
        f_combined <- as.formula(paste(f_combined, collapse = "+")) 
      } else{
        f_combined <- paste("d(",colnames(dep),")", " ~ 0 + ", "L(",colnames(dep),",1)", sep="")
        f_combined <- c(f_combined,
                        f_decomp,
                        f_dep_L,
                        f_decomp_L) 
        f_combined <- as.formula(paste(f_combined, collapse = "+")) 
      }
    } # end intercept FALSE
    # B) Combine data and place it all into a single zoo object 
    # Need to fix this part to be more general once again need to be able to handle more than one xp and xn
    #xp <- ts(xp, start = c(start(indep)[1],start(indep)[2]), end = c(end(indep)[1],end(indep)[2]), frequency = 12)
    #xn <- ts(xn, start = c(start(indep)[1],start(indep)[2]), end = c(end(indep)[1],end(indep)[2]), frequency = 12)
    data <- cbind(dep,indep,xp,xn)
    
    colnames(data) <- c(colnames(dep), colnames(indep),paste(decomp, "p", sep="") ,paste(decomp, "n", sep="") )
    model_est <- dynlm(f_combined, data = data)
    return(model_est)
  
}
