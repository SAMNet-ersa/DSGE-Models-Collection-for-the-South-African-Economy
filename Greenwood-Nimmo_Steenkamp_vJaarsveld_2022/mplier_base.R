# Function that creates the dynamic multipliers
#===============================================================================

mplier_base <- function(model,dep,decomp,k,l,h){
  # In function testing parameters
  #model <-  nardl
  
  
  #model=mod_boot;k = 2;l = 5;h = 24
  
  
# Function needs the 1) ls model object ; 2) name of dep; 3) and indep 
  # First create a vector that contains all the names 
  dep_name <- colnames(dep)
  
  var_names <- names(coefficients(model))
  if(var_names[1]=="(Intercept)"){
    var_names <- var_names[-1]
  }
  var_names <- as.data.frame(matrix(unlist(str_split(var_names,", ")),ncol=2,byrow = TRUE))  # Format to have matrix with all the names to obtain indexing for parameters  
  var_names[,1] <- as.character(var_names[,1])
  var_names[,2] <- as.character(var_names[,2])
  
  # Need to check if there are some ` ` shit in there from boot function and then remove them all 
  if(any(str_detect(var_names[,1],"`") == TRUE)){
    var_names[,1] <- str_remove(var_names[,1] ,"[`]")
    var_names[,2] <- str_remove(var_names[,2] ,"[`]")
  }
  
  var_names[,2] <- as.numeric(str_remove(var_names[,2] ,"[)]")) # Remove the ")" so that the lag orders can be used for indexing 
  
# Find indicies---------------------------------------------------------------
  
  a_index <- which(var_names == paste("L(", dep_name, sep="")) # Index the first lagged dependent in the matrix
  b_index <- which(var_names == paste("L(", decomp,"p",sep = "") | var_names == paste("L(",decomp,"n",sep = "")) # Lagged values of decomp
  
  # Need to read lags from varnames col 
  varphi_index <- which(var_names == paste("L(d(", dep_name,")", sep = "")) 
  varphi_lags <- var_names[varphi_index,2]
  
  vpi_index_p <-  which(var_names == paste("L(d(", decomp,"p",")", sep = ""))
  vpi_lags_p <- var_names[vpi_index_p,2]
  vpi_index_n <- which(var_names == paste("L(d(", decomp,"n",")", sep = ""))
  vpi_lags_n <- var_names[vpi_index_n,2]  
  
# Assign values to mplier variables-------------------------------------------
  model_params <- model$coefficients
  if(names(model_params)[1]=="(Intercept)"){
    model_params <- model_params[-1]
  }
  
  alpha <- model_params[a_index] 
  beta <- model_params[b_index]
  
  varphi <- rep(0,l-1)
  varphi[varphi_lags] <- model_params[varphi_index]
  
  vpi <- matrix(data = 0, nrow = l, ncol = 2)
  vpi[(vpi_lags_p)+1,1] <- model_params[vpi_index_p]
  vpi[(vpi_lags_n)+1,2] <- model_params[vpi_index_n] # Continue from here 

# Calculate the DM-------------------    
  vphi <- matrix(data=0, nrow = 1, ncol = l)
  vphi[1] = 1 + alpha + varphi[1]
  
  for (i in 2:(l-1)) {
    vphi[i] <- varphi[i] - varphi[i-1]
  }
  vphi[l] <- -varphi[l-1]
  
  mtheta <- matrix(data=0, nrow = (l+1), ncol = k)
  mtheta[1,] <- vpi[1,]
  mtheta[2,] <- vpi[2,] - vpi[1,] + beta
  
  for (i in 3:l) {
    mtheta[i,] <- vpi[i,] - vpi[i-1,]
  }
  mtheta[l+1,] <- -vpi[l,]
  
  mpsi <- matrix(data=0, nrow = h, ncol = k)
  mpsi[1,] <- mtheta[1,]
  
  for (i in 1:l) {
    mpsi[i+1,] <- vphi[,1:i]%*%mpsi[i:1,] + mtheta[i+1,]
  }
  
  for (i in (l+1):(h-1)) {
    mpsi[i+1,] <- vphi%*%mpsi[i:(i-l+1),]
  }
  mpsi <- apply(mpsi, 2, cumsum)
  return(mpsi)
} # end of function mplier
