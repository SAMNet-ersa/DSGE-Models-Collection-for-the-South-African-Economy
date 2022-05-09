# The GETS script

NARDL_auto_lag <- function(dep, indep, decomp, p_max, q_max, r_max=0, selection="AIC"){
  
  # Create empty matrix to store IC results in 
  ic_matrix <- matrix(data=0, nrow = q_max, ncol = p_max)
 
  # Obtain a first estimate as a benchmark to check  
  p <- c(1,2)
  q <- c(0,1)
  model_est <- estimate_nardl(dep = dep, indep = indep, decomp = decomp, p=p, qp=q , qn=q)
  model_est_AIC <- AIC(model_est)
  model_est_BIC <- BIC(model_est)
  
  # Loop for when there are non-decomposed independent variables
  if(r_max==0){
  for (i in 2:p_max) {
    for (j in 2:q_max){
      # Estimate the model
        p <- c(1,seq(2,i))
        q <- c(0,seq(1,(j-1)))
        model_est_chk <- estimate_nardl(dep = dep, indep = indep, decomp = decomp, p = p, qp=q , qn=q)  
        
        model_est_chk_AIC <- AIC(model_est_chk)
        model_est_chk_BIC <- BIC(model_est_chk)
        
        ic_matrix[i,j] <- model_est_chk_AIC 
        
        if((model_est_chk_AIC < model_est_AIC) & selection == "AIC"){ 
          model_est <- model_est_chk 
          model_est_AIC <- model_est_chk_AIC
          }
        
        if((model_est_chk_BIC < model_est_BIC) & selection == "BIC"){ 
          model_est <- model_est_chk 
          model_est_BIC <- model_est_chk_BIC
        }
        
        if(selection == "AIC"){ic_matrix[i,j] <- model_est_chk_AIC}
        if(selection == "BIC"){ic_matrix[i,j] <- model_est_chk_BIC}
      
    }# end for j
  }# end for i
  } # end if r==0
  
  if(r_max!=0){
    # Create a first model to benchamrk against
    p <- c(1,2)
    q <- c(0,1)
    r <- c(1,2)
    r_mat <- matrix(data=rep(seq(1:r_max),(ncol(indep)-1)),
                    nrow = (ncol(indep)-1), ncol = r_max, byrow = TRUE)
    
    model_est <- estimate_nardl(dep = dep, indep = indep, decomp = decomp,  p=p, qp=q , qn=q, r_mat = r_mat)
    model_est_AIC <- AIC(model_est)
    model_est_BIC <- BIC(model_est)
    
    ic_matrix <- array(data=NA,dim = c(p_max,q_max,r_max))
    
    for (i in 2:p_max) {
      for (j in 2:q_max){
        for (k in 2:r_max){
        # Estimate the model
        p <- c(1,seq(2,i))
        q <- c(0,seq(1,(j-1)))
        r_mat <- matrix(data=rep(seq(1:k),(ncol(indep)-1)),
                        nrow = (ncol(indep)-1), ncol = k, byrow = TRUE)
        
        model_est_chk <- estimate_nardl(dep = dep, indep = indep, decomp = decomp, p=p, qp=q , qn=q, r_mat = r_rmat)  
        
        model_est_chk_AIC <- AIC(model_est_chk)
        model_est_chk_BIC <- BIC(model_est_chk)
        
        
        
        if((model_est_chk_AIC < model_est_AIC) & selection == "AIC"){ 
          model_est <- model_est_chk 
          model_est_AIC <- model_est_chk_AIC
        }
        
        if((model_est_chk_BIC < model_est_BIC) & selection == "BIC"){ 
          model_est <- model_est_chk 
          model_est_BIC <- model_est_chk_BIC
        }
        
        if(selection == "AIC"){ic_matrix[i,j,k] <- model_est_chk_AIC}
        if(selection == "BIC"){ic_matrix[i,j,k] <- model_est_chk_BIC}
        
        } # end for r  
      }# end for j
    }# end for i
  } # end if r==0
  
  
  #print(ic_matrix)
  return(list(model_est,ic_matrix))
  
} # end of function
