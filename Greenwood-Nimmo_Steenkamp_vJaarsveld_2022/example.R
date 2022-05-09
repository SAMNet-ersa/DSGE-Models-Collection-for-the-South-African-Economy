# Main scripts for example------------------------------------------------------------------

# Package and function load-----------------------------------------------------------------
if (!require('pacman')) install.packages('pacman'); library('pacman')
p_load(dynlm, stringr, zoo)

source("./estimateNARDL.R")
source("./NARDL_auto_lag.R")
source("./mplier_base.R")


# Import USA data---------------------------------------------------------------------------

  USdata <- readRDS("./USdata.rds")
  
# Estimate NARDL----------------------------------------------------------------------------
  
# Convert data to zoo object
  USdata_zoo <- as.zoo(USdata[,-1], order.by = USdata$date)  

# Set function parameters
  dep <- USdata_zoo[,"UN", drop=FALSE]
  indep <- USdata_zoo[,c("IP"), drop=FALSE]
  decomp <- "IP"

# Estimate NARDL----------------------------------------------------------------------------

  nardl <-   estimate_nardl(dep = dep, indep = indep, decomp = decomp, 
                            p=c(1,11), qn=c(0,4), qp=c(0,2))
  summary(nardl)
  
# Auto lag
  NARDL_auto <- NARDL_auto_lag(dep = dep, indep = indep, decomp = decomp, p_max = 5,q_max = 5)
  
# Obtaining the dynamic multipliers---------------------------------------------------------
  
  mplier <-  mplier_base(nardl, dep = dep, decomp = decomp, k=2,l=12,h=80)  
  mplier <- rbind(c(0,0),
                  mplier)
  
  plot(mplier[,1],ylim = c(-20,30), type = "l", ylab = "dynamic multiplier", xlab = "h")  
  lines(-mplier[,2])
  
  
