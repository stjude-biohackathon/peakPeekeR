#' Test function
#'
#' Does nothing but test that we can run peakcallers from different virtual environments.
#'
#' @return A list of names of objects exposed in each module.
#' @author Jared Andrews
#'
#' @examples
#' test()
#' @export
#' @importFrom basilisk basiliskStart basiliskRun basiliskStop BasiliskEnvironment

test <- function() {
  cl <- basiliskStart(env_macs2)
  macs2.v <- basiliskRun(cl, function() {
    
    outtie <- system2('macs2', args = c('--version'), stdout = TRUE)
    
    return(outtie)
  })
  basiliskStop(cl)
  
  cl <- basiliskStart(env_macs)
  macs.v <- basiliskRun(cl, function() {
    
    outtie <- system2('macs', args = c('--version'), stdout = TRUE)
    
    return(outtie)
  })
  basiliskStop(cl)
  
  cl <- basiliskStart(env_sicer2)
  sicer2.v <- basiliskRun(cl, function() {
    
    outtie <- system2('sicer', args = c('--version'), stdout = TRUE)
    
    return(outtie)
  })
  basiliskStop(cl)
  
  list(macs2.v, macs.v, sicer2.v)
}
