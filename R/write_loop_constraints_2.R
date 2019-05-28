#'\code{write_loop_constraints_2}
#'
#'@param variables Contains the list of variables as used to formulate the ILP problem, explanations for each variable and a list of useful indices.
#'@param loops Contains the list of identified loops in the PKN
#'
#'@return This function writes the constraints preventing self-activation of nodes in the network due to positive feedback loops.
#'
#'@export

write_loop_constraints_2 <- function(variables=variables, loops = loops, pknList = pknList){
  
  library(igraph)
  gg = graph_from_data_frame(d = as.data.frame(pknList[, c(1, 3)]), directed = TRUE)
  adj = get.adjacency(graph = gg)
  
  if(!is.character(loops[[1]])){
    
    constraints = c()
    for(ii in 1:length(loops)){
      
      ll = loops[[ii]]
      
      for(vv in 1:length(variables)){
        
        var = variables[[vv]]
        for(jj in 1:(length(ll)-1)){
          
          ss = rownames(adj)[ll[jj]]
          tt = rownames(adj)[ll[jj+1]]
          
          if(jj==1){
            
            c1 = var$variables[which(var$exp==paste0("ReactionUp ", ss, "=", tt, " in experiment ", vv))]
            c2 = var$variables[which(var$exp==paste0("ReactionDown ", ss, "=", tt, " in experiment ", vv))]
            
          } else {
            
            c1 = paste0(c1, " + ", var$variables[which(var$exp==paste0("ReactionUp ", ss, "=", tt, " in experiment ", vv))])
            c2 = paste0(c2, " + ", var$variables[which(var$exp==paste0("ReactionDown ", ss, "=", tt, " in experiment ", vv))])
            
          }
          
        }
        
        c1 = paste0(c1, " < ", length(ll)-1)
        c2 = paste0(c2, " < ", length(ll)-1)
        
        constraints = c(constraints, c1, c2)
        
      }
      
    }
    
    return(constraints)
    
  } else {
    
    constraints = c()
    for(ii in 1:length(loops)){
      
      ll = loops[[ii]]
      ll = c(ll, loops[[ii]][1])
      
      for(vv in 1:length(variables)){
        
        var = variables[[vv]]
        for(jj in 1:(length(ll)-1)){
          
          # ss = rownames(adj)[ll[jj]]
          # tt = rownames(adj)[ll[jj+1]]
          ss = ll[jj]
          tt = ll[jj+1]
          
          if(jj==1){
            
            c1 = var$variables[which(var$exp==paste0("ReactionUp ", ss, "=", tt, " in experiment ", vv))]
            c2 = var$variables[which(var$exp==paste0("ReactionDown ", ss, "=", tt, " in experiment ", vv))]
            
          } else {
            
            c1 = paste0(c1, " + ", var$variables[which(var$exp==paste0("ReactionUp ", ss, "=", tt, " in experiment ", vv))])
            c2 = paste0(c2, " + ", var$variables[which(var$exp==paste0("ReactionDown ", ss, "=", tt, " in experiment ", vv))])
            
          }
          
        }
        
        c1 = paste0(c1, " < ", length(ll)-1)
        c2 = paste0(c2, " < ", length(ll)-1)
        
        constraints = c(constraints, c1, c2)
        
      }
      
    }
    
    return(constraints)
    
  }
  
}
