library(pso) #The original library that implements the Particle Swarm Optimization with many functions
source("psoptim_edited.R") #Our modified version of the psoptim function for the pso library
require(ConsRank)



#----------------------------------------
# Function useful to visualize one consensus ranking or to visualize and compare two of them
# X = matrix or dataset with the full data (used for the items name)
consensus_comparison = function(cr1, cr2=NULL, X){
  if(is.null(cr2)){
    return(t(rank2order(cr1, items=colnames(X))))
  }
  else{
    return(t(rbind(rank2order(cr1, items=colnames(X)), rank2order(cr2, items=colnames(X)))))
  }
}

##---------
breakties<-function(X){
  X<-rank2order(X)
  open<-"{"
  close<-"}"
  X<-gsub(open,"", X,fixed=TRUE)
  X<-gsub(close,"", X,fixed=TRUE)
  X<-order2rank(X)
  return(X)
}
#----------
# Function that estimate a good starting point to inizialize one particle. This function comes from the ConsRank package
findconsensusBB2 <- function(cij,FULL=FALSE) {
  N <- ncol(cij)
  X <- matrix(1,1,N)
  #I need the id later on
  idx <- seq(1:N)
  #better immediately convert cij in sign matrix
  scij <- sign(cij)
  #now, process each item (insetad of each pair of items)
  for(j in 1:N){
    #case1:
    #   if sign(cij)=1 & sign(cji)=-1, then X[i]=X[i]+1
    a <- sum( (scij[j,-seq(1:j)]==1) * (scij[-seq(1:j),j]==-1))
    X[j] <- X[j]+a
    #
    #case 2: 
    #   if sign(cij)=-1 & sign(cji)=1, then X[j]=X[j]+1
    b <- ( (scij[j,-seq(1:j)]==-1) * (scij[-seq(1:j),j]==1))
    indb=idx[-seq(1:j)][which(b==1)]
    X[indb] <- X[indb]+1 
    #
    #case 3:
    #  if sign(cij)=1 & sign(cji)=1, then X[i]=X[i]+1 and X[j]=X[j]+1
    c=( (scij[j,-seq(1:j)]==1) * (scij[-seq(1:j),j]==1))
    indc=idx[-seq(1:j)][which(c==1)] 
    #indc is the index of only X[j]'s
    X[indc] <- X[indc]+1
    X[j] <- X[j]+length(indc)
  }
  
  X<-(N+1)-X
  
  if (FULL==TRUE){
    X<-breakties(X)
  }
  
  return(X)
}

# Core function that implements the PSOPR. Input parameters name are shared with those of the ConsRank package for usability purposes, thus being user-friendly 
#INPUT:
# X: A n by m data matrix, in which there are n judges and m objects to be judged. 
#      Each row is a ranking of the objects which are represented by the columns. If X contains the rankings observed only once, the argument Wk can be used
# Wk: Optional: the frequency of each ranking in the data
# maxiter: maximum number of iterations. maxiter=1000 is the default option.

# NP: the number of population individuals (i.e. particles). NP=15 is the default option.
# L:  (this is just a shared parameter with ConsRank, equals to maxiter) generations limit, maximum number of consecutive generations without improvement. L=100 is the default option. Same of maxiter (inserted for compatibility with FASTDECOR)
# FF: inertia_weight: a constant to control the impact of the previous velocity on the new velocity. FF=0.4 is the default option. (inserted for compatibility with FASTDECOR)
# CR: the crossover range. Must be in [0,1]. CR=0.9 is the default option. (inserted for compatibility with FASTDECOR)
# PS: If PS=TRUE, on the screen some information about how many branches are processed are displayed.

# phi_p: cognitive_component constant, a term that directs the particle towards its personal best position. Default = CR
# phi_g: social_component constant, a term that directs the particle towards the global best position. Default = CR
# REPORT=5: report best achieved solution after REPORT iterations
# maxit.stagnate=5: STOP after maxit.stagnate consecutive iterations with the same best achieved solution
# lambda=NULL: deprecated, penalization term for continuous data (we work on integers now)
# cij = NULL: Optional: conbined input matrix
# reltol = default(0) means no restart. Restarting tolerance if the algorithm is stagnating
# FULL = F: allows to unfold the ties solutions founded into solutions without ties. This process can be slow and not useful
# FULL_simple: instead of providing all the possible without ties solutions, remove the ties to the solution founded.
# FULL_limit = 100: upper bound of possible without ties equivalent solutions that can be unfolded.

#OUTPUT
#A list containing:
# The consensus ranking(s)
# The associated Tau_x
# The elapsed time in second

rankagg_pso <- function(X, Wk = NULL, maxiter = L, NP = 15, L = 100, FF = 0.4, CR = 0.9, PS = FALSE, phi_p = CR, phi_g = CR, REPORT=5, maxit.stagnate=5, lambda=NULL, cij = NULL, reltol=0, FULL=FALSE, FULL_limit = 100, FULL_simple=T) { #, FULL = FALSE
  force(maxiter) #force to set maxiter = L, if maxiter have the defaul value
  force(phi_p)
  force(phi_g)
  m <- nrow(X)
  n <- ncol(X)
  
  if(is.null(cij)){
    if(is.null(Wk)){
      cij = combinpmatr(X)
    } 
    else{
      cij = combinpmatr(X, Wk)
    }
  }
  
  
  
  # Define the objective function for PSO
  obj_fun3 <- function(cr, cij, m, n, Wk=NULL, lambda=NULL) {
    
    cr = round(cr)
    
    if(!is.null(lambda)){
      lambda_terms = abs(cr - round(cr))
    }
    
    #cr = round(cr)
    
    Sij = scorematrix(cr)

    if(is.null(Wk)){
      TauX = sum(cij*Sij)/(m*n*(n-1)) #dot(cij, Sij)
    }
    else{
      TauX = sum(cij*Sij) / ( sum(Wk)* (n*(n-1)) )
    }
    
    if(!is.null(lambda)){
      TauX = TauX - sum((lambda*lambda_terms/0.5)) #fattore di penalizzazione in proporzione a quanto lontano dall'intero
    }
    
    return(-TauX)
  }
  
  # Define the bounds for PSO
  bounds <- matrix(c(rep(1, n), rep(n, n)), ncol = 2)
  
  
  # Inizialize the parameters
  swarm_size = NP
  omega = FF 
  
  #find a very good candidate
  starting_particle = as.vector(findconsensusBB2(cij))
  
  # Run PSO to find the consensus ranking
  start_time = Sys.time()
  pso_result <- psoptim(starting_particle, obj_fun3, cij = cij, m=m, n=n, Wk = Wk, lambda=lambda, lower = bounds[, 1], upper = bounds[, 2],
                        control=list(trace=as.integer(PS), maxit = maxiter, s=swarm_size, c.p=phi_p, c.g=phi_g, w=omega, maxit.stagnate = maxit.stagnate, abstol=-1, REPORT=REPORT, reltol=reltol)
                        )
  
  end_time = Sys.time()
  # Extract the consensus ranking and tau_x
  consensus_ranking <-  reordering(matrix(round(pso_result$par), nrow=length(pso_result$par)/n, byrow = T))
  
  
  #duplicated removal
  dupe_check = matrix(NA, nrow=nrow(consensus_ranking), ncol=ncol(consensus_ranking))
  for(i in 1:nrow(consensus_ranking)){
    dupe_check[i,] = t(consensus_comparison(consensus_ranking[i,], X=X))
  }
  dupe_check = !duplicated(dupe_check)
  consensus_ranking = consensus_ranking[dupe_check,]
  
  if(!is.matrix(consensus_ranking)){
    consensus_ranking = matrix(consensus_ranking, ncol=ncol(X))
  }
  
  tauX <- -pso_result$value[dupe_check]  # we invert the negative of tau_x to get the true value
  
  # Compute elapsed time
  elapsed_time <- end_time-start_time #pso_result$elapsed
  
  if(FULL == T){
    if(!FULL_simple){
        warning("FULL NO ties consensus are requested (FULL_simple = FALSE). Proceeding to unfold consesus identified with ties. This process can be slow and will not be computed in the elapsed time.")
        recursive_split = function(cr, to_return=NULL, max_solutions=Inf, counter=0){
          if(is.null(to_return)){
            to_return = matrix(NA, nrow=1, ncol=ncol(cr), byrow = T)
          }
          ties_presents = grepl("\\{|\\}", cr)
          if(sum(ties_presents)>0){
            ties = which(ties_presents)[1:2]
            
            if(counter < max_solutions){
              amount = (ties[2]-ties[1])+1
              permutated = gtools::permutations(amount, amount)
              cr = t(as.matrix(c(gsub("\\{|\\}", "", cr[,1:ties[2]]), cr[,(ties[2]+1):ncol(cr)])))
              
              return(rbind(recursive_split(t(as.matrix(cr[,c(1:(ties[1]-1),(ties[1]:ties[2])[permutated[1,]],(ties[2]+1):ncol(cr) )])), to_return, max_solutions, counter+1),
                           recursive_split(t(as.matrix(cr[,c(1:(ties[1]-1),(ties[1]:ties[2])[permutated[2,]],(ties[2]+1):ncol(cr) )])), to_return, max_solutions, counter+1)
                          ) 
                    )
            }
            else{
              cr = gsub("\\{|\\}", "", cr)
              return(cr)
            }
    
          }
          else{
            return(cr)
          }
        }
        
        for(j in 1:nrow(consensus_ranking)){
          cr = consensus_ranking[j,]
          cr = rank2order(cr)
          
          if(j == 1){
            risultato = recursive_split(cr = cr, max_solutions = FULL_limit/10)
          }
          else{
            risultato = rbind(risultato, recursive_split(cr = cr, max_solutions = FULL_limit/10))
          }
          
        }
        
        risultato = risultato[!duplicated(risultato),]
        
        if(nrow(risultato) > FULL_limit){
          print(paste("Founded", nrow(risultato), "solutions but randomly considering only", FULL_limit))
          set.seed(1)
          risultato = risultato[sample(1:nrow(risultato), FULL_limit),]
        }
        
    }
    else{
      for(j in 1:nrow(consensus_ranking)){
        cr = consensus_ranking[j,]
        cr = rank2order(cr)
        
        if(j == 1){
          risultato = gsub("\\{|\\}", "", cr)
        }
        else{
          risultato = rbind(risultato, gsub("\\{|\\}", "", cr))
        }
        
      }
    }

      risultato_tau = rep(NA, nrow(risultato))
      
      for(j in 1:nrow(risultato)){
        #print(round(j/nrow(risultato), 3))
        cr = order2rank(risultato[j,])
        cr = cr[,order(cr)]
        risultato_tau[j] = -obj_fun3(cr, cij, m, n, Wk=Wk, lambda=lambda) #mean(tau_x(data_sim, cr))
      }
      
      best_cr = risultato_tau == max(risultato_tau)
    
    # Return the results as a list
    return(list(Consensus = risultato[best_cr,], Tau = risultato_tau[best_cr], Eltime = elapsed_time, Consensus_ties = consensus_ranking, Tau_ties = tauX))
  }
  else{
    # Return the results as a list
    return(list(Consensus = consensus_ranking, Tau = tauX, Eltime = elapsed_time))
  }
  
}


####### EMD dataset #####

data("EMD")
load("EMD_results.RData")

fastEMD
quickEMD
decorEMD
psoEMD

consensus_comparison(psoEMD$Consensus[1,], psoEMD$Consensus[2,], X=EMD[,-ncol(EMD)])

# set.seed(1)
# fastEMD = consrank(EMD[,-ncol(EMD)], wk=EMD[,ncol(EMD)], algorithm = "fast")
# 
# set.seed(1)
# quickEMD = consrank(EMD[,-ncol(EMD)], wk=EMD[,ncol(EMD)], algorithm = "quick")
# 
# set.seed(1)
# decorEMD = consrank(EMD[,-ncol(EMD)], wk=EMD[,ncol(EMD)], algorithm = "decor")
# 
# set.seed(1)
# psoEMD = rankagg_pso(EMD[,-ncol(EMD)], 
#             Wk = EMD[,ncol(EMD)],
#             PS = T, 
#             maxit.stagnate = Inf, 
#             maxiter = 1000,
#             phi_p = 1.193147,
#             phi_g = 1.193147,
#             FF = 1.193147/2, 
#             NP = 18,
#             reltol = 0.10)



####### Universities dataset #####
#This code can be very slow!

load("university_data.RData") #university dataset university_soi_big

load("university_results.RData") #results cr_pso_universita_soi_big_new

cr_pso_universita_soi_big_new

View(consensus_comparison(cr_pso_universita_soi_big_new$Consensus[1,], X=university_soi_big))

# cr_pso_universita_soi_big_new = rankagg_pso(university_soi_big,
#                                             PS = T,
#                                             maxit.stagnate = Inf,
#                                             maxiter = 2,
#                                             phi_p = 1.193147,
#                                             phi_g = 1.193147,
#                                             FF = 1.193147/2,
#                                             NP = 18,
#                                             reltol = 0.10)

# dim(cr_pso_universita_soi_big_new$Consensus)



####### Footballers dataset #####
#This code can be slow!

load("serieA_footballers.RData") #calcio, calcio_bkp, calcio_att, calcio_cen, calcio_dif

# full dataset

load("footballers_full_results.RData") #results only

calcio_pso
calcio_decor
calcio_quick
calcio_fast

# set.seed(1)
# calcio_pso = rankagg_pso(calcio,
#                             PS = T, 
#                             maxit.stagnate = Inf, 
#                             maxiter =  1200,
#                             phi_p = 1.193147,
#                             phi_g = 1.193147,
#                             FF = 1.193147/4, #oppure /4
#                             NP = 18,
#                             reltol = 0.0017)
# set.seed(1)
# calcio_decor = consrank(calcio, algorithm = "decor")
# set.seed(1)
# calcio_quick = consrank(calcio, algorithm = "quick") 
# set.seed(1)
# calcio_fast = consrank(calcio, algorithm = "fast") 
# 
# confronto_calcio = (cbind(consensus_comparison(calcio_pso$Consensus[1,], 
#                            calcio_decor$Consensus, X=calcio), 
#       consensus_comparison(calcio_bb, X=calcio)))
# confronto_calcio = data.frame(confronto_calcio)
# colnames(confronto_calcio) = c("PSO", "DECOR", "FindConsensusBB")
# 
# confronto_calcio = cbind(confronto_calcio, consensus_comparison(calcio_quick$Consensus, 
#                                              calcio_fast$Consensus, X=calcio))
# colnames(confronto_calcio) = c("PSO", "DECOR", "FindConsensusBB", "QUICK", "FAST")
# View(confronto_calcio)

#strikers only (i.e. attaccanti) 

load("footballers_strikers_results.RData") #results only

calcio_att_pso
calcio_att_decor
calcio_att_quick
calcio_att_fast

# set.seed(1)
# calcio_att_pso = rankagg_pso(calcio_att,
#                          PS = T, 
#                          maxit.stagnate = Inf, 
#                          maxiter =  1500,
#                          phi_p = 1.193147,
#                          phi_g = 1.193147,
#                          FF = 1.193147/4, #oppure /4
#                          NP = 15,
#                          reltol = 0.0018)
# set.seed(1)
# calcio_att_decor = consrank(calcio_att, algorithm = "decor")
# 
# set.seed(1)
# calcio_att_quick = consrank(calcio_att, algorithm = "quick")
# 
# set.seed(1)
# calcio_att_fast = consrank(calcio_att, algorithm = "fast")
# 
# confronto_calcio_att = (cbind(consensus_comparison(calcio_att_pso$Consensus, 
#                                                    calcio_att_decor$Consensus, X=calcio_att), 
#                           consensus_comparison(calcio_att_quick$Consensus, X=calcio_att)))
# 
# confronto_calcio_att = data.frame(confronto_calcio_att)
# colnames(confronto_calcio_att) = c("PSO=0.3902", "DECOR=0.3855", "QUICK=0.3930")
# View(confronto_calcio_att)

#defenders only (i.e. difensori) 
load("footballers_defenders_results.RData") #results only

calcio_dif_pso
calcio_dif_decor
calcio_dif_quick
calcio_dif_fast

# set.seed(1)
# calcio_dif_pso = rankagg_pso(calcio_dif,
#                              PS = T, 
#                              maxit.stagnate = Inf, 
#                              maxiter =  1200,
#                              phi_p = 1.193147,
#                              phi_g = 1.193147,
#                              FF = 1.193147/4,
#                              NP = 18,
#                              reltol = 0.0017)
# set.seed(1)
# calcio_dif_decor = consrank(calcio_dif, algorithm = "decor")
# 
# set.seed(1)
# calcio_dif_quick = consrank(calcio_dif, algorithm = "quick")
# 
# set.seed(1)
# calcio_dif_fast = consrank(calcio_dif, algorithm = "fast")

#midfielders only (i.e. centrocampisti) 

load("footballers_midfielders_results.RData") #results only

calcio_cen_pso
calcio_cen_decor
calcio_cen_quick
calcio_cen_fast

# set.seed(1)
# calcio_cen_pso = rankagg_pso(calcio_cen,
#                              PS = T, 
#                              maxit.stagnate = Inf, 
#                              maxiter =  1200,
#                              phi_p = 1.193147,
#                              phi_g = 1.193147,
#                              FF = 1.193147/4, #oppure /4
#                              NP = 18,
#                              reltol = 0.0017)
# set.seed(1)
# calcio_cen_decor = consrank(calcio_cen, algorithm = "decor")
# 
# set.seed(1)
# calcio_cen_quick = consrank(calcio_cen, algorithm = "quick")
# 
# set.seed(1)
# calcio_cen_fast = consrank(calcio_cen, algorithm = "fast")


##### Simulated data ####


load("simulated_data_results.RData") #simulazione
View(simulazione)


# this full code can be very, very slow!
# library(PerMallows)
# soluzione_sim = c(1:100) #this is the central tendency of the distribution
# set.seed(1)
# data_sim = rgmm(20, soluzione_sim, rep(0.0125, length(soluzione_sim)-1))
# round(colMeans(data_sim))
# 
# set.seed(1)
# decor_sim = consrank(data_sim, algorithm = "decor", full = T)
# set.seed(1)
# quick_sim = consrank(data_sim, algorithm = "quick", full = T)
# set.seed(1)
# fast_sim = consrank(data_sim, algorithm = "fast", full = T)
# 
# set.seed(1)
# pso_sim = rankagg_pso(data_sim,
#             PS = T, 
#             maxit.stagnate = Inf, 
#             maxiter = 200,
#             phi_p = 1.193147,
#             phi_g = 1.193147,
#             FF = 1.193147/2, 
#             NP = 18,
#             reltol = 0.0017)
# 
# # this process can be very, very slow!
# library(RankAggSIgFUR)
# 
# start = Sys.time()
# sigfur_sim = sigfur(t(data_sim), 2:3, 10, 2:3, 30)
# sigfur_tempo = Sys.time()-start
# 
# #composition of the experiment parameters
# theta_values = c(0.0125, 0.1, 0.8)
# items_values = c(100, 500, 1000)
# algorithm = c("PSO", "DECOR", "QUICK", "FAST", "SIGFUR")
# simulazione = data.frame(algorithm = rep(algorithm, each = length(theta_values)*length(items_values)),
#            theta = rep(rep(theta_values, each = length(items_values)), length(algorithm)),
#            items_values = rep(items_values, length(algorithm)*length(theta_values)),
#            tau = NA,
#            time_s = NA
#            )
# #removal of some very slow combinations
# simulazione = simulazione[!(simulazione$items_values == 1000 & simulazione$algorithm != "PSO"),]
# rownames(simulazione) = 1:nrow(simulazione)
# 
# simulazione$items_values[simulazione$items_values == 500 & simulazione$algorithm == "SIGFUR"] = 250
# 
# 
# #addin
# simulazione = rbind(simulazione, data.frame(algorithm = rep(algorithm[2:4], each = length(theta_values)),
#                          theta = rep(rep(theta_values, each = 1), length(algorithm[2:4])),
#                          items_values = rep(250, length(algorithm[2:4])*length(theta_values)),
#                          tau = NA,
#                          time_s = NA
# ))
# 
# simulazione = rbind(simulazione, data.frame(algorithm = rep("PSO", each = length(theta_values)),
#                               theta = rep(rep(theta_values, each = 1), 1),
#                               items_values = rep(250, 1*length(theta_values)),
#                               tau = NA,
#                               time_s = NA
# ))
# 
# for(i in 34:nrow(simulazione)){
#   print(round(i/nrow(simulazione)*100))
#   soluzione_sim = 1:simulazione$items_values[i]
#   set.seed(1)
#   data_sim = rgmm(20, soluzione_sim, rep(simulazione$theta[i], length(soluzione_sim)-1))
#   
#   print(paste("alg =", simulazione$algorithm[i], "theta =", simulazione$theta[i], "n =", simulazione$items_values[i]))
#   set.seed(1)
#   if(!(simulazione$algorithm[i] %in% c("QUICK", "FAST") & simulazione$items_values[i] == 500)){
#     if(simulazione$algorithm[i] == "PSO"){
#       pso_sim = rankagg_pso(data_sim,
#                             PS = F, 
#                             maxit.stagnate = Inf, 
#                             maxiter = 200,
#                             phi_p = 1.193147,
#                             phi_g = 1.193147,
#                             FF = 1.193147/2, 
#                             NP = 18,
#                             reltol = 0.0017)
#       simulazione$tau[i] = pso_sim$Tau[1]
#       simulazione$time_s[i] = as.numeric(pso_sim$Eltime, units="secs")
#     }
#     if(simulazione$algorithm[i] == "DECOR"){
#       decor_sim = consrank(data_sim, algorithm = "decor", full = T, ps=F)
#       if(length(decor_sim$Tau) == 1){
#         simulazione$tau[i] = decor_sim$Tau[1]
#       }
#       else{
#         simulazione$tau[i] = decor_sim$Tau[1,1]
#       }
#       simulazione$time_s[i] = unname(decor_sim$Eltime)
#     }
#     if(simulazione$algorithm[i] == "QUICK"){
#       quick_sim = consrank(data_sim, algorithm = "quick", full = T, ps=F)
#       if(length(quick_sim$Tau) == 1){
#         simulazione$tau[i] = quick_sim$Tau[1]
#       }
#       else{
#         simulazione$tau[i] = quick_sim$Tau[1,1]
#       }
#       simulazione$time_s[i] = unname(quick_sim$Eltime)
#     }
#     if(simulazione$algorithm[i] == "FAST"){
#       fast_sim = consrank(data_sim, algorithm = "fast", full = T, ps=F)
#       if(length(fast_sim$Tau) == 1){
#         simulazione$tau[i] = fast_sim$Tau[1]
#       }
#       else{
#         simulazione$tau[i] = fast_sim$Tau[1,1]
#       }
#       simulazione$time_s[i] = unname(fast_sim$Eltime)
#     }
#     if(simulazione$algorithm[i] == "SIGFUR"){
#       start = Sys.time()
#       sigfur_sim = sigfur(t(data_sim), 2:3, 10, 2:3, 30)
#       sigfur_tempo = Sys.time()-start
#       
#       simulazione$tau[i] = sigfur_sim$tau
#       simulazione$time_s[i] = as.numeric(sigfur_tempo, units="secs")
#     }
#     print(paste("tau =", simulazione$tau[i], "time_s =", simulazione$time_s[i]))
#   }
#   else{
#     print("SKIP")
#   }
# }
# 
# bkp_simulazione = simulazione
# 
# 
# for(j in 1:nrow(consensus_tied_list)){
#   cr = consensus_tied_list[j,]
#   cr = rank2order(cr)
#   if(j == 1){
#     risultato = recursive_split(cr = cr, max_solutions = 100)
#   }
#   else{
#     risultato = rbind(risultato, recursive_split(cr = cr, max_solutions = 100))
#   }
# 
# }
# 
# risultato = risultato[!duplicated(risultato),]
# 
# for(j in 1:nrow(risultato)){
#   print(j/nrow(risultato))
#   cr = order2rank(risultato[j,])
#   cr = cr[,order(cr)]
#   if(j == 1){
#     risultato_tau = mean(tau_x(data_sim, cr))
#   }
#   else{
#     risultato_tau = c(risultato_tau, mean(tau_x(data_sim, cr)))
#   }
# }
# #1024 equivalent solutions without ties with the same tau = 0.9868819
# 
# which.max(risultato_tau)
# 
# 
# # addin simulation of PSO with FULL=T
# for(i in 1:nrow(simulazione)){
#   print(round(i/nrow(simulazione)*100))
#   soluzione_sim = 1:simulazione$items_values[i]
#   set.seed(1)
#   data_sim = rgmm(20, soluzione_sim, rep(simulazione$theta[i], length(soluzione_sim)-1))
#   
#   print(paste("alg =", simulazione$algorithm[i], "theta =", simulazione$theta[i], "n =", simulazione$items_values[i]))
#   set.seed(1)
#   if(simulazione$algorithm[i] == "PSO"){
#       pso_sim = rankagg_pso(data_sim,
#                             PS = F, 
#                             maxit.stagnate = Inf, 
#                             maxiter = 200,
#                             phi_p = 1.193147,
#                             phi_g = 1.193147,
#                             FF = 1.193147/2, 
#                             NP = 18,
#                             reltol = 0.0017,
#                             FULL = T)
#       simulazione$tau[i] = pso_sim$Tau[1]
#       simulazione$time_s[i] = as.numeric(pso_sim$Eltime, units="secs")
#       print(paste("tau =", simulazione$tau[i], "tau_ties =", pso_sim$Tau_ties[1], "time_s =", simulazione$time_s[i], "Equivalent consensus founded =", nrow(pso_sim$Consensus)))
#   }
#   else{
#     print("SKIP")
#   }
# }
# 
# # addin simulation with QUICK N=500
# for(i in 19:nrow(simulazione)){
#   print(round(i/nrow(simulazione)*100))
#   soluzione_sim = 1:simulazione$items_values[i]
#   set.seed(1)
#   data_sim = rgmm(20, soluzione_sim, rep(simulazione$theta[i], length(soluzione_sim)-1))
#   
#   print(paste("alg =", simulazione$algorithm[i], "theta =", simulazione$theta[i], "n =", simulazione$items_values[i]))
#   set.seed(1)
#   if(simulazione$algorithm[i] == "QUICK" & simulazione$items_values[i] == 500){
#     quick_sim = consrank(data_sim, algorithm = "quick", full = T, ps=F)
#     if(length(quick_sim$Tau) == 1){
#       simulazione$tau[i] = quick_sim$Tau[1]
#     }
#     else{
#       simulazione$tau[i] = quick_sim$Tau[1,1]
#     }
#     simulazione$time_s[i] = unname(quick_sim$Eltime)
#     print(paste("tau =", simulazione$tau[i], "time_s =", simulazione$time_s[i]))
#   }
#   else{
#     print("SKIP")
#   }
# }
# 
# 




