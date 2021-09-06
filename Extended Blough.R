rm(list=ls())
options(scipen=999)
library(dplyr)
library(ggplot2)

#### define parameters ####
#Inputs 
beta_plus <- 0.2
beta_minus <- 0.25
lamb_plus <- 1
lamb_minus <- 0
sd <- 10 #activation function sigma

# Experimental design
name <- 'Jess exemplars_sd20' #name (for plotting)
stim_num_d1 <- 100 # how many stimuli on dimension 1? -> sets dimension size/stimuli distance
stim_inc_d1 <- 0.25 # how many elements between stimuli? -> sets smoothness of plots
stim_num_d2 <- 100 # how many stimuli on dimension 2?
stim_inc_d2 <- 0.25 # how many elements between stimuli?
dim_lims <- 5 #set boundaries of stimulus dimension, increase to prevent activation functions cutting off for stimuli at extremes

g1_name <- "positive only"
g1_training_stimuli_names <- c("CS1+","CS2+","CS3+","CS4+")
g1_training_trialtypes <- c(1,1,1,1) # 1 = reinforce, 0 = non-reinforce
g1_d1_training_stimuli_elements <- c(50,75,50,25)
g1_d2_training_stimuli_elements <- c(75,50,25,50)
g1_N <- 0 #how many participants?

g2_name <- "positive-negative"
g2_training_stimuli_names <- c("CS1+","CS2+","CS3+","CS4+","CS1-","CS2-","CS3-","CS4-")
g2_training_trialtypes <- c(1,1,1,1,0,0,0,0) 
g2_d1_training_stimuli_elements <- c(50,75,50,25,50,85,50,15)
g2_d2_training_stimuli_elements <- c(75,50,25,50,85,50,15,50)
g2_N <- 1

#test stimuli
test_stimuli_names <- c("CS1+","CS2+","CS3+","CS4+",
                        "CS1-","CS2-","CS3-","CS4-",
                        "T1","T2","T3","T4","T5","T6","T7","T8")
d1_test_stimuli_elements <- c(50,75,50,25,
                              50,85,50,15,
                              50,60,50,40,50,95,50,5)
d2_test_stimuli_elements <- c(75,50,25,50,
                              85,50,15,50,
                              60,50,40,50,95,50,5,50)

stimreps_per_block <- 2
num_blocks <- 6
partial_reinf <- FALSE # if TRUE, randomly turns CS+ trials into CS- 
partial_reinf_rate <- 0.75

#Misc.
digits = 20 # rounding for elemental activation function, change to set how much spreading activation each stimuli has/use if elemental activation function for stimuli at extremes getting cut off too much


####### MODEL #######
#set up 1st stimulus dimension
d1_ll <- 1 #lower limit
d1_ul <- stim_num_d1 #upper limit

d1_elems <- seq(from = d1_ll-dim_lims,
                to   = d1_ul+dim_lims,
                by   = (stim_inc_d1))

#dimension 2
d2_ll <- 1
d2_ul <- stim_num_d2
d2_elems <- seq(from = d2_ll-dim_lims,
                to   = d2_ul+dim_lims,
                by   = (stim_inc_d2))

# function to calculate activation gradient for stimulus on a dimension
elem_activation <- function(stim,dimension) {
  if (dimension == 1) {
    return(round(dnorm(d1_elems,stim,sd),digits=digits))
  }
  
  if (dimension == 2) {
    return(round(dnorm(d2_elems,stim,sd),digits=digits))
  }
}

#Training stimuli frame
g1_training_stimuli <- data.frame(matrix(c(g1_training_stimuli_names,
                                g1_training_trialtypes,
                                g1_d1_training_stimuli_elements,
                                g1_d2_training_stimuli_elements),
                                nrow=length(g1_training_stimuli_names),
                                ncol=4),
                                stringsAsFactors = FALSE)
colnames(g1_training_stimuli) <- c("stimname","trialtype","d1_e","d2_e")

g1_training_stimuli$trialtype <- as.numeric(g1_training_stimuli$trialtype)
g1_training_stimuli$d1_e <- as.numeric(as.character(g1_training_stimuli$d1_e))
g1_training_stimuli$d2_e <- as.numeric(as.character(g1_training_stimuli$d2_e))

g2_training_stimuli <- data.frame(matrix(c(g2_training_stimuli_names,
                                g2_training_trialtypes,
                                g2_d1_training_stimuli_elements,
                                g2_d2_training_stimuli_elements),
                                nrow=length(g2_training_stimuli_names),
                                ncol=4),
                                stringsAsFactors = FALSE)
colnames(g2_training_stimuli) <- c("stimname","trialtype","d1_e","d2_e")
g2_training_stimuli$trialtype <- as.numeric(g2_training_stimuli$trialtype)
g2_training_stimuli$d1_e <- as.numeric(as.character(g2_training_stimuli$d1_e))
g2_training_stimuli$d2_e <- as.numeric(as.character(g2_training_stimuli$d2_e))

test_stimuli <- data.frame(matrix(c(paste0("S",1:length(d1_test_stimuli_elements)),
                       d1_test_stimuli_elements,
                       d2_test_stimuli_elements),
                       nrow=length(d1_test_stimuli_elements),
                       ncol=3))

colnames(test_stimuli) <- c("stimname","d1_e","d2_e")
test_stimuli$d1_e <- as.numeric(as.character(test_stimuli$d1_e))
test_stimuli$d2_e <- as.numeric(as.character(test_stimuli$d2_e))

##experiment design
trainingdata <- list()
testdata <- data.frame()

##loop training for every participant
for (i in 1:(g1_N+g2_N)) {
  pid <- i
  #set up training schedule
  if (pid>g1_N) {
    group <- 2
    training_stimuli <- g2_training_stimuli
  } else {
    group <- 1
    training_stimuli <- g1_training_stimuli
  }
  
  #set random trial orders + outcomes (1/0 shock/no-shock)
  
  training_trials_order <- sample(c(rep(as.character(training_stimuli$stimname)
                                        ,num_blocks*stimreps_per_block)))
  
  #set partial reinforcement
  if (partial_reinf == TRUE){
    trial_outcome_order[sample(which(trial_outcome_order==1),(1-partial_reinf_rate)*stimreps_per_block*num_blocks,replace=FALSE)] <- 0
    
  }
  
  #set up training data frames
  d1_training <- data.frame(matrix(data = 0, nrow = length(training_trials_order), ncol = length(d1_elems)))
  colnames(d1_training) <- paste0("e",c(1:length(d1_elems)))
  rownames(d1_training) <- paste0("trial",c(1:length(training_trials_order)))
  
  d2_training <- data.frame(matrix(data = 0, nrow = length(training_trials_order), ncol = length(d2_elems)))
  colnames(d2_training) <- paste0("e",c(1:length(d2_elems)))
  rownames(d2_training) <- paste0("trial",c(1:length(training_trials_order)))
  
  ## Simulations
  # function to calculate total associative strength of a stimulus compound on a trial
  stimulus_v <- function(stim_d1,stim_d2,trialnum) { #stim_d1/stim_d2 are the peak elements of each stimulus
    d1_V <- elem_activation(stim_d1,1) * d1_training[trialnum,]
    d2_V <- elem_activation(stim_d2,2) * d2_training[trialnum,]
    return(sum(d1_V) + sum(d2_V))
  }
  
  # training
  #trial 1
  #reinforcement
  if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order[1]] == 1){
    d1_training[1,] <-  elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order[1]],1)*beta_plus*(lamb_plus - 0)
    d2_training[1,] <-  elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order[1]],2)*beta_plus*(lamb_plus - 0)
  }
  
  #non-reinforcement
  if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order[1]] == 0){
    d1_training[1,] <-  elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order[1]],1)*beta_minus*(lamb_minus - 0)
    d2_training[1,] <-  elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order[1]],2)*beta_minus*(lamb_minus - 0)
  }
  
  #rest of the trials
  for (j in 2:length(training_trials_order)) {
    
    #calculate stimulus V for trial
    trial_V <- stimulus_v(
      training_stimuli$d1_e[training_stimuli$stimname == training_trials_order[j]],
      training_stimuli$d2_e[training_stimuli$stimname == training_trials_order[j]],
      j-1)
    
    #if reinforcement trial
    if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order[j]] == 1){
      d1_training[j,] <-  d1_training[j-1,] + elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order[j]],1)*beta_plus*(lamb_plus - trial_V)
      d2_training[j,] <-  d2_training[j-1,] + elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order[j]],2)*beta_plus*(lamb_plus - trial_V)
    }
    
    #if non-reinforcement trial
    if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order[j]] == 0){
      d1_training[j,] <-  d1_training[j-1,] + elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order[j]],1)*beta_minus*(lamb_minus - trial_V)
      d2_training[j,] <-  d2_training[j-1,] + elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order[j]],2)*beta_minus*(lamb_minus - trial_V)
    }
  }
  
  #save training data
  participant_trainingdata <- list(d1_training,d2_training)
  trainingdata[i] <- participant_trainingdata

  # sum associative strengths to simulate test per participant
  # set up data frame
  gen_test <- data.frame(matrix(data = 0, ncol = length(test_stimuli$stimname)))
  colnames(gen_test) <- test_stimuli$stimname
  for (j in test_stimuli$stimname) {
    gen_test[,j] <- stimulus_v(test_stimuli$d1_e[test_stimuli$stimname == j],
                               test_stimuli$d2_e[test_stimuli$stimname == j],
                               length(training_trials_order)) #associative strength of each stimulus at the end of training (the weights on the last trial)
    gen_test$group <- group
    gen_test$pid <- pid
    gen_test$exp <- name
  }
  
  testdata <- rbind(testdata,gen_test)
}


##calculate group means
testdata_means <- data.frame(matrix(data = NA, ncol = length(test_stimuli$stimname)))
colnames(testdata_means) <- test_stimuli_names

for (i in 1:length(unique(testdata$group))){
  tempdata <- data.frame()
  tempdata <- rbind(tempdata,colMeans(testdata[testdata$group==i,1:length(test_stimuli$stimname)])) 
  colnames(tempdata) <- test_stimuli_names
  testdata_means[i,] <- tempdata
}

write.table(testdata,file=paste0("simdata_",name,"_",(Sys.Date()),".txt"),sep="\t",row.names=FALSE)


# plot
jpeg(filename = paste0("Simulation Gen Test","-",name,"-",g1_name,"-",g2_name,"_",Sys.Date(),".jpeg"),
     width = 2200, height = 1500, units = "px", pointsize = 12,
     quality = 1000,res=200,
     bg = "white")

plot(as.numeric(testdata_means[1,]),
     type="l", col="hotpink", lwd=3, ylim=c(-1,1),
     xlab="",ylab="",axes=FALSE,cex.lab=2,ann = FALSE)
lines(as.numeric(testdata_means[2,]),
      type="l", col="steelblue", lwd=3)

par(mar=c(5,5,3,2))

points(as.numeric(testdata_means[1,]),
       pch = 16, col="hotpink",cex=3,lwd=3)
points(as.numeric(testdata_means[2,]),
       pch = 15, col="steelblue",cex=3,lwd=3)

legend("topright",c(g1_name,g2_name), 
       col=c("hotpink","steelblue"),pch=c(16,15), inset=(c(-0.015,-0.03)),pt.cex=3, lwd=3,cex=2, bty="n")

mtext(text = "stimulus",side = 1,line = 3,cex=2)
mtext(text = "associative strength",side = 2,line = 3,cex=2)

axis(side=1,at=0:stim_num_d1,pos=0,
     labels=c(0:stim_num_d1),cex.axis=1.5)

axis(side=2,labels=c("-0.5","0.0","0.5","1.0"),at=c(-0.5,0,0.5,1),cex.axis=2)


dev.off()

#### graph for Ng et al. (2021) ####
# plot activation functions
jpeg(filename = paste0("Simulation activation function","-",name,"-",g1_name,"-",g2_name,"_",Sys.Date(),".jpeg"),
     width = 2200, height = 1300, units = "px", pointsize = 12,
     quality = 1000,res=200,
     bg = "white")

act_functions <- data.frame(CS_plus = elem_activation(6,1),
                            CS1_minus = elem_activation(8,1), 
                            CS2_minus = elem_activation(4,1), 
                            CS_extreme1 = elem_activation(11,1),
                            CS_extreme2 = elem_activation(1,1))

plot(act_functions$CS_plus,type="p",
     col=hsv(.489,1,.75),pch=16,
     cex=1,
     xlab="",ylab="",axes=FALSE,cex.lab=2)
points(act_functions$CS1_minus,type="p",
       col=hsv(.523,1,.75),pch=16,
       cex=1)
points(act_functions$CS2_minus,type="p",
       col=hsv(.447,1,.75),pch=16,
       cex=1)
points(act_functions$CS_extreme1,type="p",
       col=hsv(.583,1,.75),pch=16,
       cex=1)
points(act_functions$CS_extreme2,type="p",
       col=hsv(.396,1,.75),pch=16,
       cex=1)


lines(act_functions$CS_plus,
      col=hsv(.489,1,.75),
      lwd=1.5)
lines(act_functions$CS1_minus,
      col=hsv(.523,1,.75),
      lwd=1.5)
lines(act_functions$CS2_minus,
      col=hsv(.447,1,.75),pch=16,
      cex=1.5)
lines(act_functions$CS_extreme1,
      col=hsv(.583,1,.75),
      lwd=1.5)
lines(act_functions$CS_extreme2,
      col=hsv(.396,1,.75),
      lwd=1.5)

par(mar=c(5,5,3,2))


legend("topleft",c("GS1","CS1-","CS+","CS2-","GS11"),
       col=c(hsv(.396,1,.75),hsv(.447,1,.75),hsv(.489,1,.75),hsv(.523,1,.75),hsv(.583,1,.75)),pch=c(16), inset=(.02),pt.cex=1,lwd=2,cex=2,bty="n")

mtext(text = "element",side = 1,line = 3,cex=2)
mtext(text = "activation",side = 2,line = 3,cex=2)

axis(side=1,pos=0,labels=c("","GS1","CS1-","CS+","CS2-","GS11",""),at=c(0,21,32,41,50,61,81),cex.axis=2)
axis(side=2,cex.axis=2)


dev.off()