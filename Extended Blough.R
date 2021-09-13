rm(list=ls())
options(scipen=999)
library(dplyr)
library(ggplot2)
library(reshape2)

#### define parameters ####
#Inputs 
beta_plus <- 0.2
beta_minus <- 0.25
lamb_plus <- 1
lamb_minus <- 0
sd <- 1

# Experimental design
name <- 'Exp1A' 
stim_num_d1 <- 11 
stim_inc_d1 <- 0.25
stim_num_d2 <- 11 
stim_inc_d2 <- 0.25 
dim_lims <- 5

g1_name <- "Single Positive"
g1_training_stimuli_names <- c("CS+")
g1_training_trialtypes <- c(1) 
g1_d1_training_stimuli_elements <- c(6)
g1_d2_training_stimuli_elements <- c(1)
g1_N <- 46 

g2_name <- "Double Negative"
g2_training_stimuli_names <- c("CS1-","CS+","CS2-")
g2_training_trialtypes <- c(0,1,0) 
g2_d1_training_stimuli_elements <- c(4,6,8)
g2_d2_training_stimuli_elements <- c(1,1,1)
g2_N <- 50 

#test stimuli
test_stimuli_names <- c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11")
d1_test_stimuli_elements <- c(1:11)
d2_test_stimuli_elements <- c(rep(1,11))

#experiment design
stimreps_per_block <- 2 
num_blocks <- 6

#saving data
savetraining <-  TRUE
savetest <-  TRUE

#Misc.
digits <-  20


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

test_stimuli <- data.frame(matrix(c(test_stimuli_names,
                       d1_test_stimuli_elements,
                       d2_test_stimuli_elements),
                       nrow=length(test_stimuli_names),
                       ncol=3))

colnames(test_stimuli) <- c("stimname","d1_e","d2_e")
test_stimuli$d1_e <- as.numeric(as.character(test_stimuli$d1_e))
test_stimuli$d2_e <- as.numeric(as.character(test_stimuli$d2_e))

##experiment design
trainingdata <- data.frame()
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
  training_trials_order <- data.frame(stimname = sample(c(rep(as.character(training_stimuli$stimname)
                                                              ,num_blocks*stimreps_per_block))))
  
  #set up training data frames
  d1_training <- data.frame(matrix(data = 0, nrow = length(training_trials_order$stimname), ncol = length(d1_elems)))
  colnames(d1_training) <- paste0("d1_e",c(1:length(d1_elems)))
  rownames(d1_training) <- paste0("trial",c(1:length(training_trials_order$stimname)))
  
  d2_training <- data.frame(matrix(data = 0, nrow = length(training_trials_order$stimname), ncol = length(d2_elems)))
  colnames(d2_training) <- paste0("d2_e",c(1:length(d2_elems)))
  rownames(d2_training) <- paste0("trial",c(1:length(training_trials_order$stimname)))
  
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
  if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order$stimname[1]] == 1){
    d1_training[1,] <-  elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order$stimname[1]],1)*beta_plus*(lamb_plus - 0)
    d2_training[1,] <-  elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order$stimname[1]],2)*beta_plus*(lamb_plus - 0)
  }
  
  #non-reinforcement
  if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order$stimname[1]] == 0){
    d1_training[1,] <-  elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order$stimname[1]],1)*beta_minus*(lamb_minus - 0)
    d2_training[1,] <-  elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order$stimname[1]],2)*beta_minus*(lamb_minus - 0)
  }
  
  #rest of the trials
  for (j in 2:length(training_trials_order$stimname)) {
    
    #calculate stimulus V for trial
    trial_V <- stimulus_v(
      training_stimuli$d1_e[training_stimuli$stimname == training_trials_order$stimname[j]],
      training_stimuli$d2_e[training_stimuli$stimname == training_trials_order$stimname[j]],
      j-1)
    
    #if reinforcement trial
    if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order$stimname[j]] == 1){
      d1_training[j,] <-  d1_training[j-1,] + elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order$stimname[j]],1)*beta_plus*(lamb_plus - trial_V)
      d2_training[j,] <-  d2_training[j-1,] + elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order$stimname[j]],2)*beta_plus*(lamb_plus - trial_V)
    }
    
    #if non-reinforcement trial
    if (training_stimuli$trialtype[training_stimuli$stimname == training_trials_order$stimname[j]] == 0){
      d1_training[j,] <-  d1_training[j-1,] + elem_activation(training_stimuli$d1_e[training_stimuli$stimname == training_trials_order$stimname[j]],1)*beta_minus*(lamb_minus - trial_V)
      d2_training[j,] <-  d2_training[j-1,] + elem_activation(training_stimuli$d2_e[training_stimuli$stimname == training_trials_order$stimname[j]],2)*beta_minus*(lamb_minus - trial_V)
    }
  }
  
  #save training data
  d1_training$pid <- pid
  d1_training$group <- group
  d1_training$trialnum <- 1:length(d1_training$d1_e1)
  d2_training$pid <- pid
  d2_training$group <- group
  d2_training$trialnum <- 1:length(d2_training$d2_e1)
  participant_trainingdata <- merge(d1_training,d2_training, by = c("pid","group","trialnum"))
  
  if (pid == 1) {
    trainingdata <- participant_trainingdata
  } else {
    trainingdata <- rbind(trainingdata,participant_trainingdata)
  } 

  # sum associative strengths to simulate test per participant
  # set up data frame
  gen_test <- data.frame(matrix(data = 0, ncol = length(test_stimuli$stimname)))
    colnames(gen_test) <- test_stimuli$stimname
  for (j in test_stimuli$stimname) {
    gen_test[,j] <- stimulus_v(test_stimuli$d1_e[test_stimuli$stimname == j],
                               test_stimuli$d2_e[test_stimuli$stimname == j],
                               length(training_trials_order$stimname)) #associative strength of each stimulus at the end of training (the weights on the last trial)
    gen_test$group <- group
    gen_test$pid <- pid
    gen_test$exp <- name
  }
  
  testdata <- rbind(testdata,gen_test)
}


##calculate group means
#change testdata to long format
testdata <- melt(testdata,id = c("group","pid","exp"),variable.name = "stimname",value.name = "response")
testdata_means <- aggregate(testdata$response,list(testdata$group,testdata$stimname),mean)
colnames(testdata_means) <- c("group","stimname","response")
testdata_means$group <- factor(testdata_means$group)

#Save data
setwd("../Simulated Data")
if (savetraining == TRUE) {
  write.table(trainingdata,file=paste0("EB_sim_trainingdata_",name,"_",(Sys.Date()),".txt"),sep="\t",row.names=FALSE)  
}
if (savetest == TRUE) {
  write.table(testdata,file=paste0("EB_sim_testdata_",name,"_",(Sys.Date()),".txt"),sep="\t",row.names=FALSE)  
}

#Plotting
testdata_means
setwd("../Plots")

jpeg(filename = paste0(name,"_",format(Sys.time(), "%d-%b-%Y %H.%M.%S"),".jpeg"),
     width = 4500, height = 1800, units = "px", pointsize = 12,
     quality = 1000,res=200,
     bg = "white")

ggplot(testdata_means, 
       aes(y = response, 
           x = factor(stimname),
           group = group,
           color = group
       )) + 
  geom_line(lwd=2) +  
  geom_point(lwd=7) +
  ylim (-1,1) +
  xlab("stimulus") + ylab("associative strength") +
  theme_gray(base_size = 30) +
  scale_color_manual(values=c('steelblue','hotpink')) +
  scale_x_discrete(labels = c(levels(testdata_means$stimname)))

dev.off()
