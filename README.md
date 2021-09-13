###Extended Blough Model

This code is to simulate training and test data as predicted by our Extended Blough Model. There are several parameters that can be changed by users to predict generalization for a variety of experimental designs. This code was designed to simulate training and test for 2 groups, but can be used to simulate a single group. As an example, the initial values have been inputed to simulate Experiment 1A. 

Details for each parameter are as follows:

#Inputs 
beta_plus <- Learning rate for reinforced trials
beta_minus <- Learning rate for non-reinforced trials
lamb_plus <- Learning assymptote for reinforced trials
lamb_minus <- Learning assymptote for non-reinforced trials
sd <- Standard deviation of Gaussian function for elemental activation

### Experimental design
name <- name of experiment/design for filenames when saving/plotting
stim_num_d1 <- sets size of dimension 1 (d1) in stimulus units 
stim_inc_d1 <- sets how many elemental units between stimuli on d1. Increasing elements smooth out plots, but will take longer to run simulations
stim_num_d2 <- sets size of dimension 2 (d2) in stimulus units 
stim_inc_d2 <- sets how many elemental units between stimuli on d2
dim_lims <- set boundaries of both stimulus dimensions. Increase limits to prevent activation functions cutting off for stimuli at ends of stimulus dimensions

g1_name <- name of group 1 (g1) for filenames when saving/plotting
g1_training_stimuli_names <- vector of names of g1 stimuli for saving/plotting
g1_training_trialtypes <- vector of numeric values representing trial types for g1 training stimuli. 1 is reinforcement and 0 is non-reinforcement. 
g1_d1_training_stimuli_elements <- vector of numeric values representing each training stimulus' most activated element on d1
g1_d2_training_stimuli_elements <- vector of numeric values representing each training stimulus' most activated element on d2
g1_N <- sets how many participants to simulate for g1. As this code randomizes training trials, there can be variation between participants that may be smoothed out by repeated simulations. Larger N's may take longer to simulate however. Note: users can set to 0 if they only wish to simulate 1 group. 

g2_name <- name of group 2 (g2) for filenames when saving/plotting
g2_training_stimuli_names <- vector of names of g2 stimuli for saving/plotting
g2_training_trialtypes <- vector of numeric values representing trial types for g1 training stimuli.
g2_d1_training_stimuli_elements <- vector of numeric values representing each training stimulus' most activated element on d1
g2_d2_training_stimuli_elements <- vector of numeric values representing each training stimulus' most activated element on d2
g2_N <- sets how many participants to simulate for g2

#test stimuli
test_stimuli_names <- vector of names of test stimuli for saving/plotting
d1_test_stimuli_elements <- vector of numeric values representing each test stimulus' most activated element on d1
d2_test_stimuli_elements <- vector of numeric values representing each test stimulus' most activated element on d1

#experiment design
stimreps_per_block <- sets how many times a training stimulus is repeated per training trial block
num_blocks <- sets how many training blocks to simulate

#saving data
savetraining <-  set to TRUE to save a .txt file of training data
savetest <-  set to TRUE to save .txt file of test data

#misc.
digits <-  rounding for elemental activation function. Decrease digits to limit the spreading activation of each stimulus, or use if elemental activation functions substantially cut off for stimuli at ends of dimensions


### inputs for other simulations
# Experiment 1B (partial reinforcement)
We simulated a partial reinforcement design, by repeating the CS+ and CS-s in our vectors of training stimuli and setting the trial type to 0 for one of the repeated CS+s. 

# Experimental design
name <- 'Exp1A' 
stim_num_d1 <- 11 
stim_inc_d1 <- 0.25
stim_num_d2 <- 11 
stim_inc_d2 <- 0.25 
dim_lims <- 5

g1_name <- "Single Positive"
g1_training_stimuli_names <- c("CS+","CS+","CS+","CS+)
g1_training_trialtypes <- c(1,1,1,0) 
g1_d1_training_stimuli_elements <- c(6,6,6,6)
g1_d2_training_stimuli_elements <- c(1,1,1,1)
g1_N <- 56 

g2_name <- "Double Negative"
g2_training_stimuli_names <- c("CS1-","CS1-","CS1-","CS1-","CS+","CS+","CS+","CS+","CS2-","CS2-","CS2-","CS2-")
g2_training_trialtypes <- c(0,0,0,0,1,1,1,0,0,0,0,0) 
g2_d1_training_stimuli_elements <- c(4,4,4,4,6,6,6,6,8,8,8,8)
g2_d2_training_stimuli_elements <- c(1,1,1,1,1,1,1,1,1,1,1,1)
g2_N <- 46 

#test stimuli
test_stimuli_names <- c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11")
d1_test_stimuli_elements <- c(1:11)
d2_test_stimuli_elements <- c(rep(1,11))

#experiment design
stimreps_per_block <- 1
num_blocks <- 3

#saving data
savetraining <-  TRUE
savetest <-  TRUE

#Misc.
digits <-  20

#Experiment 2 (shape-change)
These inputs are essentially the same as Experiment 1A's, except we simulated the shape-change at test by inputting different elements on d2 between training and test stimuli. The elemental distance between the training and test shape was equal to the elemental distance between GS1 and GS11 on the colour dimension 

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
g2_N <- 48

#test stimuli
test_stimuli_names <- c("S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11")
d1_test_stimuli_elements <- c(1:11)
d2_test_stimuli_elements <- c(rep(11,11)) 

#experiment design
stimreps_per_block <- 2 
num_blocks <- 6
