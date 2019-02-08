##################################################################
# Title: Normalizing CyTOF data to machine controls
# The purpose of this script is to normalize the signal intensities of several CyTOF runs across each other so that samples can be 
# compared between runs
##################################################################

##### Experiment design: 
#1. Samples were acquired over 4 days, from 3 separate staining experiments, labeled MM183, MM184, MM185-1, MM185-2
#2. The same sample (MC0067) was thawed, stained, and run as a machine control for each run
#3. For all experiments, Live Singlet CD45+ CD3- CD19- CD14- CD34- CD56+ NK cells were manually gated in Cytobank
#4. NK cell populations from each of the machine controls were exported from CyTOF as csv files with signal intensities per cell


##### Set working directory and load libraries

setwd("/Users/mmack/Box/Documents/Experiments/MM-185 Pre-Post dupi NK CyTOF/Combining and normalizing")
library(Biobase)
library(ggplot2)
library(RColorBrewer)

##################################################################
##### Step 1: Import machine control files, extract medians #####
##################################################################

# import files to run day datasets
mm183.mc <- read.csv("MM183_MC0067.fcs_spill_applied_CD56+ CD3-_MC0067.fcs.csv",header=T, row.names=1)
mm184.mc <- read.csv("MM184_MC0067.fcs_spill_applied_CD56+ CD3-_MC0067.fcs.csv",header=T, row.names=1)
mm185.mc1 <- read.csv("MM185_d1_MC0067.fcs_spill_applied_CD56+ CD3-_MC0067.fcs.csv",header=T,row.names=1)
mm185.mc2 <- read.csv("MM185_d2_MC0067.fcs_spill_applied_CD56+ CD3-_MC0067.fcs.csv",header=T,row.names=1)

# remove event length and time variables
drops <- c("Time", "Event_length")
mm183.mc <- mm183.mc[!(colnames(mm183.mc) %in% drops)]
mm184.mc <- mm184.mc[!(colnames(mm184.mc) %in% drops)]
mm185.mc1 <- mm185.mc1[!(colnames(mm185.mc1) %in% drops)]
mm185.mc2 <- mm185.mc2[!(colnames(mm185.mc2) %in% drops)]

# transform so that the individual events are on the columns and channels are the rows
mm183.mc <- t(mm183.mc)
mm184.mc <- t(mm184.mc)
mm185.mc1 <- t(mm185.mc1)
mm185.mc2 <- t(mm185.mc2)

# create a row medians column for each channel
mm183.mc <- data.frame(mm183.mc, medians = rowMedians(mm183.mc, na.rm = TRUE))
mm184.mc <- data.frame(mm184.mc, medians = rowMedians(mm184.mc, na.rm = TRUE))
mm185.mc1 <- data.frame(mm185.mc1, medians = rowMedians(mm185.mc1, na.rm = TRUE))
mm185.mc2 <- data.frame(mm185.mc2, medians = rowMedians(mm185.mc2, na.rm = TRUE))

# check that this worked out okay [should end in a 'medians' column]
tail(colnames(mm183.mc))

##################################################################
##### Step 2: Create 'delta' vectors for each MC compared to the 
#####         MM183 MC (chosen standard) 
##################################################################

# I have chosen a ratio approach (ie divide by the MM183 median for each channel)
#eg., if mm183 channel X is 504 and mm184 channel X is 312, then the norm vector value for channel X would be 312/504 = 0.6190476
#Then all the channel X values for each event on mm184 run day would be divided by 0.6190476 when normalized

# Covert all 0's to 1.0 in order to avoid diving by zero -- this way all channels with 0 values will just not be normalized
mm183.mc[mm183.mc == 0] <- 1.0
mm184.mc[mm184.mc == 0] <- 1.0
mm185.mc1[mm185.mc1 == 0] <- 1.0
mm185.mc2[mm185.mc2 == 0] <- 1.0

# 'nv' for normalization vector
mm183.nv <- rep(1, nrow(mm183.mc)) # all 1's because it IS the baseline/norm
mm184.nv <- (mm184.mc$medians/mm183.mc$medians)
mm185.nv1 <- (mm185.mc1$medians/mm183.mc$medians)
mm185.nv2 <- (mm185.mc2$medians/mm183.mc$medians)

# divide each event in the MC dataset by the nv
normed_mm183.mc <- mm183.mc/mm183.nv
normed_mm184.mc <- mm184.mc/mm184.nv
normed_mm185.mc1 <- mm185.mc1/mm185.nv1
normed_mm185.mc2 <- mm185.mc2/mm185.nv2

##### generate some plots to show to effect of normalization on medians
# X axis: channels
# Y axis: medians
# line plots for each day's mc and normed_mc

# first make a df of the values I want to plot (medians)
df <- data.frame(channels = rownames(mm183.mc), mm183.mc$medians, mm184.mc$medians, mm185.mc1$medians, mm185.mc2$medians, normed_mm183.mc$medians, normed_mm184.mc$medians, normed_mm185.mc1$medians, normed_mm185.mc2$medians)
row.names(df) <- df$channels
df$channels <- NULL
write.csv(df, "normed_mc_medians_all.csv")

# remove channels that are not associated with markers in this panel, just for plotting
keeps <- c(1, 12, 14:46, 54)
df2 <- df[keeps,]

# rearrange df to make ggplot-ready
channel <- rep(row.names(df2), ncol(df2))
sample <- c(rep("mm183.mc", nrow(df2)), rep("mm184.mc", nrow(df2)), rep("mm185.mc1", nrow(df2)), rep("mm185.d2", nrow(df2)), rep("normed_mm183.mc", nrow(df2)), rep("normed_mm184.mc", nrow(df2)), rep("normed_mm185.mc1", nrow(df2)), rep("normed_mm185.mc2", nrow(df2)))
medians <- c(df2$mm183.mc.medians, df2$mm184.mc.medians, df2$mm185.mc1.medians, df2$mm185.mc2.medians,
             df2$normed_mm183.mc.medians, df2$normed_mm184.mc.medians, df2$normed_mm185.mc1.medians, df2$normed_mm185.mc2.medians)
normed <- c(rep(0, 144), rep(1, 144))
plt <- data.frame(sample = sample, channel = channel, medians = medians, normed = normed)

plt.raw <- plt[normed == 0, ]
plt.norm <- plt[normed == 1, ]

# Plot each median value for each channel as a line plot to visualize the effect of normalization
raw <- ggplot(plt.raw, aes(channel, medians, group=sample, color = sample)) +
  geom_line() +
  xlab('Channels') +
  ylab('Medians') + 
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle=45)) + 
  ggtitle("Raw medians for MC0067 across runs")

norm <- ggplot(plt.norm, aes(channel, medians, group=sample, color = sample)) +
  geom_line() +
  xlab('Channels') +
  ylab('Medians') +
  theme_bw() +
  theme(axis.text.x = element_text(size=8, angle=45)) +
  ggtitle("Normalized medians for MC0067 across runs")

ggsave("Raw_medians_MCs.pdf", raw, device = 'pdf')
ggsave("Normed_medians_MCs.pdf", norm, device = 'pdf')


##################################################################
##### Step 3: Normalize the experimental samples in each run #####
##################################################################


# Divide each sample by the nv for that day's run

mc.normalize <- function(samplefile, nv){
  ## function takes an experimental sample file and normalizes it according to it's normalization vector (Step 2)
  # 'samplefile' = .txt file exported from cytobank, pre-gated on population of choice
  # 'nv' = normalized vector determine from previous step
  drops <- c("Time", "Event_length")
  samplefile <- samplefile[!(colnames(samplefile) %in% drops)]
  samplefile <- t(samplefile)
  normed_sample <- samplefile/nv
  normed_sample <- t(normed_sample)
  return(normed_sample)
}

# usage: 
# samplefile <- read.csv()
# normed_sample <- mc.normalize(samplefile, nv)

# test that this works with my mc samples fom earlier
mm185.mc2 <- read.csv("MM185_d2_MC0067.fcs_spill_applied_CD56+ CD3-_MC0067.fcs.csv",header=T,row.names=1)
normed2_mm185.mc2 <- mc.normalize(mm185.mc2, mm185.nv2)
write.csv(normed2_mm185.mc2, file = paste0("./Exp files/normed_files/test.csv"), row.names=FALSE)

# Note that this function does not remove any 'unused' channels, it does this to all channels collected on the instrument except for time and event_length

##### Load all the files for this project
# read in each day's files separately as a list so that the same norm vector can be applied to the appropriate datasets using lapply

setwd("/Users/mmack/Box/Documents/Experiments/MM-185 Pre-Post dupi NK CyTOF/Combining and normalizing/Exp files/mm183")
mm183_files <- list.files(".")
mm183_data <- lapply(mm183_files, read.delim, row.names=1)

setwd("/Users/mmack/Box/Documents/Experiments/MM-185 Pre-Post dupi NK CyTOF/Combining and normalizing/Exp files/mm184")
mm184_files <- list.files(path=".")
mm184_data <- lapply(mm184_files, read.delim, row.names=1)

setwd("/Users/mmack/Box/Documents/Experiments/MM-185 Pre-Post dupi NK CyTOF/Combining and normalizing/Exp files/mm185_1")
mm185.1_files <- list.files(path=".")
mm185.1_data <- lapply(mm185.1_files, read.delim, row.names=1)

setwd("/Users/mmack/Box/Documents/Experiments/MM-185 Pre-Post dupi NK CyTOF/Combining and normalizing/Exp files/mm185_2")
mm185.2_files <- list.files(path=".")
mm185.2_data <- lapply(mm185.2_files, read.delim, row.names=1)

##### Now apply the mc.normalize function to each set of datafiles
# Using another function that applies 'mc.normalize' and then writes out each normed data.frame to a .csv file

normalize_and_write <- function(datalist, nv, datafiles){
  # datalist is a list of data frames containing the raw data for each sample
  # nv the appropriate normalization vector for this day's mc
  # list of file names corresponding to the datalist
  
  # step 1 - normalize the data 
  normed_datalist <- lapply(datalist, mc.normalize, nv)
  
  
  # step 2 - name the elements of the list so we can use it to name the files
  names(normed_datalist) <- paste0("normed_",datafiles)
  
  
  # step 3 - write out as individual files to the normed_files folder
  for(i in (1:length(normed_datalist))){
    write.csv(normed_datalist[[i]], file = paste0("./Exp files/normed_files/",names(normed_datalist[i]), ".csv"))
  }
}

normalize_and_write(mm183_data, mm183.nv, mm183_files)

normalize_and_write(mm184_data, mm184.nv, mm184_files)

normalize_and_write(mm185.1_data, mm185.nv1, mm185.1_files)

normalize_and_write(mm185.2_data, mm185.nv2, mm185.2_files)


# Okay! Now these .csv files can be imported into Cytobank and converted back to fcs files using their software. 


sessionInfo()


