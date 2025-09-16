# raw statistical analysis for simple, untransformed, quantitative data (like qPCR) 

library(rstatix)

# remove all the objects from the R session
rm(list=ls())#read input file, containing culture data


input<-fread("Acant_DNA_all_for_r.csv", sep=',', header=TRUE, fill=TRUE)
input.long <- melt(input, measure.vars = c("0","5","11","20","35"))

#############################################################################################################################
############################ Statistical analysis ###############################################################
#############################################################################################################################
###### test for normality ######

SC <- input.long[input.long$Sample == "SC", c(1,2,3)]
inputAc <-  AC[AC$variable == "35",c(1,2,3)]
shapiro.test(inputAc$value) #  normal

AC <- input.long[input.long$Sample == "Ac", c(1,2,3)]
inputAc <-  AC[AC$variable == "35",c(1,2,3)]
shapiro.test(inputAc$value) #  normal

E25 <- input.long[input.long$Sample == "E25", c(1,2,3)]
inputE <-  E25[E25$variable == "35",c(1,2,3)]
shapiro.test(inputE$value) #  normal

UV7 <- input.long[input.long$Sample == "UV7", c(1,2,3)]
inputUV7 <-  UV7[UV7$variable == "20",c(1,2,3)]
shapiro.test(inputUV7$value) #  normal


#########################################
##################################
# non-parametric (Kruskal + Wilcoxon)
##################################

input.long$variable <- factor(input.long$variable, levels =c("0","5","11","20","35"))
input.long <- input.long %>%
  filter(Sample != "SC")
library(dplyr)
library(rstatix)
library(tidyr)
library(purrr)

# Function to process one timepoint with Kruskal-Wallis + pairwise Wilcoxon
process_timepoint_nonparam <- function(input.long, variable) {
  # Kruskal-Wallis test
  kruskal_res <- kruskal.test(value ~ Sample, data = input.long)
  kruskal_p <- kruskal_res$p.value
  
  # Pairwise Wilcoxon test (with BH correction)
  pwilcox <- pairwise.wilcox.test(input.long$value, input.long$Sample, p.adjust.method = "BH")
  
  # Convert matrix of results into tidy dataframe
  pw_df <- as.data.frame(as.table(pwilcox$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(Sample1 = Var1,
           Sample2 = Var2,
           Wilcox_pvalue = Freq) %>%
    mutate(
      Timepoint = variable,
      Kruskal_pvalue = kruskal_p,
      Wilcox_sig = case_when(
      	Wilcox_pvalue < 0.0001 ~ "****",
        Wilcox_pvalue < 0.001 ~ "***",
        Wilcox_pvalue < 0.01 ~ "**",
        Wilcox_pvalue < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    select(Timepoint, Sample1, Sample2, Kruskal_pvalue, Wilcox_pvalue, Wilcox_sig)
  
  return(pw_df)
}

# Apply to each timepoint
results_samplediff_nonparam <- input.long %>%
  group_by(variable) %>%
  group_split() %>%
  map_df(~process_timepoint_nonparam(.x, unique(.x$variable)))

write.csv(results_samplediff_nonparam, "Stats_wilcox_results_sample_differences_all.csv", row.names = FALSE)
#############################################################################################################################

# Function to process one timepoint with Kruskal-Wallis + pairwise Wilcoxon
process_sample_nonparam <- function(input.long, Sample) {
  # Kruskal-Wallis test
  kruskal_res <- kruskal.test(value ~ variable, data = input.long)
  kruskal_p <- kruskal_res$p.value
  
  # Pairwise Wilcoxon test (with BH correction)
  pwilcox <- pairwise.wilcox.test(input.long$value, input.long$variable, p.adjust.method = "BH")
  
  # Convert matrix of results into tidy dataframe
  pw_df <- as.data.frame(as.table(pwilcox$p.value)) %>%
    filter(!is.na(Freq)) %>%
    rename(Time1 = Var1,
           Time2 = Var2,
           Wilcox_pvalue = Freq) %>%
    mutate(
      Sample = Sample,
      Kruskal_pvalue = kruskal_p,
      Wilcox_sig = case_when(
      	Wilcox_pvalue < 0.0001 ~ "****",
        Wilcox_pvalue < 0.001 ~ "***",
        Wilcox_pvalue < 0.01 ~ "**",
        Wilcox_pvalue < 0.05 ~ "*",
        TRUE ~ "ns"
      )
    ) %>%
    select(Sample, Time1, Time2, Kruskal_pvalue, Wilcox_pvalue, Wilcox_sig)
  
  return(pw_df)
}

# Apply to each timepoint
results_time_nonparam <- input.long %>%
  group_by(Sample) %>%
  group_split() %>%
  map_df(~process_sample_nonparam(.x, unique(.x$Sample)))



write.csv(results_time_nonparam, "wilcox_Stats_results_change_over_time.csv", row.names = FALSE)
#############################################################################################################################
##################################
# parametric (ANOVA + tukey hsd)
##################################
# Function to process one timepoint
process_timepoint <- function(input.long, variable) {
  # Run ANOVA
  stat <- aov(value ~ Sample, data = input.long)
  anova_p <- summary(stat)[[1]][["Pr(>F)"]][1]
  
  # Tukey HSD with rstatix
  tukey <- tukey_hsd(stat)
  
  # Add timepoint & ANOVA p-value
  tukey %>%
    select(group1, group2, p.adj, p.adj.signif) %>%
    mutate(Timepoint = variable,
           Anova_pvalue = anova_p) %>%
    rename(Sample1 = group1,
           Sample2 = group2,
           TukeyHSD_pvalue = p.adj,
           TukeyHSD_sig = p.adj.signif) %>%
    select(Timepoint, Sample1, Sample2, Anova_pvalue, TukeyHSD_pvalue, TukeyHSD_sig)
}

# Apply for each timepoint
results_samplediff <- input.long %>%
  group_by(variable) %>%
  group_split() %>%
  map_df(~process_timepoint(.x, unique(.x$variable)))

write.csv(results_samplediff, "Stats_results_sample_differences_all.csv", row.names = FALSE)
#############################################################################################################################

process_sample <- function(input.long, Sample) {
  # ANOVA across timepoints
  stat <- aov(value ~ variable, data = input.long)
  anova_p <- summary(stat)[[1]][["Pr(>F)"]][1]
  
  # Tukey HSD
  tukey <- tukey_hsd(stat)
  
  # Build results
  tukey %>%
    select(group1, group2, p.adj, p.adj.signif) %>%
    mutate(Sample = Sample,
           Anova_pvalue = anova_p) %>%
    rename(Time1 = group1,
           Time2 = group2,
           TukeyHSD_pvalue = p.adj,
           TukeyHSD_sig = p.adj.signif) %>%
    select(Sample, Time1, Time2, Anova_pvalue, TukeyHSD_pvalue, TukeyHSD_sig)
}

results_time <- input.long %>%
  group_by(Sample) %>%
  group_split() %>%
  map_df(~process_sample(.x, unique(.x$Sample)))

write.csv(results_time, "Stats_results_change_over_time.csv", row.names = FALSE)
#############################################################################################################################






########################################################################################################################
###### Statistical analysis by hand #################################################################################
########################################################################################################################
###### 
input0 <-  input.long[input.long$variable == "T0",c(1,2,3)]
input5 <-  input.long[input.long$variable == "T5",c(1,2,3)]
input11 <-  input.long[input.long$variable == "T11",c(1,2,3)]
input20 <-  input.long[input.long$variable == "T20",c(1,2,3)]
input35 <-  input.long[input.long$variable == "35",c(1,2,3)]

setkey(input5, Sample)
input5 <- input5[!"SC"]
setkey(input11, Sample)
input11 <- input11[!"SC"]
setkey(input20, Sample)
input20 <- input20[!"SC"]


SC <- input.long[input.long$Sample == "SC", c(1,2,3)]
AC <- input.long[input.long$Sample == "Ac", c(1,2,3)]
E25 <- input.long[input.long$Sample == "E25", c(1,2,3)]
UV7 <- input.long[input.long$Sample == "UV7", c(1,2,3)]

######## Compare means of multiple groups using ANOVA ###################################################################
stat0 <- aov(value ~ Sample, data = input0)
summary(stat0)
TukeyHSD(stat0)
rstatix::tukey_hsd(stat0)

stat5 <- aov(value ~ Sample, data = input5)
summary(stat5)
TukeyHSD(stat5)
rstatix::tukey_hsd(stat5)

stat11 <- aov(value ~ Sample, data = input11)
summary(stat11)
TukeyHSD(stat11)
rstatix::tukey_hsd(stat11)


stat20 <- aov(value ~ Sample, data = input20)
summary(stat20)
TukeyHSD(stat20)
rstatix::tukey_hsd(stat20)

stat35 <- aov(value ~ Sample, data = input35)
summary(stat35)
TukeyHSD(stat35)
rstatix::tukey_hsd(stat35)


statAc <- aov(value ~ variable, data = AC)
summary(statAc)
TukeyHSD(statAc)
rstatix::tukey_hsd(statAc)

statE <- aov(value ~ variable, data = E25)
summary(statE)
TukeyHSD(statE)
rstatix::tukey_hsd(statE)

statU <- aov(value ~ variable, data = UV7)
summary(statU)
TukeyHSD(statU)
rstatix::tukey_hsd(statU)


######## Test for normaility of the data ################################################################################

shapiro.test(input0$value) # not normal
shapiro.test(input5$value) # not normal
shapiro.test(input11$value) # not normal
shapiro.test(input20$value) # not normal
shapiro.test(input35$value) # not normal

shapiro.test(AC$value) # not normal
shapiro.test(E25$value) # not normal
shapiro.test(UV7$value) # not normal


shapiro.test(input.long$value) # not normal

######## Compare means of multiple groups using Kruskal-wallis (for non-normal data) #################################

sum(rank(input20$value)[input20$Sample=="UV7"]) - ((7*8)/2)
sum(rank(input0$value)[input0$Sample=="Ac"])

kruskal.test(value ~ Sample, data = input0)
pairwise.wilcox.test(input0$value, input0$Sample, p.adjust.method = "BH", conf.int = TRUE)

setkey(input5, Sample)
input5_b <- input5[!"SC"]
input5_b <- input5[!"Ac"]

wilcox.test(value ~ Sample, data = input5_b, conf.int = TRUE)

kruskal.test(value ~ Sample, data = input5)
pairwise.wilcox.test(input5$value, input5$Sample, p.adjust.method = "BH")

kruskal.test(value ~ Sample, data = input11)
pairwise.wilcox.test(input11$value, input11$Sample, p.adjust.method = "BH")

kruskal.test(value ~ Sample, data = input.long20)
pairwise.wilcox.test(input20$value, input20$Sample, p.adjust.method = "BH")

kruskal.test(value ~ Sample, data = input35)
pairwise.wilcox.test(input35$value, input35$Sample, p.adjust.method = "BH")

##### Sample change overtime ####
kruskal.test(value ~ variable, data = SC)
pairwise.wilcox.test(SC$value, SC$variable, p.adjust.method = "BH")

kruskal.test(value ~ variable, data = AC)
pairwise.wilcox.test(AC$value, AC$variable, p.adjust.method = "BH")

kruskal.test(value ~ variable, data = E25)
pairwise.wilcox.test(E25$value, E25$variable, p.adjust.method = "BH")

kruskal.test(value ~ variable, data = UV7)
pairwise.wilcox.test(UV7$value, UV7$variable, p.adjust.method = "BH")


#############################################################################################################################
############################               OLD              #################################################################
############################ Statistical analysis in a loop #################################################################
#############################################################################################################################
input.long$variable <- factor(input.long$variable, levels =c("0","5","11","20","35"))
#################### Sort the data by Timepoint#################################
#input.long <- input.long[order(variable),] # Order input.long by the variable column
timepoints <- c("0","5","11","20","35") # Define timepoints as a vector of strings
#timepoints <- c("0","35") # Define timepoints as a vector of strings

# Define the sample subsets
#sample_subsets <- list(c("SC", "Ac"),c("SC", "E25"),c("SC", "UV7"),c("Ac", "E25"), c("Ac", "UV7"), c("E25", "UV7"))
sample_subsets <- list(c("Ac", "E25"), c("Ac", "UV7"), c("E25", "UV7"))

# Initialize an empty data.table to store the results
results_time <- data.table(Timepoint = integer(),
                           Sample1 = character(),    # Sample 1 or 2 referr to position 1 or 2 in each item of the sample_subsets list
                           Sample2 = character(),
                           Anova_pvalue = numeric(),
                           TukeyHSD_pvalue = numeric())

# Loop over timepoints and sample subsets
for (i in seq_along(timepoints)) {
  for (subset in sample_subsets) {
    
    # Extract the samples in the subset
    subset_data <- input.long[variable == timepoints[i] & Sample %in% subset]
    # Check that the subset contains both samples and at least two unique values for "Sample"
    if (length(unique(subset_data$Sample)) != 2 | nlevels(factor(subset_data$Sample)) < 2) next
    # Check that the subset contains both samples and has at least 2 levels
    # Check if the subset contains any non-NA values for each sample
    #if (all(sapply(subset, function(x) all(is.na(subset_data$value[subset_data$Sample == x]))))) {
    #  next
    #}
    print(subset_data)
    # Check that the subset contains both samples
    #if (length(unique(subset_data$Sample)) < 2) next
    # Perform the ANOVA
    anova_result <- aov(value ~ Sample, data = subset_data)
    # Perform Tukey's HSD test
    tukeyhsd_result <- TukeyHSD(anova_result)
    # Extract the p-value from the ANOVA and TukeyHSD results
    anova_pvalue <- summary(anova_result)[[1]][["Pr(>F)"]][1]
    tukeyhsd_pvalue <- tukeyhsd_result[[1]][, "p adj"][1]
    # Add a row to the results data.table with the summary statistics and significance level
    if (tukeyhsd_pvalue < 0.0001) {
      sig <- "****"
    } else if (tukeyhsd_pvalue < 0.001) {
      sig <- "***"
    } else if (tukeyhsd_pvalue < 0.01) {
      sig <- "**"
    } else if (tukeyhsd_pvalue < 0.05) {
      sig <- "*"
    } else {
      sig <- "ns"
    }
    # Add a row to the results data.table with the summary statistics
    results_time <- rbind(results_time, data.table(Timepoint = timepoints[i], Sample1 = subset[1], Sample2 = subset[2], Anova_pvalue = anova_pvalue, TukeyHSD_pvalue = tukeyhsd_pvalue, TukeyHSD_sig = sig), fill = TRUE)
  }
}
# Save results as a CSV file
write.csv(results_time, "Stats_results_sample_differences_all.csv", row.names = FALSE)


#################### Sort the data by Sample #################################
samples <- c("SC","Ac", "E25", "UV7") # Define samples as a vector of strings
samples <- c("Ac", "E25", "UV7") # Define samples as a vector of strings

# Define the timepoint subsets
timepoint_subsets <- list(c("0", "5"), c("0", "11"), c("0", "20"), c("0", "35"), c("5", "11"), c("5", "20"), c("5", "35"), c("11", "20"), c("11", "35"), c("20", "35"))
timepoint_subsets <- list(c("0", "35"))

# Initialize an empty data.table to store the results
results_sample <- data.table(Sample = character(),
                             Sample1 = character(),  # Sample 1 or 2 referr to position 1 or 2 in each item of the timepoints_subsets list
                             Sample2 = character(),
                             Anova_pvalue = numeric(),
                             TukeyHSD_pvalue = numeric())

# Loop over samples and timepoint_subsets
for (i in seq_along(samples)) {
  for (subset in timepoint_subsets) {
    # Extract the timepoints in the subset
    subset_data <- input.long[Sample == samples[i] & variable %in% subset]
    
    # Check that the subset contains both samples
    if (length(unique(subset_data$variable)) != 2) next
    # Perform the ANOVA
    anova_result <- aov(value ~ variable, data = subset_data)
    # Perform Tukey's HSD test
    tukeyhsd_result <- TukeyHSD(anova_result)
    # Extract the p-value from the ANOVA and TukeyHSD results
    anova_pvalue <- summary(anova_result)[[1]][["Pr(>F)"]][1]
    tukeyhsd_pvalue <- tukeyhsd_result[[1]][, "p adj"][1]
    # Add a row to the results data.table with the summary statistics and significance level
    if (tukeyhsd_pvalue < 0.0001) {
      sig <- "****"
    } else if (tukeyhsd_pvalue < 0.001) {
      sig <- "***"
    } else if (tukeyhsd_pvalue < 0.01) {
      sig <- "**"
    } else if (tukeyhsd_pvalue < 0.05) {
      sig <- "*"
    } else {
      sig <- "ns"
    }
    # Add a row to the results data.table with the summary statistics
    results_sample <- rbind(results_sample, data.table(Sample = samples[i], Sample1 = subset[1], Sample2 = subset[2], Anova_pvalue = anova_pvalue, TukeyHSD_pvalue = tukeyhsd_pvalue, TukeyHSD_sig = sig), fill = TRUE)
  }
}

# Save results as a CSV file
write.csv(results_sample, "Stats_results_change_over_time.csv", row.names = FALSE)


