library(ggplot2)
library(magrittr)
library(ggpubr)
library(data.table)
library(RColorBrewer)
library(dplyr)
library(stringr)
library(tidyverse)
library(rstatix)
library(vegan)
library("ggdensity")
library(ggbreak)

# remove all the objects from the R session
rm(list=ls())#read input file, containing culture data


input<-fread("Bact_cDNA_all_for_r.csv", sep=',', header=TRUE, fill=TRUE)

input.long <- melt(input, measure.vars = c("0","5","11","20","35"))

input_summary <- input.long %>% # the names of the new data frame and the data frame to be summarised
  group_by(Sample, variable) %>%   # the grouping variable
  summarise(mean_value = mean(value, na.rm = TRUE),  # calculates the mean of each group
            sd_value = sd(value, na.rm = TRUE), # calculates the standard deviation of each group
            n_value = sum(!is.na(value)),  # calculates the sample size per group
            SE_value = sd(value, na.rm = TRUE)/sqrt(n())) # calculates the standard error of each group

#input_summary$variable <- factor(input_summary$variable, levels =c("0","5","11","20","35"))
#input.long$variable <- factor(input.long$variable, levels =c("0","5","11","20","35"))

# Convert 'Sample' column to a factor with specified levels
input_summary$Sample <- factor(input_summary$Sample, levels =c("SC","Ac", "E25","UV7"))
input.long$Sample <- factor(input.long$Sample, levels =c("SC","Ac", "E25","UV7"))
# Create a new column for the variable in the Amoeba data, to later have a continuous x axis
#input_summary$variable <- c(0,1,2,4,10,0,1,2,4,10,0,1,2,4,10)

# Create a new column for the variable in the Bacteria data, to later have a continuous x axis
input_summary$variable <- c(0,5,11,20,35,0,5,11,20,35,0,5,11,20,35,0,5,11,20,35)
input.long$variable <- c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,20,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35)
# Remove rows with NA in 'value'
#input.long <- na.omit(input.long, cols = "value")

########################################################################################################################
########################################################################################################################
########          GGPLOT        ########################################################################################
########################################################################################################################

########################################################################
####################### Lineplot #######################################
########################################################################
mb <- unique(as.numeric(c(1,5) %o% 10 ^ (9:13))) # y axis breaks



#image=
  # Creates a ggplot object with input_summary data, mapping variable to x-axis, mean_value to y-axis, Sample to group and Sample to colour with width of 0.75 units
  ggplot(input_summary, aes(x=variable, y=mean_value, group=Sample,  colour=Sample, width = 0.75))  + 
  
  # Adds a line layer to the plot with x and y mappings from the ggplot call, dodging the lines so they don't overlap using a width of 0.2, and setting the size of the line to 1.5 units
  geom_line(aes(x=variable, y=mean_value), position=position_dodge(width=0.2), size=1.5)+
  
  # Adds a point layer to the plot with Sample mapping to shape, dodging the points so they don't overlap using a width of 0.2, and setting the size of the points to 5 units
  geom_point(aes(shape=Sample), position=position_dodge(width=0.2), size=2)+
  geom_jitter(data=input.long, aes(x=variable, y=value,  colour=Sample, shape=Sample), size=3, width=1.2, height=0.2)+
  
  #geom_boxplot () +#(aes(x= variable, y=mean_value, fill=Sample, ymin = mean_value - sd_value, ymax = mean_value + sd_value), size =2, stat="identity",position=position_dodge(width = 0.8)) + # , width=0.5, position=position_dodge(1)
  #geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(), binwidth = 0.2) +#geom_jitter(shape=5, position=position_jitter(0.1)) +
  #geom_col(position = position_dodge(width = 0.8)) +  #The function geom_col uses the value of the y variable (mean_value) as the height of the bars
  #geom_smooth(method=lm, level=0.95)+
  
  #Adding error bars to the plot with specified range
  geom_errorbar(aes(ymin = mean_value, ymax = mean_value + sd_value),  width=2.0, size=0.5, position=position_dodge(width=0.2))+
  #scale_fill_manual(values = c("grey15","steelblue3", "orange2","firebrick3"))+
  scale_color_manual(values = c("grey15","steelblue3", "orange2","firebrick3"))+
  ylab("copies rRNA *g-1") +
  #Setting the y-axis to use a logarithmic scale with 7 log breaks
  scale_y_continuous(trans='log2', breaks = scales::log_breaks(n = 5))+# breaks =c(0,50000, 1000000,20000000,250000000))+
  scale_x_continuous(breaks =c(0,5,11,20,35))+
  #scale_y_continuous(breaks =mb)+
  #facet_grid(. ~ variable, scales="free") + 
  #guides(size = "none") +
  #annotate("text", x = 1.81, y = 250, label =paste ("**")) +
  theme_bw() +  # makes a very simplified background 
  # set font sizes for axis titles and labels
  theme(axis.title.x  = element_text(angle=90, vjust=2.2, size=25)) +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=25)) +
  theme(axis.text.y  = element_text(angle=90, vjust=0.5, size=20)) +
  theme(axis.title.y  = element_text(angle=90, vjust=2.2, size=25)) +
  #labs(y = "µg N/g dry soil/day", size =15) +
  ggtitle("Bacterial rRNA") +
  theme(plot.title = element_text(color="black", size=20, face="bold.italic"))+
  labs(x= "")+
  labs(fill="Sample") +
  theme(legend.key.size = unit(1.0, 'cm'))+
  #theme(legend.position="bottom")+
  theme(legend.position="none")+
  theme(legend.text = element_text(size=10))+
  theme(legend.title = element_blank())#+ ### to remove only the Legend title
#    scale_y_break(c(10000000, 1000000000), scales = 2)  # Include both below and above the break

ggsave(file="test.svg", plot=image, width=6, height=5)
###################################################################################################
####################### BOXPLOT ###################################################################
###################################################################################################
shape <- c(16,17, 15,3)

#image=
ggplot(input.long, aes(variable, value,  fill=Sample, width = 0.75))  + 
  scale_x_continuous(breaks = c(0,5,11, 20,35)) +
  geom_boxplot(
    data = input.long,
    aes(x = variable, y = value, fill = Sample, group = Sample),
    position = position_dodge(width = 2),
    outlier.shape = NA,
    alpha = 0.5, #fill = NA,
    width = 1.5
  )+
  geom_jitter(
    data = input.long, 
    aes(x = variable, y = value, colour = Sample, shape = Sample), 
    position = position_jitterdodge(dodge.width = 2, jitter.width = 0.2), 
    size = 3
  )+
  scale_shape_manual(values = shape)+
  scale_color_manual(values = c("grey15","steelblue3", "orange2","firebrick3")) +
  scale_fill_manual(values = c("grey15","steelblue3", "orange2","firebrick3"))+
  ylab("copies total / g soil") +
  scale_y_continuous(trans='log2', breaks = scales::log_breaks(n = 5))+
  facet_grid(. ~ variable, scales="free") + 
  guides(size = "none") +
  theme_bw() +  
  ggtitle("Bacterial rRNA") +
  theme(axis.text.x  = element_blank()) +
  #theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.x  = element_text(angle=90, vjust=2.2, size=15)) +
  theme(axis.text.y  = element_text(angle=90, vjust=0.5, size=15)) +
  theme(axis.title.y  = element_text(angle=90, vjust=2.2, size=15)) +
  labs(x= "")+
  labs(fill="Sample") +
  theme(legend.key.size = unit(1.6, 'cm'))+
  theme(legend.position="bottom")+
  theme(legend.position="none")+
  theme(legend.title = element_blank()) ### to remove only the Legend title




ggsave(file="boxtest.svg", plot=image, width=8, height=4)

###### Boxplot not facets ######
shape <- c(16,17, 15,3)

image=
  ggplot(input.long, aes(variable, value,  fill=Sample, width = 1))  + 
  geom_line(
    data = input_summary,
    aes(x=variable, y=mean_value, colour =Sample), position=position_dodge(width=3), 
    size=0.8)+
  scale_x_continuous(breaks = c(0,5,11, 20,35)) +
    scale_y_continuous(
      trans = 'log2',
      breaks = c(1e7, 1e8, 1e9, 1e10, 1e11)#,
     # labels = scales::label_number(scale_cut = scales::cut_short_scale())
    )+
  #scale_y_continuous(trans='log2', breaks = scales::log_breaks(1e7,1e11))+
  #scale_y_continuous(breaks = c(1e7,1e8,1e9,1e10,1e11)) +
  #scale_y_continuous(limits = c(1e7,1e11), breaks = seq(1e7,1e11, by=100))+
  geom_boxplot(
    data = input.long,
    aes(x = variable, y = value, fill = Sample, group = interaction(Sample, variable)),
    position = position_dodge(width = 3),
    outlier.shape = NA,
    alpha = 0.5, #fill = NA,
    width = 2
  )+
  geom_jitter(
    data = input.long, 
    aes(x = variable, y = value, colour = Sample, shape = Sample), 
    position = position_jitterdodge(dodge.width = 3, jitter.width = 1.5), 
    size = 2
  )+
  scale_shape_manual(values = shape)+
  scale_color_manual(values = c("grey15","steelblue3", "orange2","firebrick3")) +
  scale_fill_manual(values = c("grey15","steelblue3", "orange2","firebrick3"))+
  ylab("copies total / g soil") +
  #facet_grid(. ~ variable, scales="free") + 
  guides(size = "none") +
  theme_bw() +  
  ggtitle("Bacterial rRNA") +
  theme(axis.text.x  = element_blank()) +
  theme(axis.text.x  = element_text(angle=0, vjust=0.5, size=15)) +
  theme(axis.title.x  = element_text(angle=90, vjust=2.2, size=15)) +
  theme(axis.text.y  = element_text(angle=90, vjust=0.5, size=15)) +
  theme(axis.title.y  = element_text(angle=90, vjust=2.2, size=15)) +
  labs(x= "")+
  labs(fill="Sample") +
  theme(legend.key.size = unit(1.6, 'cm'))+
  theme(legend.position="bottom")+
  theme(legend.position="none")+
  theme(legend.title = element_blank()) ### to remove only the Legend title




ggsave(file="boxtest.svg", plot=image, width=8, height=4)



#############################################################################################################################
############################ Statistical analysis in a loop #################################################################
#############################################################################################################################

#################### Sort the data by Timepoint#################################
#input.long <- input.long[order(variable),] # Order input.long by the variable column
timepoints <- c("0","5","11","20","35") # Define timepoints as a vector of strings

# Define the sample subsets
sample_subsets <- list(c("SC", "Ac"),c("SC", "E25"),c("SC", "UV7"),c("Ac", "E25"), c("Ac", "UV7"), c("E25", "UV7"))

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
    # Check that the subset contains both samples
    if (length(unique(subset_data$Sample)) != 2) next
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
write.csv(results_time, "Stats_results_sample_differences.csv", row.names = FALSE)


#################### Sort the data by Sample #################################
samples <- c("SC","Ac", "E25", "UV7") # Define samples as a vector of strings

# Define the timepoint subsets
timepoint_subsets <- list(c("0", "5"), c("0", "11"), c("0", "20"), c("0", "35"), c("5", "11"), c("5", "20"), c("5", "35"), c("11", "20"), c("11", "35"), c("20", "35"))

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

########################################################################################################################
###### Statistical analysis by hand #################################################################################
########################################################################################################################
###### 
input0 <-  input.long[input.long$variable == "T0",c(1,2,3)]
input5 <-  input.long[input.long$variable == "T5",c(1,2,3)]
input11 <-  input.long[input.long$variable == "T11",c(1,2,3)]
input20 <-  input.long[input.long$variable == "T20",c(1,2,3)]
input35 <-  input.long[input.long$variable == "T35",c(1,2,3)]

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

statSC <- aov(value ~ variable, data = SC)
summary(statSC)
TukeyHSD(statSC)
rstatix::tukey_hsd(statSC)

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
shapiro.test(input11$value) #  normal
shapiro.test(input20$value) # not normal
shapiro.test(input35$value) #  normal

shapiro.test(AC$value) # not normal
shapiro.test(E25$value) # not normal
shapiro.test(UV7$value) # not normal
shapiro.test(SC$value) # not normal
shapiro.test(Prot$value) # not normal


shapiro.test(input.long$value) # not normal

######## Compare means of multiple groups using Kruskal-wallis (for non-normal data) #################################

kruskal.test(value ~ Sample, data = input0)
pairwise.wilcox.test(input0$value, input0$Sample, p.adjust.method = "BH")

kruskal.test(value ~ Sampe, data = input5)
pairwise.wilcox.test(input5$value, input5$Sample, p.adjust.method = "BH")

kruskal.test(value ~ Sample, data = input11)
pairwise.wilcox.test(input11$value, input11$Sample, p.adjust.method = "BH")

kruskal.test(value ~ Sample, data = input.long20)
pairwise.wilcox.test(input20$value, input20$Sample, p.adjust.method = "BH")

ekruskal.test(value ~ Sample, data = input35)
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

########## check for normal distribution of your data  #################
#the density plot provides a visual judgment about whether the distribution is bell shaped.
ggdensity(input.long$value, 
          main = "Density plot",
          xlab = "bla")

#Q-Q plot: Q-Q plot (or quantile-quantile plot) draws the correlation between a given sample and the normal distribution
ggqqplot(input.long$value)

#shapiro.test() can be used to perform the Shapiro-Wilk test of normality for one variable (univariate)
shapiro.test(input.long$value)
#p-value > 0.05 implying that the distribution of the data are not significantly different from normal distribution.





