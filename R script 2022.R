# Install required packages if not already installed
install.packages("readxl")
install.packages("vegan")
install.packages("ggplot2")
install.packages("reshape2")

# Load libraries
library(readxl)
library(vegan)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(reshape2)

#Call my data file
data <- read.csv("C:/Users/gomez.428/OneDrive - The Ohio State University/Documents/Cesia/2022 ant data.csv")
data <- read.csv("./2022 ant data.csv")

data

str(data) # need to change categorical variables to factors
data$Site <- as.factor(data$Site)
data$Plot_no <- as.factor(data$Plot_no)
data$Year <- as.factor(data$Year)
data$Collection <- as.factor(data$Collection)
str(data)

# Inspect data
head(data)

data <- as.data.frame(data)

# Data Preparation
data_abundance <- data %>%
  dplyr::select(-Site, -Plot, -Year, -Collection, -Date_set, -Date_collected)

# Convert all columns to numeric for abundance calculation
data_abundance <- data_abundance %>%
  mutate(across(everything(), as.numeric))

# Diversity Analysis (Shannon and Simpson Index)
shannon_index <- diversity(data_abundance, index = "shannon")
simpson_index <- diversity(data_abundance, index = "simpson")

library(hillR)
data$shan <- hill_taxa(data[,8:20], q = 1, MARGIN = 1) # shannon

hist(data$shan)
dotchart(data$shan, group = data$Site, pch = 19)
dotchart(data$shan, group = data$Year, pch = 19)
boxplot(shan ~ Site, data = data)
stripchart(shan ~ Site, data = data, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

boxplot(shan ~ Site * Year, data = data)
stripchart(shan ~ Site * Year, data = data, pch = 19, add = TRUE,
           vertical = TRUE, method = "jitter", jitter = 0.2)

### Run models
library(lme4)
library(car)
library(lmerTest)
library(emmeans)

mod.shan <- lmer(shan ~ Site * Year + (Collection|Plot_no), data = data)
summary(mod.shan)
Anova(mod.shan, type = "III")
plot(mod.shan, pch = 19)
qqnorm(residuals(mod.shan))
qqline(resid(mod.shan))
# No differences in genera diversity

# Adding diversity indices to the original data
analysis_data <- data %>%
  mutate(Shannon_Index = shannon_index,
         Simpson_Index = simpson_index)

# Compare Diversity Metrics Across Sites
shannon_model <- aov(Shannon_Index ~ Site, data = analysis_data)
simpson_model <- aov(Simpson_Index ~ Site, data = analysis_data)

# Summarize and Display ANOVA results
summary(shannon_model)
summary(simpson_model)

# Plotting Diversity Indices
p1 <- ggplot(analysis_data, aes(x = Site, y = Shannon_Index, fill = Site)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Shannon Diversity Index by Site", y = "Shannon Index")

p2 <- ggplot(analysis_data, aes(x = Site, y = Simpson_Index, fill = Site)) +
  geom_boxplot() +
  theme_minimal() +
  labs(title = "Simpson Diversity Index by Site", y = "Simpson Index")

# Display plots
p1
p2

# Abundance Analysis (PERMANOVA)
ant_matrix <- as.matrix(data_abundance)  # Create a matrix of only genera counts

# Remove rows where all genera counts are zero
rows_to_keep <- rowSums(ant_matrix) > 0
ant_matrix_clean <- ant_matrix[rows_to_keep, ]

# Also clean the original dataset to match the filtered matrix
data_abundance_clean <- data[rows_to_keep, ]  

# Create a distance matrix with the cleaned data
dist_matrix <- vegdist(ant_matrix_clean, method = "bray")

# Run PERMANOVA
permanova <- adonis2(dist_matrix ~ Site, data = data_abundance_clean)
print(permanova)


# Abundance Comparison for Each Genus
abundance_data <- data_abundance_clean %>%  # Use the cleaned dataset
  pivot_longer(cols = Aphaenogaster:Aconthomyops, names_to = "Genus", values_to = "Count")

# Plotting Abundance by Site
p3 <- ggplot(abundance_data, aes(x = Site, y = Count, fill = Site)) +
  geom_boxplot() +
  facet_wrap(~ Genus, scales = "free_y") +
  theme_minimal() +
  labs(title = "Ant Genera Abundance Across Sites", y = "Abundance") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display Plot
print(p3)


# Sum the counts of each genus across all plots to see what genus of ant is most abundant
genus_abundance <- abundance_data %>%
  group_by(Genus) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE)) %>%
  arrange(desc(Total_Count))

# Display the result
print(genus_abundance)

# Summarize the abundance of each genus
genus_abundance <- abundance_data %>%
  group_by(Genus) %>%
  summarise(Total_Count = sum(Count, na.rm = TRUE)) %>%
  arrange(desc(Total_Count))

# Plotting the bar graph
p4 <- ggplot(genus_abundance, aes(x = reorder(Genus, -Total_Count), y = Total_Count, fill = Genus)) +
  geom_bar(stat = "identity") +
  theme_minimal() +
  labs(title = "Total Abundance of Ant Genera Across All Sites",
       x = "Genus",
       y = "Total Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Display the plot
print(p4)

# Prepare the data for NMDS
ant_matrix_nmds <- data_abundance_clean %>%
  dplyr::select(where(is.numeric)) %>%  # Select only numeric columns
  as.matrix()  # Convert to matrix

# Remove rows with zero counts for NMDS calculation
ant_matrix_nmds <- ant_matrix_nmds[rowSums(ant_matrix_nmds) > 0, ]

# Perform NMDS (Bray-Curtis Distance)
set.seed(123)
nmds_result <- metaMDS(ant_matrix_nmds, distance = "bray", k = 2, trymax = 100)

# Extract NMDS scores (site scores)
nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))

# Add Site information to NMDS scores
nmds_scores$Site <- data_abundance_clean$Site[rowSums(ant_matrix_nmds) > 0]

# Plot the NMDS Graph
p5 <- ggplot(nmds_scores, aes(x = NMDS1, y = NMDS2, color = Site)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "NMDS Plot of Ant Community Structure",
       x = "NMDS1",
       y = "NMDS2") +
  theme(legend.position = "right")

# Display the NMDS plot
print(p5)



###
colSums(data[,8:20])
rowSums(data[,8:20]) # have rows without any counts (i.e., no ants collected in a trap at a site), so we need to remove those before moving on

dat <- data[rowSums(data[,8:20])>0,]
rowSums(dat[,8:20])

dis.matrix <- vegdist(dat[,8:20], method = "bray")
dis.matrix

# run the nonmetric multidimensional scaling model
nmds.ants <- metaMDS(dis.matrix, trymax = 500, autotransform = TRUE, k = 2)
nmds.ants # stress is quality of fit
stressplot(nmds.ants)
plot(nmds.ants) # basic plot with no treatment distinctions

adonis2(dis.matrix ~ dat$Site * dat$Year, permutations = 999)
library(pairwiseAdonis)
pairwise.adonis(dis.matrix, dat$Site)

# plot the NMDS model
ordiplot(nmds.ants, disp = "sites", type = "n", xlim = c(-2, 2), ylim = c(-1, 1.5))
points(nmds.ants, dis = "sites", select = which(dat$Site=="Forest"), pch = 17, cex = 2, col = "#73D055FF")
points(nmds.ants, dis = "sites", select = which(dat$Site=="Salvaged"), pch = 16, cex = 2, col = "#481567FF")
points(nmds.ants, dis = "sites", select = which(dat$Site=="Windthrow"), pch = 15, cex = 2, col = "#2D708EFF")
levels(dat$Site)
ordiellipse(nmds.ants, dat$Site, draw = "lines", col = c("#73D055FF", "#481567FF", "#2D708EFF"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("bottomleft", legend = c("Forest", "Salvaged", "Windthrow"),
       pch = c(17, 16, 15), cex = 1.5, bty = "n", col = c("#73D055FF", "#481567FF", "#2D708EFF"))


pairwise.adonis(dis.matrix, dat$Year)

# plot the NMDS model
ordiplot(nmds.ants, disp = "sites", type = "n", xlim = c(-2, 2), ylim = c(-1, 1.5))
points(nmds.ants, dis = "sites", select = which(dat$Year=="2015"), pch = 17, cex = 2, col = "#22A884FF")
points(nmds.ants, dis = "sites", select = which(dat$Year=="2022"), pch = 16, cex = 2, col = "#FDE725FF")
levels(dat$Year)
ordiellipse(nmds.ants, dat$Year, draw = "lines", col = c("#22A884FF", "#FDE725FF"), 
            lwd = 3, kind = "sd", conf = 0.90, label = FALSE)

legend("bottomleft", legend = c("2015", "2022"),
       pch = c(17, 16, 15), cex = 1.5, bty = "n", col = c("#22A884FF", "#FDE725FF"))
