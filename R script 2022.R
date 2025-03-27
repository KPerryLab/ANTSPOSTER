# Install required packages if not already installed
install.packages("readxl")
install.packages("vegan")
install.packages("ggplot2")
install.packages("reshape2")

# Load libraries
library(readxl)
library(vegan)
library(ggplot2)
library(readxl)
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggpubr)
library(reshape2)

#Call my data file
data <- read.csv("C:/Users/gomez.428/OneDrive - The Ohio State University/Documents/Cesia/2022 ant data.csv")

data

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