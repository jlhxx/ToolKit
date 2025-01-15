#multiple linear regression

#linear regression sex and age
# Load necessary libraries
library(dplyr)
library(ggplot2)
library(broom)  # For tidy() function

# Extract the OTU (abundance) table
otu_table <- otu_table(gut_filtered_abundance)

# Convert to a data frame
otu_df <- as.data.frame(otu_table)
otu_df <- t(otu_df)  # Transpose to have samples as rows

# Extract sample metadata
sample_data <- sample_data(gut_filtered_abundance)

# Convert to a data frame
sample_df <- as.data.frame(sample_data)

# Extract taxonomy table
taxonomy_table <- tax_table(gut_filtered_abundance)

# Convert to a data frame
taxonomy_df <- as.data.frame(taxonomy_table)

# Ensure that row names of taxonomy_df match OTU names in otu_df
species_names <- taxonomy_df$Species

# Set the column names of otu_df to the species names
colnames(otu_df) <- species_names

# Combine sample data and OTU data
combined_df <- cbind(sample_df, otu_df)

# Check the first few rows to ensure correct merging
head(combined_df)
######################################

# Define the trait you want to regress on
trait_name <- "Fecal_pH"

# List of species names
species_names <- colnames(otu_df)

# Initialize a list to store results
results_list <- list()

# Loop through each species and perform regression adjusted for Timepoint, Condition, Sex, and Age
for (species in species_names) {
  # Create formula to include Timepoint, Condition, and possible interaction term
  formula <- as.formula(paste(species, "~", trait_name, "+ Sex + Age + Timepoint + Condition + Timepoint:Condition"))
  
  # Perform linear regression
  lm_result <- lm(formula, data = combined_df)
  
  # Store results in a tidy format
  tidy_result <- tidy(lm_result)
  tidy_result$species <- species
  results_list[[species]] <- tidy_result
}

# Combine all results into a single data frame
results_df <- do.call(rbind, results_list)

# Perform FDR correction on the p-values
results_df <- results_df %>%
  group_by(species) %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr"))

# Define the FDR p-value threshold
fdr_threshold <- 0.05

# Filter the results to get only significant species based on FDR-corrected p-values
significant_results <- results_df %>% 
  filter(p.adjusted < fdr_threshold) %>%
  select(species, term, estimate, p.value, p.adjusted)

# Format p-values for plotting
format_p_value <- function(p) {
  if (p < 0.001) {
    return("p.adj < 0.001")
  } else {
    return(sprintf("p.adj = %.3f", p))
  }
}

significant_results <- significant_results %>%
  mutate(p.label = sapply(p.adjusted, format_p_value))

# Plot the estimates of significant species after FDR correction
ggplot(significant_results, aes(x = reorder(species, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = p.adjusted < fdr_threshold), width = 0.5) +
  geom_text(aes(label = p.label), 
            vjust = -0.5, size = 3, nudge_y = 0.1, check_overlap = TRUE) +
  coord_flip() +
  labs(title = "Significant Species Regression Estimates (FDR Corrected)",
       x = "Species",
       y = "Estimate",
       fill = "Significant") +
  theme_minimal()
#####################################
# Initialize a list to store results
results_list <- list()

# Loop through each species and perform regression adjusted for Timepoint, Condition, Sex, Age, and their interaction
for (species in species_names) {
  # Create formula to include Timepoint, Condition, and their interaction term
  formula <- as.formula(paste(species, "~", trait_name, "+ Sex + Age + Timepoint + Condition + Timepoint:Condition"))
  
  # Perform linear regression
  lm_result <- lm(formula, data = combined_df)
  
  # Store results in a tidy format
  tidy_result <- tidy(lm_result)
  tidy_result$species <- species
  results_list[[species]] <- tidy_result
}

# Combine all results into a single data frame
results_df <- do.call(rbind, results_list)

# Perform FDR correction on the p-values
results_df <- results_df %>%
  group_by(species) %>%
  mutate(p.adjusted = p.adjust(p.value, method = "fdr"))

# Define the FDR p-value threshold
fdr_threshold <- 0.05

# Filter the results to get only significant species based on FDR-corrected p-values
significant_results <- results_df %>% 
  filter(p.adjusted < fdr_threshold) %>%
  select(species, term, estimate, p.value, p.adjusted)

# Format p-values for plotting
format_p_value <- function(p) {
  if (p < 0.001) {
    return("p.adj < 0.001")
  } else {
    return(sprintf("p.adj = %.3f", p))
  }
}

significant_results <- significant_results %>%
  mutate(p.label = sapply(p.adjusted, format_p_value))

# Plot the estimates of significant species after FDR correction
ggplot(significant_results, aes(x = reorder(species, estimate), y = estimate)) +
  geom_bar(stat = "identity", aes(fill = p.adjusted < fdr_threshold), width = 0.5) +
  geom_text(aes(label = p.label), 
            vjust = -0.5, size = 3, nudge_y = 0.1, check_overlap = TRUE) +
  coord_flip() +
  labs(title = "Significant (FDR) Species Regression Estimates Fecal pH (+ Sex + Age + Timepoint + Condition + Timepoint:Condition)",
       x = "Species",
       y = "Estimate",
       fill = "Significant") +
  theme_minimal() +
  theme(plot.title = element_text(size = 10))  # Adjust size as needed
