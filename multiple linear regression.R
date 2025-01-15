# Multiple linear regression
# Data initially stored as a Phyloseq obj
setwd()

# Load libraries
library(phyloseq)
library(dplyr)
library(ggplot2)
library(broom)

# Extract the OTU abundance data & convert to a df
otu_table <- otu_table(gut_filtered_abundance)
otu_df <- as.data.frame(otu_table)
otu_df <- t(otu_df)  # Transpose to have samples as rows

# Extract sample metadata
sample_data <- sample_data(gut_filtered_abundance)
sample_df <- as.data.frame(sample_data)

# Extract taxonomy table
taxonomy_table <- tax_table(gut_filtered_abundance)
taxonomy_df <- as.data.frame(taxonomy_table)

# Ensure that row names of taxonomy_df match OTU names in otu_df
species_names <- taxonomy_df$Species
colnames(otu_df) <- species_names

# Combine sample data and OTU data
combined_df <- cbind(sample_df, otu_df)
head(combined_df)
#################################################

# Define the trait for regression
trait_name <- "Saliva_pH"

# List of species names
species_names <- colnames(otu_df)

# Initialize a list to store xthe results
results_list <- list()

# Loop through each species and perform regression adjusted for Sex and Age
for (species in species_names) {
  # Create formula to include Timepoint, Condition, and possible interaction term
  formula <- as.formula(paste(species, "~", trait_name, "+ Sex + Age"))
  
  # Perform regression using lm
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

# Define FDR threshold and filter results to get significant species based on FDR-corrected p-values
fdr_threshold <- 0.05
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
