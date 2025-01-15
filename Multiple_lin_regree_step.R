# Mulitple linear regression using step
# This code performs a multiple linear regression analysis with and without stepwise selection on abundance data from a microbiome dataset. 

# Extract the OTU (abundance) table
otu_table <- otu_table(gut_filtered_abundance)

# Convert to a data frame and transpose so that sample names are rows, gene/species names are columns
otu_df <- as.data.frame(t(otu_table))

# Extract sample metadata and convert to a data frame
sample_df <- as.data.frame(sample_data(gut_filtered_abundance))

# Extract the taxonomy table and convert to a data frame
taxonomy_df <- as.data.frame(tax_table(gut_filtered_abundance))

# Ensure that the row names of taxonomy_df match OTU names in otu_df
species_names <- taxonomy_df$Species

# Set the column names of otu_df to the species names
colnames(otu_df) <- species_names

# Combine sample data and OTU data
combined_df <- cbind(sample_df, otu_df)

# Check the first few rows to ensure correct merging
head(combined_df)

# Focus on one species
species_of_interest <- "s__Bifidobacterium_longum"  # Replace with name of one species of interest

# Create desired formula, e.g. below:
full_formula <- as.formula(paste(species_of_interest, "~", "Dietary_sugar", "+ Timepoint + Condition + Timepoint:Condition"))

# Perform linear regression without stepwise selection
lm_no_step <- lm(full_formula, data = combined_df)

# Perform stepwise selection
lm_step <- step(lm_no_step, direction = "both", trace = 0)

# Summarize the results for both models (step and no step)
summary_no_step <- summary(lm_no_step)
summary_step <- summary(lm_step)

# Comparison
summary_no_step
summary_step

# Store the results in tidy format
tidy_no_step <- broom::tidy(lm_no_step)
tidy_step <- broom::tidy(lm_step)

# Add model information
tidy_no_step$model <- "No Stepwise"
tidy_step$model <- "Stepwise"

# Combine the results for comparison
comparison_df <- rbind(tidy_no_step, tidy_step)

# View the results
print(comparison_df)

# Save the results to a CSV file
write.csv(comparison_df, file = "comparison_species_X.csv", row.names = FALSE)
