# Load required packages
require(dplyr)

merge_dataframe = function(directory){
  
  # This function merges the output of the trait matching analyses for all levels of mutualistic selection.
  # The function requires the following parameters:
  # directory: directory where trait matching results are stored
  
  # Create empty dataframe
  full_dataframe = NULL
  
  # List files in the provided directory
  list_of_files = list.files(directory)
  
  # Initialize counter
  counter = 1
  
  # Set the level of mutualistic selection
  odd_mutualistic_selection = 0.15
  even_mutualistic_selection = 0.1
  
  # Loop through the files in the directory
  for (file in list_of_files) {
    
    # Read the current result file
    current_dataframe = read.csv(paste(directory,file,sep = ''),row.names = 1)
    
    # Add column containing the level of mutualistic selection
    if (counter %% 2 != 0) {
      current_dataframe$mutualistic_selection = odd_mutualistic_selection
      # Update the level of mutualistic selection
      odd_mutualistic_selection = odd_mutualistic_selection + 0.1
    }
    
    if (counter %% 2 == 0) {
      current_dataframe$mutualistic_selection = even_mutualistic_selection
      # Update the level of mutualistic selection
      even_mutualistic_selection = even_mutualistic_selection + 0.1
    }
    
    # Join the dataframes
    full_dataframe = rbind(full_dataframe,current_dataframe)
    # Update counter
    counter = counter + 1
  }
  
  # Write merged file
  write.csv(full_dataframe,paste(directory,'merged_results.csv',sep = ''))
  # Print
  print('done')
}

# Define the directory where trait matching results are stored.
# Note that information on number and ratio of species in the networks are specified in the path (processed_output_30_30).
# If you want to work with results from other network sizes or ratios, change the numbers (processed_output_X_X)
trait_matching_directory = 'Results/processed_output_30_30/trait_matching/'

# Run function with the specified trait matching directory
merge_dataframe(trait_matching_directory)
