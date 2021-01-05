# Load required functionss
source('Code/measure_network_properties/network_roles_function.R')
source('Code/measure_network_properties/membership_function.R')

run_network_roles = function(){
  # This function computes the roles of species in the networks

  # Default parameter values
  min_connectance = 0.05
  max_connectance = 0.49
  connectance_by = 0.02
  min_replica = 1
  max_replica = 20

  # Initialize empty dataframe
  final_df = NULL

  # Loop through connectance and replica values
  for (connectance in seq(min_connectance,max_connectance,connectance_by)) {
    for (replica in min_replica:max_replica) {
      # Compute network roles for the current network with given connectance and replica
      current_result = compute_network_roles(connectance, replica)
      # Append output to the dataframe
      final_df = rbind(final_df,current_result)
    }
  }

  # Save dataframe
  write.csv(final_df,'Results/processed_output/network_roles/network_roles_niche.csv')
}

# Run the function
run_network_roles()
