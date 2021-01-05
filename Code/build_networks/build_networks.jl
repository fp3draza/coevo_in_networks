# Define packages to use
using Distributed, Distributions, Random, DelimitedFiles

@everywhere function niche_networks(num_sp_in_A,num_sp_in_B,empirical_connectance)
    #=
    This function builds bipartite networks with a specified size, defined as the  
    sum of species belonging to guild A and guild B, and with a specified 
    target connectance following the niche algorithm. It requires the following parameters:
	- num_sp_in_A: int, the number of species belonging to guild A
	- num_sp_in_B: int, the number of species belonging to guild B
	- empirical_connectance: float, the target connectance value
    =#

    # Initialize empty matrix
    niche_matrix =  Matrix{Int}(undef, num_sp_in_A, num_sp_in_B)

    # Initialize uniform distribution
    axis_distribution = Uniform(0,1)

    # Position species in group A along niche axis
    axis_a = rand(axis_distribution,num_sp_in_A)
    axis_b = rand(axis_distribution,num_sp_in_B)

    # Define Beta
    beta_value = (1/(2*empirical_connectance)) - 1

    # Define Beta distribution
    beta_distribution = Beta(1,beta_value)

    # Draw from Beta
    beta_vector = rand(beta_distribution,num_sp_in_A)

    # Niche range for species in axis A
    range_vector = axis_a .* beta_vector

    # Calculate niche halved
    range_vector_halved = range_vector/2
    range_vector_halved_minus_one = 1 .- (range_vector/2)

    # Determine maximum value for range
    max_vec = map(min,axis_a,range_vector_halved_minus_one)

    # Initialize vector for center values for range
    range_center = zeros(num_sp_in_A)

    # Draw range center
    for i in 1:num_sp_in_A

        # Determine the center of range for each species
        range_center[i] = rand(range_vector_halved[i]:0.1:max_vec[i])

	# Draw interactions if niche of species in axis B overlaps with the diet range of species in guild A
        niche_matrix[i,:] = (axis_b .> (range_center[i] - range_vector_halved[i])) .&
                            (axis_b .< (range_center[i] + range_vector_halved[i]))

    end

    # Return matrix
    return niche_matrix

end


@everywhere function build_niche_networks(connectance)
    #=
    This function builds niche networks for a range of connectance values.
    As currently defined, the function will build networks with 60 species,
    building 20 replica networks per specified connectance values.
    The networks will be stored as both csv and txt files for further analyses.
    The function requires the following parameters:
	- connectance: float, the target connectance value 
    =#

    # Parameters
    
    # Obtain current directory, assumes you are running from the directory where this file is stored
    current_directory = @__DIR__
    
    # Define where networks will be stored
    folder_to_save_csv = string(current_directory,"/../../Data/networks/")
    folder_to_save_txt = string(current_directory,"/../../Data/for_modular/")
    
    # Determine network size and replicas
    num_sp_in_A = 30
    num_sp_in_B = 30
    replicas = 20

    # Cycle through the replcias for each connectance value
    for replica in 1:replicas
        # Create a network following the niche model
        niche_matrix = niche_networks(num_sp_in_A,num_sp_in_B,connectance)
        # If a network is disconnected
        while any(sum(niche_matrix,dims = 2) .== 0) || any(sum(niche_matrix,dims = 1) .== 0)
	    # Redraw the network
            niche_matrix = niche_networks(num_sp_in_A,num_sp_in_B,connectance)
        end
	# Move to folder where network will be stored
        cd(folder_to_save_csv)
        # Save network as csv file
        writedlm(string("C_",connectance,"_R_",replica,".csv"),niche_matrix,',')
        # Move to folder where network will be stored
        cd(folder_to_save_txt)
        # Save network as txt file
        writedlm(string("C_",connectance,"_R_",replica,".txt"),niche_matrix,' ')
    end
    return true
end

# This following section should only be run if the required folders have not been created

# Fill this section in to place where files should be created
#folder_to_save = ""

# Create directory for storing networks
#mkdir(string(folder_to_save,"/networks/"))
# Create directory for storing networks for MODULAR analysis
#mkdir(string(folder_to_save,"/for_modular/"))

# Create networks by running the function 
# Currently this will build networks ranging in target connectance from 0.05 to 0.49 by 0.02
# Note this will run in parallel, each connectance value will run as a separate process
result = pmap(build_niche_networks,collect(0.05:0.02:0.49))
