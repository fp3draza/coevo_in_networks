# Define packages to use
using Random, Distributions, CSV, Distributed, Statistics, DelimitedFiles, DataFrames

@everywhere function run_simulations_on_HPC(m)
    #=
    This function simulates coevolution at a specified mutualistic selection strength.
    It requires the following parameters:
	- m: float, value of mutualistic selection strenght
    =#

    # Define parameters
   
    # Number of replica simulations to perform per network
    n_sim = 5
    # Alpha parameter
    alpha = 0.2
    # Phi mean, standard deviation and distribution
    phi_mean = 0.5
    phi_sd = 0.01
    phi_distribution = Normal(phi_mean, phi_sd)
    # mutualistic selection strength mean, standard deviation and distribution
    m_mean = m
    m_sd = 0.01
    m_distribution = Normal(m_mean, m_sd)
    # Trait min and max values and distribution
    theta_min = 0
    theta_max = 10
    theta_distribution =  Uniform(theta_min, theta_max)
    # Epsilon
    epsilon = 0.000001
    # Max number of generations
    t_max = 10000

    # Folder to store results
    # This will obtain the current directory, it assumes you are running this script form where it is located
    current_directory = @__DIR__
    # Define folder to store results
    folder =  string(current_directory,"/../../Results/uniform/")
    # Define folder where networks are stored
    network_directory = string(current_directory,"/../../Data/networks_20_20/")
   
    # Define metadata to store 
    alpha_string = replace(string(alpha),"." => "")
    phi_string = replace(string(phi_mean),"." => "")
    m_string = replace(string(m_mean),"." => "")
    theta_min_string = replace(string(theta_min),"." => "")
    theta_max_string = replace(string(theta_max),"." => "")

    # Create folder to store results and move there
    folder_to_store = string(folder,"all_networks","_m_",m_string)
    mkdir(folder_to_store)
    cd(folder_to_store)

    # List files in network directory and make a list
    net_files = readdir(network_directory)
    net_names = net_files

    # Cycle through the list of networks
    for i in 1:length(net_files)

        # Create a directory to store results of current network
       mkdir(replace(net_files[i],".csv"=>""))

        # Read network
        mat = readdlm(string(network_directory,net_files[i]),',',)[1:end,1:end]
       
        # Calculate network size
        n_row = size(mat)[1]
        n_col = size(mat)[2]
        n_sp = n_row + n_col

        # Build square adjacency matrix
        f = vcat(hcat(zeros(n_row,n_col),mat),hcat(transpose(mat),zeros(n_row,n_col)))

        ## Cycle through the number of simulations to perform per network
        for k in 1:n_sim

            # Sample phi values
            phi = rand(phi_distribution,n_sp)
	    # Make sure no values are negative
            while any(phi.<0)
                phi = rand(phi_distribution,n_sp)
            end
	
            # Sample m values
            m = rand(m_distribution,n_sp)
   	    # Make sure no values are negative
            while any(m.<0)
                m = rand(m_distribution,n_sp)
            end

            # Sample theta (environmental optima)
            theta = rand(theta_distribution,n_sp) 
            # Sample initial trait values
            init = rand(theta_distribution,n_sp) 

            ## Run simulations
            z_list = coevolution_model(n_sp,f,phi,alpha,theta,init,m,epsilon,t_max)

            # Extract results
            z_initial = z_list[1,:]
            z_final = z_list[size(z_list)[1],:]
            
            # Build dataframe
            z = DataFrame(transpose(hcat(theta,z_initial,z_final)))
            z = rename(z, vcat([Symbol("R$i") for i in 1:n_row],[Symbol("C$i") for i in 1:n_col]))
            z = insert!(z,1,["theta","z_initial","z_final"],:names)

            # Write result
            CSV.write(string(folder,"all_networks_m_",m_string,"/",replace(net_names[i],".csv"=>""),"/",
            net_names[i],"_m_",m_string,"_alpha",alpha_string,"_phi",phi_string,"_theta",theta_min_string,
            "-",theta_max_string,"_sim",k,".csv"),z)
        end
    end
end

@everywhere function coevolution_model(n_sp, f, phi, alpha, theta, init, m, epsilon, t_max)
  #=
  This function simulates coevolution in a given mutualistic network.
  It requires the following parameters:
	- n_sp: int, the total number of species in the network
	- f: matrix, the square adjacency matrix representing the mutualistic network
	- h: vector, the heritability value of traits
 	- alpha: float, the sensitivity of selection
 	- theta: vector, the environmental optima of trait values for all species
  	- init: vector, the initial trait values for all species
	- m: vector, the strenght of mutualistic selection
	- epsilon: float, a number used to determine when equilibrium is reached
	- t_max: int, the maximum number of timesteps to simuulate
  =#
  
  # Initialize empty matrix
  z_mat = zeros(t_max,n_sp)
  # Fill first row of matrix with initial trait values
  z_mat[1,:] = init

  # Initialize counter
  up_to = 0

  # Cycle throught the generations
  for t = 1:(t_max - 1)

    # Obtain current trait values
    z = z_mat[t,:]

    # Calculate difference in trait values between partners
    z_dif = transpose(f.*z) - f.*z

    # Calculate evolutionary effect of interactions
    q = f.*(exp.(-alpha.*(z_dif.^2)))
    q_n = q./sum(q,dims = 2)
    q_m = q_n.* m
    sel_dif = q_m .* z_dif
    r_mut = phi .* sum(sel_dif,dims=2)

    # Calculate evolutionary effect of environment
    r_env = phi .* (1 .- m) .* (theta .- z)

    # Update trait values given selection imposed by partners and environment
    z_mat[t+1,:] = z .+ r_mut .+ r_env

    # Calculate the difference in trait values between current and previous generation
    dif = mean(abs.(z .- z_mat[t+1,:]))

    # Update counter
    up_to = t

    # Determine if equilibrium has been reached
    if dif < epsilon
      break
    end
  end

  # Return matrix with trait values
  return z_mat[1:up_to,:]
end

# Run the simulations
# Currently, this will simulate coevolution for m ranging from
# 0.1 to 0.95 by 0.05.
# Note that this is currently implemented to run in parallel. Each mutualistic selection
# value will run as its own process.
result = pmap(run_simulations_on_HPC,collect(0.1:0.05:0.95))

