using Plots, Random, Distributions, Colors, LaTeXStrings, Printf, TickTock, JLD2

include("auxiliary_functions_v28.jl") # Pull in functions 
Random.seed!(12345678)

###################################
## Pull data from benchmark case ##
###################################

file = jldopen("scenarios_spatial.jld2","r")
iso_data = file["scenario_isotropic"]
close(file)

########################################
## Convert data to initial conditions ##
########################################

# List of initial conditions X0 for all patients (20 in this case)
X0s = data_to_X0s(iso_data)

# mutation parameters (NOTE: number of elements of u is number of driver mutations minus 1)
s = 0.004           # Selective Advantage for Driver mutations
u = [3.4e-5,3.4e-5,3.4e-5,3.4e-5,3.4e-5,3.4e-5,3.4e-5,3.4e-5,3.4e-5] # driver mutation probabilities 
# NOTE: number of elements of u is number of driver mutations - 1

# spatial distribution
L = 1.
σ₁ = 1*L
σ₂ = 1*L
dx = mx -> Normal(mx,σ₁)  # spatial distribution x
dy = my -> Normal(my,σ₂)  # spatial distribution y

# fertility distribution
α = 1/2.
β = 1.
shift = 10

# Different types of fertility functions
g = t -> 1/t
f = (x,y,t) -> 1-exp(-β*(x^2+y^2)*g(t))


treatment_time = 1 # time in months


function f_chemo(x,y,t;treatment_time=treatment_time,k=k)
    chemo_day = treatment_time*30
    if (mod(t,chemo_day)<=k-1)
        return  f(x,y,t)
    else
        return 1
    end
end

# time parameters
T = 10.        # Total time in Years
k = 3.               # Time step in days - Generation time

# other parameters
extinction = true  # Allow process to go extinct/condition on survival
Nrec = 0 # Count number of recurrences



###############
## Pull Data ##
###############
# Number of runs or patients to simulate

# Runes benchmark() Nav number of runs and retunrns in array "data"
# data = [[X,Y,L,Nrec,Nm,Nt], ... ]
tick()
# data, Nav = PullDataBenchmarkTherapy(X0s,s,u,dx,dy,f,T,k,extinction, Nrec)
tock()

file = jldopen("scenarios_treatments.jld2","r")
data = file["scenario_chemo_treatment"]
close(file)
#####################
## POST PROCESSING ##
#####################

savefig=true

nx = 101 # Number of x bins 
ny = 101 # Number of y bins

# post_process_chemo(data,X0s,nx,ny,f_chemo,savefig)

