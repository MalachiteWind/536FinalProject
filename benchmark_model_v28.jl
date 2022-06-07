using Plots, Random, Distributions, Colors, LaTeXStrings, Printf, TickTock, JLD2


#####################################
## Name data that you want to save ##
#####################################

scenario = "scenario_treatment_8_benchmark_case" # File name

include("auxiliary_functions_v28.jl") # Pull in functions 

# ===== BENCHMARK RUN ===== #

Random.seed!(12345678)

# initial distribution
# [x_position, y_position, #_driver_mutations]
X₀ = [[0.,0.,1]]               # Initial tumor cell population (with driver mutation)
#X₀ = [[0.,0.,1]]

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
β = 1/5.
shift = 10
# Different types of fertility functions
# f = (x,y,t) -> 1/(1+exp(-α*(x+shift))) 
# f = (x,y,t) -> 1
f = (x,y,t) -> 1-exp(-β*(x^2))
# f = (x,y,t) -> exp(-β*x^2)
# f = (x,y,t) -> 1 - exp(-(β/abs(log(t+1)))*(x^2+y^2))
# h = (x,y,t) -> max(cos(pi*x/100)*cos(pi*y/100),0.2)

treatment_time = 1 # time in months
# prob_chemo = 0.8
# treatment_start = 5 # start treatment after 5 years

function f_chemo(x,y,t;treatment_time=treatment_time,k=k,treatment_start=treatment_start)
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
extinction = false  # Allow process to go extinct/condition on survival
Nrec = 0 # Count number of recurrences

###############
## Pull Data ##
###############
# Number of runs or patients to simulate
Nav = 50

# Runes benchmark() Nav number of runs and retunrns in array "data"
# data = [[X,Y,L,Nrec,Nm,Nt], ... ]
tick()
data = PullData(X₀,s,u,dx,dy,f_chemo,T,k,extinction,Nrec, Nav)
tock()

##################################
## SAVE DATA OF TO SCENARIO.JLD ##
##################################

## Save data and give file name

# ONLY FOR CREATING 
# save("scenarios_fertility.jld2",scenario,data)

# ADD new variables to save

# f = jldopen("scenarios.jld2", "r+")
# write(f, scenario, data)
# close(f)



## Loading data
# See what file contians

# load("scenarios.jld2")

# Access certain file

# data = load("scenarios.jld2")[scenario]



#############################################################################
## All data minipulations is done in the post_process_analysis_v23.jl file ##
#############################################################################







