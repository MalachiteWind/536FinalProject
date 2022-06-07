
using Plots, Random, Distributions, Colors, LaTeXStrings, Printf, TickTock, JLD2

include("auxiliary_functions_v28.jl") # Pull in functions

###########################
## POST~PROCESS~ANALYSIS ##
###########################

## Load data from scenarios.jld ##
# scenario = "scenario_isotropic"
# f = jldopen("scenarios_spatial.jld2","r")
# data = f[scenario]
# close(f)
# data = load("scenarios.jld2")[scenario]



# Grab specific data from paintient i ∈ [1,Nav]
i=6
X, Y, L, Nrec, Nm, Nt = Grab_Patient_data(data,i)

##############
## PLOTTING ##
##############

# fn="patient$i" # For saving images
# fn="average_tumour"
fn=0

nx = 101 # Number of x bins 
ny = 101 # Number of y bins

###############
## Patient i ##
###############

# heat_map_nxny(X,Y,L,nx,ny,[],[],fn)
# analysis_plots_time(Nm,Nt,X₀,fn) # (2 plot)
# analysis_plots_nxny(X,Y,L,nx,ny,f_chemo,T,[],[],fn) # (3 plots) Return Scatter & fertility, just scatter, just fertility


###############
## Average ##
###############
# Xranges=[-100.,100.]   # Range to be determined by surviving colonies (Look at patient data to get idea)
# Yranges=[-100.,100.]
# XY_bin_av,Xinfo,Yinfo,Ll,Nm_av,Nt_av = average_run_nxny2(Nav,nx,ny,Xranges,Yranges,data)
# make_heat_map(XY_bin_av,Xinfo,Yinfo,Ll,fn) # (1 plot) Average heatmap
# analysis_plots_time(Nm_av,Nt_av,X₀,fn) # (2 plots) Mother cell linearage and number of k types cells






# Δx = 3# bin width in x direction
# Δy = 5 # bin width in y direction
# heat_map_dxdy(X,Y,L,Δx,Δy,[],[],fn) # (1 plot) Heatmap (binned by fixed interval length Δx and Δy)
