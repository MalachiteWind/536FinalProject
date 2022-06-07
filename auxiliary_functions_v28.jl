#########################
## AUXILIARY FUNCTIONS ##
#########################
using Plots, Random, Distributions, Colors, LaTeXStrings, Printf, TickTock, JLD2, FileIO

function grab_values(network)
    
    output = values(network)
    n= length(output)
    X = zeros(n)
    Y = zeros(n)
    L = [] #Array(Int64,n)
    i=1
    
    for (x,y,l) in output
        X[i]=x
        Y[i]=y
        append!(L,[convert(Int64,l)])
        i+=1
    end
    
    return X, Y, L
    
end

function get_scatter_arrays(X,Y,L)
    
    # Ike changed this 6/5/22
    N = size(L)[1]
    M = maximum(L) # number of labels in data
    
    Xl = []
    Yl = []
    Ll = []
    for j=1:M
        Xl_j = []
        Yl_j = []
        for i=1:N
            if L[i] == j
                append!(Xl_j,[X[i]])
                append!(Yl_j,[Y[i]])
            end
        end
        if size(Xl_j)[1] > 0
            append!(Xl,[Xl_j])
            append!(Yl,[Yl_j])
            append!(Ll,[j])
        end
     end
    
    return Xl, Yl, Ll
    
end

function linspace(a,b,n)
    # this function creates an array starting at a, ending at b, with n points (so intervals 1/(n-1))
    # a and b should be floats, n should be an integer
    
    h = (b-a)/(n-1)
    x = zeros(n)
    x[1] = a
    x[n] = b
    
    for i=1:(n-2)
        x[i+1] = a+i*h
    end
    
    return x

end

function discretize_xy_nxny(X,Y,nx,ny,Xranges=[],Yranges=[])
    # Xranges = [xmin,xmax], Yranges = [ymin,ymax]
    
    if Xranges == []
        xmin = minimum(X)
        xmax = maximum(X)
    else
        xmin = Xranges[1]
        xmax = Xranges[2]
    end
        
    if Yranges == []
        ymin = minimum(Y)
        ymax = maximum(Y)
    else
        ymin = Yranges[1]
        ymax = Yranges[2]
    end

    # length scales
    Lx = xmax-xmin
    Ly = ymax-ymin
    
    # 20% margin in plots to not have the outer cells on the boundary
    cx = Lx/5.
    cy = Ly/5.
    
    x = linspace(xmin-cx, xmax+cx, nx)
    y = linspace(ymin-cy, ymax+cy, ny)
    
    return x,y,cx,cy
    
end

function discretize_xy_dxdy(X,Y,Δx,Δy,Xranges=[],Yranges=[])
        
    if Xranges == []
        xmin = minimum(X)
        xmax = maximum(X)
    else
        xmin = Xranges[1]
        xmax = Xranges[2]
    end
        
    if Yranges == []
        ymin = minimum(Y)
        ymax = maximum(Y)
    else
        ymin = Yranges[1]
        ymax = Yranges[2]
    end
    
    # length scales
    Lx = xmax-xmin
    Ly = ymax-ymin
    
    # 20% margin in plots to not have the outer cells on the boundary
    cx = Lx/5.
    cy = Ly/5.
    
    # number of grid points
    nx = convert(Int64, ceil((Lx+2*cx)/Δx))+1
    ny = convert(Int64, ceil((Ly+2*cy)/Δy))+1
    
    x = linspace(xmin-cx, xmax+cx, nx)
    y = linspace(ymin-cy, ymax+cy, ny)
    
    return x,y,cx,cy
    
end

function data_binning_nxny(X,Y,L,nx,ny,Xranges=[],Yranges=[])

    Xl, Yl, Ll = get_scatter_arrays(X,Y,L)
    
    M = size(Ll)[1]
    
    # length ranges for cells 
    x,y,cx,cy = discretize_xy_nxny(X,Y,nx,ny,Xranges,Yranges)
    xmin = x[1]
    xmax = x[end]
    ymin = y[1]
    ymax = y[end]

    Δx = (xmax-xmin)/(nx-1)
    Δy = (ymax-ymin)/(ny-1)

    XYbin = []
    Xinfo = []
    Yinfo = []
    
    for j=1:M
        
        XYbin_j = zeros(nx-1,ny-1)
        N = size(Xl[j])[1]
        
        for i=1:N
            
            k=1
            l=1
            xi = X[i]
            yi = Y[i]
            
            while (!(xi > xmin+(k-1)*Δx && xi <= xmin+k*Δx) && k < nx)
                k = k+1
            end
            
            while (!(yi > ymin+(l-1)*Δy && yi <= ymin+l*Δy) && l < ny)
                l = l+1
            end

            if k > nx-1
                println("ERROR [data_binning_dxdy]: while loop for xi exceeds indices, k set to nx-1")
                k = nx-1
            end
            
            if l > ny-1
                println("ERROR [data_binning_dxdy]: while loop for yi exceeds indices, l set to ny-1")
                l = ny-1
            end
            
            XYbin_j[k,l] = XYbin_j[k,l] + 1
            
        end
        
        append!(XYbin,[XYbin_j])
        append!(Xinfo,[[xmin,xmax,Δx]])
        append!(Yinfo,[[ymin,ymax,Δy]])
        
    end

    return XYbin,Xinfo,Yinfo,Ll

end

function data_binning_dxdy(X,Y,L,Δx,Δy,Xranges=[],Yranges=[])
    
    Xl, Yl, Ll = get_scatter_arrays(X,Y,L)

    M = size(Ll)[1]
    
    # length ranges for cells 
    x,y,cx,cy = discretize_xy_dxdy(X,Y,Δx,Δy,Xranges,Yranges)
    xmin = x[1]
    xmax = x[end]
    ymin = y[1]
    ymax = y[end]

    nx = size(x)[1]
    ny = size(y)[1]
    
    XYbin = []
    Xinfo = []
    Yinfo = []
    
    for j=1:M
        
        XYbin_j = zeros(nx-1,ny-1)
        N = size(Xl[j])[1]
        
        for i=1:N
            
            k=1
            l=1
            xi = X[i]
            yi = Y[i]
            
            while (!(xi > xmin+(k-1)*Δx && xi < xmin+k*Δx) && k < nx)
                k = k+1
            end
            
            while (!(yi > ymin+(l-1)*Δy && yi < ymin+l*Δy) && l < ny)
                l = l+1
            end
            
            if k > nx-1
                println("ERROR [data_binning_dxdy]: while loop for xi exceeds indices, k set to nx-1")
                k = nx-1
            end
            
            if l > ny-1
                println("ERROR [data_binning_dxdy]: while loop for yi exceeds indices, l set to ny-1")
                l = ny-1
            end
            
            XYbin_j[k,l] = XYbin_j[k,l] + 1
            
        end
        
        append!(XYbin,[XYbin_j])
        append!(Xinfo,[[xmin,xmax,Δx]])
        append!(Yinfo,[[ymin,ymax,Δy]])
        
    end

    return XYbin,Xinfo,Yinfo,Ll

end

function analysis_plots_nxny(X,Y,L,nx,ny,f,T,Xranges=[],Yranges=[],fn="")
    # this function creates a plots to analyze corresponding to the output network of, for example, the benchmark run
    # one plot given is a scatter plot of the types of cells
    # another plot given is a probability distribution for the vertility
    
    Xl, Yl, Ll = get_scatter_arrays(X,Y,L)

    #############################
    # ===== SCATTER PLOT ==== #
    #############################
    
    r = scatter(Xl[1],Yl[1], size = (720, 400), color=Ll[1], label="type "*string(Ll[1]))
    
    if size(Ll)[1] > 1
        
        for i in Ll[2,end]
            
            scatter!(Xl[i],Yl[i], size = (720, 400), color=i, label="type "*string(i))
            
        end
        
    end

    if fn == 0
        display(r)
    else
        png(r,"spatial_analysis_plot_S_"*fn)
    end
    
    ######################################################
    # ===== SCATTER PLOT + VERTILITY PROB. FUNCTION ==== #
    ######################################################

    xplot,yplot,cx,cy = discretize_xy_nxny(X,Y,nx,ny,Xranges,Yranges)
    
    function ff(x,y)
        f(x,y,T)
    end
    p = plot(contourf(xplot, yplot, ff; levels = linspace(0,1,20)),size = (800, 400))
    
    # ===== SCATTER PLOT ===== #
    

    scatter!(Xl[1],Yl[1], size = (800, 400), color=Ll[1], label="type "*string(Ll[1]))
    
    if size(Ll)[1] > 1
        
        for i in Ll[2,end]
            
            scatter!(Xl[i],Yl[i], size = (800, 400), color=i, label="type "*string(i))
            
        end
        
    end
    
    xlims!((minimum(X)-cx,maximum(X)+cx))
    ylims!((minimum(Y)-cy,maximum(Y)+cy))

    if fn == 0
        display(p)
    else
        png(p,"spatial_analysis_plot_SV_"*fn)
    end
    
    
        
    ######################################################
    # ===== VERTILITY PROBABILITY DISTRIBUTION PLOT ==== #
    ######################################################
    
    # ===== VERTILITY PROBABILITY DISTRIBUTION PLOT ===== #
    
    q = plot(contourf(xplot, yplot, ff; levels = linspace(0,1,20)),size = (800, 400))
        
    if fn == 0
        display(q)
    else
        png(q,"spatial_analysis_plot_V_"*fn)
    end
        
end

#######################
## Ike added funcion ##
#######################

function heatmap_grid(XY_bin,Xinfo,Yinfo,i) # To create the appropriate tick marsk for make_heat_map()
    nx = size(XY_bin[i][:,1])[1]
    ny = size(XY_bin[i][1,:])[1]

    xmin = Xinfo[i][1]
    xmax = Xinfo[i][2]

    ymin = Yinfo[i][1]
    ymax = Yinfo[i][2]

    x_space = (xmax-xmin)/nx
    y_space = (ymax-ymin)/ny

    x_axis = 0:nx
    y_axis = 0:ny

    x_axis = x_axis*x_space
    y_axis = y_axis*y_space

    x_axis = x_axis .+xmin
    y_axis = y_axis .+ymin

    return x_axis, y_axis
end


function make_heat_map(XY_bin,Xinfo,Yinfo,Ll,fn="")

    #######################
    ## Ike added funcion ##
    #######################

    x1,y1 = heatmap_grid(XY_bin,Xinfo,Yinfo,1)    

    s = heatmap(x1,y1,transpose(XY_bin[1]),xlabel=latexstring(@sprintf("x~\\mathrm{(in~terms~of~} \\Delta x = ")*string(round(Xinfo[1][3],digits=3))*")"), 
        ylabel=latexstring(@sprintf("y~\\mathrm{(in~terms~of~} \\Delta y = ")*string(round(Yinfo[1][3],digits=3))*")"),
        title=latexstring(@sprintf("\\mathrm{type~}")*string(Ll[1])))
    # s = transpose(s)
    if fn == 0
        display(s)
    else
        savefig(s,"spatial_analysis_plot_H_new"*string(Ll[1])*"_"*fn)
    end
        
    if size(Ll)[1] > 1
        
        for i in Ll[2,end]
            #######################
            ## Ike added funcion ##
            #######################

            xi,yi = heatmap_grid(XY_bin,Xinfo,Yinfo,i) 

            s = heatmap(xi,yi,transpose(XY_bin[i]),xlabel=latexstring(@sprintf("x~\\mathrm{(in~terms~of~} \\Delta x = ")*string(round(Xinfo[i][3],digits=3))*")"), 
            ylabel=latexstring(@sprintf("y~\\mathrm{(in~terms~of~} \\Delta y = ")*string(round(Yinfo[i][3],digits=3))*")"),
            title=latexstring(@sprintf("\\mathrm{type~}")*string(Ll[i])))
            
            if fn == 0
                display(s)
            else
                savefig(s,"spatial_analysis_plot_H_"*string(Ll[i])*"_"*fn)
            end
                
        end
        
    end
    
end

function heat_map_nxny(X,Y,L,nx,ny,Xranges=[],Yranges=[],fn="")

    XY_bin,Xinfo,Yinfo,Ll = data_binning_nxny(X,Y,L,nx,ny,Xranges,Yranges)

    make_heat_map(XY_bin,Xinfo,Yinfo,Ll,fn)
    
end

function heat_map_dxdy(X,Y,L,Δx,Δy,Xranges=[],Yranges=[],fn="")

    XY_bin,Xinfo,Yinfo,Ll = data_binning_dxdy(X,Y,L,Δx,Δy,Xranges,Yranges)
    
    make_heat_map(XY_bin,Xinfo,Yinfo,Ll,fn)
    
end

function analysis_plots_time(Nm,Nt,X₀,fn="")

    N = size(Nm)[2]
    Ntime = linspace(0,T,N)

    #     println(size(Nm))
    #     println(size(Nt))

    m1=size(Nm)[1]
    p = scatter(Ntime,Nm[1,1:N], size = (720, 400), color=1, label="mother cell "*string([X₀[1][1],X₀[1][2]]), legend=:topleft)
    for i=2:m1
        scatter!(Ntime,Nm[i,1:N], size = (720, 400), color=i, label="mother cell "*string([X₀[i][1],X₀[i][2]]), legend=:topleft)
    end
    xlabel!(latexstring(@sprintf("t")))
    ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
 
    if fn == 0
        display(p)
    else
        savefig(p,"spatial_analysis_plot_T1_"*fn)
    end
    
    Ndriver=size(Nt)[1]



    q = scatter(Ntime,Nt[1,1:N], size = (720, 400), color=1, label="type "*string(1), legend=:topleft)
    for i=2:Ndriver
        ####################
        ## IKE ADDED THIS ##
        ####################
        # scatter!(Ntime,Nt[i,1:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft)
        if any(Nt[i,:].>0.)
            vec = Nt[i,:].>0.
            index = findmax(vec)[2]
            scatter!(Ntime[index:end],Nt[i,index:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft)
        end
    end
    xlabel!(latexstring(@sprintf("t")))
    ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
    
    if fn == 0
        display(q)
    else
        savefig(q,"spatial_analysis_plot_T2_"*fn)
    end
    
    # ===== LOG VERSIONS ===== #
    
    #     ind_m = (Nm[1,1:N] .> 0)
    #     ind_t = (Nt[1,1:N] .> 0)
    
    Nm_log = Nm
    Nt_log = Nt
    for i=1:m1
        for j=1:N
            if Nm_log[i,j] == 0
                Nm_log[i,j] = NaN
            end
        end
    end
    for i=1:Ndriver
        for j=1:N
            if Nt_log[i,j] == 0
                Nt_log[i,j] = NaN
            end
        end
    end

    
    
    p_log = scatter(Ntime,Nm_log[1,2:N], size = (720, 400), color=1, label="mother cell "*string([X₀[1][1],X₀[1][2]]), legend=:topleft, yaxis=:log)
    for i=2:m1
        scatter!(Ntime,Nm_log[i,2:N], size = (720, 400), color=i, label="mother cell "*string([X₀[i][1],X₀[i][2]]), legend=:topleft, yaxis=:log)
    end
    xlabel!(latexstring(@sprintf("t")))
    ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
 
    if fn == 0
        display(p_log)
    else
        savefig(p_log,"spatial_analysis_plot_T1_log_"*fn)
    end
    
    q_log = scatter(Ntime,Nt_log[1,1:N], size = (720, 400), color=1, label="type "*string(1), legend=:topleft, yaxis=:log)
    for i=2:Ndriver
        ####################
        ## IKE ADDED THIS ##
        ####################
        # scatter!(Ntime,Nt_log[i,1:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft, yaxis=:log)
        if any(Nt_log[i,:].>0.)
            vec = Nt_log[i,:].>0.
            index = findmax(vec)[2]
            scatter!(Ntime[index:end],Nt_log[i,index:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft, yaxis=:log)
        end
    end
    xlabel!(latexstring(@sprintf("t")))
    ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
    
    if fn == 0
        display(q_log)
    else
        savefig(q_log,"spatial_analysis_plot_T2_log_"*fn)
    end
        
end


### OLD FUNCTION ##
function average_run_nxny(Nav,nx,ny,Xranges,Yranges,X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
    
    N = convert(Int64, ceil(T*365.25/k))      # Number of instances 
    Nav = float(Nav)
    
    # initialization
    X,Y,L,Nrec,Nm,Nt = benchmark(X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
    Xl, Yl, Ll = get_scatter_arrays(X,Y,L)
    XY_bin,Xinfo,Yinfo = data_binning_nxny(X,Y,L,nx,ny,Xranges,Yranges)
    Nlabels = size(XY_bin)[1]
    for j=1:Nlabels
        XY_bin[j] = XY_bin[j]/Nav
    end
    Nm = Nm/Nav
    Nt = Nt/Nav
    
    for i=2:Nav
        X,Y,L,Nrec,Nm_i,Nt_i = benchmark(X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
        Xl, Yl, Ll = get_scatter_arrays(X,Y,L)
        XY_bin_i,_,_ = data_binning_nxny(X,Y,L,nx,ny,Xranges,Yranges)
        for j=1:Nlabels
            XY_bin[j] = XY_bin[j] + XY_bin_i[j]/Nav
        end
        Nm = Nm + Nm_i/Nav
        Nt = Nt + Nt_i/Nav
    end
    
    return XY_bin,Xinfo,Yinfo,Ll,Nm,Nt
    
end


#####################
## BENCHMARK MODEL ##
#####################

function benchmark(X₀,s,u,dx,dy,f,T,k,extinction=false,Nrec=0)
    
    # ===== INITIALIZATION ===== #

    N = convert(Int64, ceil(T*365.25/k))      # Number of instances 
    m1 = size(X₀)[1]                       # Starting populations
    #    Nrec = 0                              # number of recurrent calls
    
    # initialize driver mutation probability array p 
    # array p contains Ndriver arrays of the form [b,d,m], 
    # where b is the birth rate for that mutation, d the death rate and m the mutation rate to the next type
    Ndriver = size(u)[1]+1
    b = zeros(Ndriver)
    d = zeros(Ndriver)
    p = []
    for i=1:(Ndriver-1)
        d[i] = 0.5*(1.0-s)^i#/2.
        b[i] = 1.0-d[i] #1.0-((1-s)^i/2.)
        append!(p,[[b[i]*(1.0-u[i]),d[i],b[i]*u[i]]])
    end
    d[Ndriver] = (1-s)^Ndriver/2.
    b[Ndriver] = 1-d[Ndriver]
    append!(p,[[b[Ndriver],d[Ndriver],0]])

    # ===== MAIN ROUTINE ===== #
    
    # initialize dictorary with X_0 data
    network = Dict()
    for i=1:m1
        network[string(i)*".0"] = X₀[i]
    end
    
    Nm = zeros(m1,N+1)
    Nt = zeros(Ndriver,N+1)
    
    for j=1:m1
        Nm[j,1] = 1
    end
    
    for j=1:Ndriver
        for i=1:m1
            if X₀[i][3] == j
                Nt[j,1] = Nt[j,1] + 1
            end
        end
    end
               
    tt=0
    
    for j=1:N # Iterate over time
        
        network_update = Dict() # network to put cells born at new time step into
        tt = tt+k
        
        for (key,value) in network # Iterate over each cell in network
            
            tk=convert(Int64,value[3])             # (driver mutation) type

            #################################
            ## Check fertility probability ##
            #################################
            ###########################
            ### IKE ADDED THIS PART ###
            ###########################
            fert_prob = f(value[1],value[2],tt)
            d_temp = 1/2*(1-s)^value[3]
            birth = (1-d_temp)*fert_prob
            death = 1- birth
            
            prob_vector = [birth*(1-u[convert(Int64,value[3])]),death,birth*u[convert(Int64,value[3])]]
            
            ###### I changed p[tk] to prob_vector
            bk,dk,mk = rand(Multinomial(1, prob_vector)) # draw whether birth, death or mutation occurs for the cell regarded (only one of them is 1, obviously)
        #             println(bk,dk,mk)
            
            # Update cell population
            if bk == true # If birth or mutation occurs, update network_update accordingly
                
                x = rand(dx(value[1]),1)[1]         # x-coordinate 
                y = rand(dy(value[2]),1)[1]         # y-coordinate
                
                new_key = key*"."*"b $j"     # cell name for dictionary
                
                # check whether new cell survives on the new spot (compare to vertility probability)
                pk = rand(Uniform(0., 1.))
                if f(x,y,tt) > pk
                    network_update[new_key] = [x,y,tk] # Add new cell to network
                end
                
             #                 println(split(new_key,"."))
             #                 println(split(new_key,".")[1])
                im = split(new_key,".")[1]
                im = parse(Int64, string(split(new_key,".")[1])) # mother cell index
             #                 println(im+2)
                Nm[im,j+1] = Nm[im,j+1] + 1
                Nt[tk,j+1] = Nt[tk,j+1] + 1

            elseif mk == true
                delete!(network,key)         # delete old network entry
                new_key = key*"."*"m $j"     # cell name for dictionary
                network_update[new_key] = [value[1],value[2],tk+mk]
                
                im = split(new_key,".")[1]
                im = parse(Int64, string(split(new_key,".")[1])) # mother cell index
                Nm[im,j+1] = Nm[im,j+1] + 1
                Nt[tk+mk,j+1] = Nt[tk+mk,j+1] + 1
                
            else # Otherwise death occurs
                delete!(network,key) # Death occurs
            end

            # Check if population has gone extinct
            # If extinction has occured and extinction is set to false then 
            # rerun simuluation until surviving process occurs.
            if ((isempty(network)==true) && (extinction==false))
                Nrec = Nrec + 1 
                return benchmark(X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
            end
            
        end
        
        network = merge(network,network_update)
        
    end
    
    X,Y,L = grab_values(network)
    
    return X,Y,L,Nrec,Nm,Nt

end


##########################
## ACQUIRE PATIENT DATA ##
##########################

# Runes benchmark() Nav number of runs and retunrns in array "data"
# data = [[X,Y,L,Nrec,Nm,Nt], ... ]
function PullData(X₀,s,u,dx,dy,f,T,k,extinction,Nrec, Nav)
    data = Array{Any}(nothing, Nav)
    for i=1:Nav
        X,Y,L,Nrec,Nm,Nt = benchmark(X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
        data[i]=[X,Y,L,Nrec, Nm, Nt]
        println("number of recurrent calls in run = "*string(Nrec))
    end
    data
end

function PullDataBenchmarkTherapy(X0s, s, u, dx,dy,f,T,k,extinction, Nrec)
    Nav = size(X0s)[1]
    data = Array{Any}(nothing, Nav)
    for i = 1:Nav
        X0_i = X0s[i]
        X,Y,L,Nrec,Nm,Nt = benchmark(X0_i,s,u,dx,dy,f,T,k,extinction,Nrec)
        data[i] = [X,Y,L,Nrec,Nm,Nt]
        println("number of recurrent calls in run = $Nrec")
    end
    return data, Nav
end

# Grab patient i data
function Grab_Patient_data(data,i)
    X = data[i][1]
    Y = data[i][2]
    L = data[i][3]
    Nrec = data[i][4]
    Nm = data[i][5]
    Nt = data[i][5]
    return X, Y, L, Nrec, Nm, Nt
end


##########################
## AVERAGE PATIENT DATA ##
##########################

####################
## UPDATED BY IKE ##
####################

function average_run_nxny2(Nav,nx,ny,Xranges,Yranges,data)
    # Rewrite average run so that we don't run into edge case where first patient has less type i cells
    N = convert(Int64, ceil(T*365.25/k))      # Number of instances 
    # Nav = float(Nav)

    XY_bins_plural = Array{Any}(nothing, Nav)       # To store all XY_bins
    Nlabels_plural = Array{Any}(nothing, Nav)       # To store number of driver mutations per patient
    Nm_plural = Array{Any}(nothing, Nav)            # To store all Nm
    Nt_plural = Array{Any}(nothing, Nav)            # To store all Nt
 
    Ll = 0
    Xinfo=0
    Yinfo=0

    Nlabels_max = 0
    Nlabels_max_index = 0
    for i=1:Nav
        X_i,Y_i,L_i,Nrec_i,Nm_i,Nt_i = data[i] # Grab data from patient
        _,_,Ll_i = get_scatter_arrays(X_i,Y_i,L_i)
        XY_bin_i,Xinfo_i,Yinfo_i = data_binning_nxny(X_i,Y_i,L_i,nx,ny,Xranges,Yranges) # Grab binning infor

        XY_bins_plural[i] = XY_bin_i    # Store values
        Nm_plural[i] = Nm_i
        Nt_plural[i] = Nt_i

        temp_Nlables = size(XY_bin_i)[1] # Pick a patient that had the most driver mutations
        Nlabels_plural[i] = temp_Nlables


        if temp_Nlables > Nlabels_max 
            Nlabels_max = temp_Nlables      # Max driver mutations of patients
            Nlabels_max_index = i           # First patient that acquire max driver mutations of group
            Ll = Ll_i
            Xinfo = Xinfo_i
            Yinfo = Yinfo_i
        end
    end
    
    ### calculate average 

    for i=1:Nav
        if i != Nlabels_max_index
            for j=Nlabels_plural[i]
                XY_bins_plural[Nlabels_max_index][j] += XY_bins_plural[i][j]
            end
        end
    end

  
    XY_bin = XY_bins_plural[Nlabels_max_index]/Nav
    Nm = sum(Nm_plural)/Nav
    Nt = sum(Nt_plural)/Nav
    Nm = Nm/Nav
    #     Nt = Nt/Nav
    
    #     # initialization
    #     X,Y,L,Nrec,Nm,Nt = data[1] #benchmark(X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
    #     _, _, Ll = get_scatter_arrays(X,Y,L)
    #     XY_bin,Xinfo,Yinfo = data_binning_nxny(X,Y,L,nx,ny,Xranges,Yranges)
    #     Nlabels = size(XY_bin)[1]
    #     for j=1:Nlabels
    #         XY_bin[j] = XY_bin[j]/Nav
    #     end
    #     Nm = Nm/Nav
    #     Nt = Nt/Nav
        
    #     for i=2:Nav
    #         println(i)
    #         X,Y,L,Nrec,Nm_i,Nt_i = data[convert(Int64,i)] #benchmark(X₀,s,u,dx,dy,f,T,k,extinction,Nrec)
    #         # println(size(data[convert(Int64,i)])[1])
    #  #       _, _, Ll = get_scatter_arrays(X,Y,L)
    #         XY_bin_i,_,_ = data_binning_nxny(X,Y,L,nx,ny,Xranges,Yranges)
    #         # println(XY_bin_i)
    #         Nlabels = size(XY_bin_i)[1]
    #         for j=1:Nlabels
    #             # println(j)
    #             XY_bin[j] = XY_bin[j] + XY_bin_i[j]/Nav
    #         end
    #         Nm = Nm + Nm_i/Nav
    #         Nt = Nt + Nt_i/Nav
    #     end
    
    return XY_bin,Xinfo,Yinfo,Ll,Nm,Nt
end

##########################
## THERAPY TO BENCHMARK ##
##########################

# for converting data from bench mark case to initial conditions. 
# So we can rerun with chemotherapy at some frequence.

function data_to_X0s(data)
    Nav = size(data)[1]
    X0s = Array{Any}(nothing, Nav)
    for i=1:Nav
        pop_size = size(data[i][1])[1]
        X0_i = Array{Any}(nothing, pop_size)
        for j=1:pop_size
            xj = data[i][1][j]
            yj = data[i][2][j]
            lj = convert(Int64,data[i][3][j])
            X0_i[j] = [xj,yj,lj]
        end
        X0s[i]=X0_i
    end
    return X0s
end

function analysis_plots_time_therapy(Nm,Nt,X₀,fn="")

    N = size(Nm)[2]
    Ntime = linspace(0,T,N)

    #     println(size(Nm))
    #     println(size(Nt))

    m1=size(Nm)[1]
    # p = scatter(Ntime,Nm[1,1:N], size = (720, 400), color=1, label="mother cell "*string([X₀[1][1],X₀[1][2]]), legend=:topleft)
    # for i=2:m1
    #     scatter!(Ntime,Nm[i,1:N], size = (720, 400), color=i, label="mother cell "*string([X₀[i][1],X₀[i][2]]), legend=:topleft)
    # end
    # xlabel!(latexstring(@sprintf("t")))
    # ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
 
    # if fn == 0
    #     display(p)
    # else
    #     savefig(p,"spatial_analysis_plot_T1_"*fn)
    # end
    
    Ndriver=size(Nt)[1]



    q = scatter(Ntime,Nt[1,1:N], size = (720, 400), color=1, label="type "*string(1), legend=:topleft)
    for i=2:Ndriver
        ####################
        ## IKE ADDED THIS ##
        ####################
        # scatter!(Ntime,Nt[i,1:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft)
        if any(Nt[i,:].>0.)
            vec = Nt[i,:].>0.
            index = findmax(vec)[2]
            scatter!(Ntime[index:end],Nt[i,index:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft)
        end
    end
    xlabel!(latexstring(@sprintf("t")))
    ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
    
    if fn == 0
        display(q)
    else
        savefig(q,"spatial_analysis_plot_T2_"*fn)
    end
    
    # ===== LOG VERSIONS ===== #
    
    #     ind_m = (Nm[1,1:N] .> 0)
    #     ind_t = (Nt[1,1:N] .> 0)
    
    Nm_log = Nm
    Nt_log = Nt
    for i=1:m1
        for j=1:N
            if Nm_log[i,j] == 0
                Nm_log[i,j] = NaN
            end
        end
    end
    for i=1:Ndriver
        for j=1:N
            if Nt_log[i,j] == 0
                Nt_log[i,j] = NaN
            end
        end
    end

    
    
    # p_log = scatter(Ntime,Nm_log[1,2:N], size = (720, 400), color=1, label="mother cell "*string([X₀[1][1],X₀[1][2]]), legend=:topleft, yaxis=:log)
    # for i=2:m1
    #     scatter!(Ntime,Nm_log[i,2:N], size = (720, 400), color=i, label="mother cell "*string([X₀[i][1],X₀[i][2]]), legend=:topleft, yaxis=:log)
    # end
    # xlabel!(latexstring(@sprintf("t")))
    # ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
 
    # if fn == 0
    #     display(p_log)
    # else
    #     savefig(p_log,"spatial_analysis_plot_T1_log_"*fn)
    # end
    
    q_log = scatter(Ntime,Nt_log[1,1:N], size = (720, 400), color=1, label="type "*string(1), legend=:topleft, yaxis=:log)
    for i=2:Ndriver
        ####################
        ## IKE ADDED THIS ##
        ####################
        # scatter!(Ntime,Nt_log[i,1:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft, yaxis=:log)
        if any(Nt_log[i,:].>0.)
            vec = Nt_log[i,:].>0.
            index = findmax(vec)[2]
            scatter!(Ntime[index:end],Nt_log[i,index:N], size = (720, 400), color=i, label="type "*string(i), legend=:topleft, yaxis=:log)
        end
    end
    xlabel!(latexstring(@sprintf("t")))
    ylabel!(latexstring(@sprintf("N_{\\textrm{cells}}")))
    
    if fn == 0
        display(q_log)
    else
        savefig(q_log,"spatial_analysis_plot_T2_log_"*fn)
    end
        
end

function post_process_chemo(data,X0s,nx,ny,f_chemo,savefig)
    Nav = size(data)[1]
    success_count = Nav
    for i=1:Nav
        if isempty(data[i][1]) == false
            if savefig==true
                fn = "patient$i" 
            else
                fn =0
            end
            success_count -=1    # Reduce the successes

            X, Y, L, Nrec, Nm, Nt = Grab_Patient_data(data,i)
            X₀ = X0s[i]

            println("patient$i")
            heat_map_nxny(X,Y,L,nx,ny,[],[],fn)
            analysis_plots_nxny(X,Y,L,nx,ny,f_chemo,T,[],[],fn) # (3 plots) Return Scatter & fertility, just scatter, just fertility
            # analysis_plots_time_therapy(Nm,Nt,X₀,fn) # (2 plot)
        end
    end
    println("$success_count of $Nav patients tumour free.")
end

