include("mcs.jl")
include("plot_chloride_IC.jl")
include("plot_pdf.jl")
#=
_____________________________________________________________________________________________________________
 Program Name  : Analysis Probability of Corrosion using Monte Carlo Simulation                         
 Version       : v0                                                                                     
 © 𝐀𝐛𝐝𝐮𝐥 𝐑𝐨𝐜𝐡𝐢𝐦 (abdul.rochim.civeng@gmail.com)                                                                                                      
_____________________________________________________________________________________________________________
=#

function plot_mcs(params_object)
#=____________________________________  Plot probability of corrosion  ____________________________________=#
    @unpack handle, num_year, t_used = Params_Data()
    obj_C_mcs_cycle_atm, obj_C_mcs_cycle_splash, obj_C_mcs_cycle_sub = params_object
    mcs(params_object)  #mcs() invoked
    if is_real_time_plotting == 0
        p_10percent = plot([0.0, num_year], [0.1, 0.1], c=:magenta, w=:0.75, s=:dash, label=nothing)
        p_50percent = plot(p_10percent, [0.0, num_year], [0.5, 0.5], c=:magenta, w=:0.75, s=:dash, label=nothing)

        p_mcs_atm_75 = plot(p_50percent, t_used/handle, obj_C_mcs_cycle_atm.c_cp_mcs[1], c=:red, w=:1, 
        legend=:outertopright, label = "Atmosperic")

        p_mcs_splash_75 = plot!(p_mcs_atm_75, t_used/handle, obj_C_mcs_cycle_splash.c_cp_mcs[1],
        c=:blue, w=:1, legend=:outertopright, label = "Splash/Tidal")

        p_mcs_sub_75 = plot!(p_mcs_splash_75, t_used/handle, obj_C_mcs_cycle_sub.c_cp_mcs[1],
        c=:black, w=:1, legend=:outertopright, label = "Submersion")

        plot(p_mcs_sub_75, xlabel="t [year]", ylabel="Prob. of Corrosion",
        title = "Probability of Corrosion for Different Zones")
    end       
end #end of plot_mcs()

function main()
    @unpack plotting_choice = Params_Data()
    if plotting_choice == 1
        try
            plot_chloride(params_object_mcs)
        catch e
            println("Plot of Chloride content does not appears...")
        end
    elseif plotting_choice == 2
        try
            plot_pdf_(params__, params_object_mcs)    
        catch e
            println("\nYou enter wrong data.\nIt will not show the PDF.
You should NOT choose Live plotting!")
        end
    elseif plotting_choice == 3  
        plot_mcs(params_object_mcs)
    else
        println("The choices for plotting are 1 for plot chloride concentration, 2 for plot pdf and 3 for plot mcs!")
        @assert plotting_choice == 1 || plotting_choice == 2 || plotting_choice == 3 "You must enter 1 or 2 or 3!"
    end
end
main()
#=_________________________________________________________________________________________________________=#