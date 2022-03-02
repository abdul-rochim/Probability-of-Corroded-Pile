#include("mcs.jl")
include("chlorideIonConcentration.jl")
#using BenchmarkTools
import .Chloride_IC

function plot_chloride_choice()
    xf = XLSX.readxlsx("input/data.xlsx")
    sh = xf["rv_mcs"]
    return sh
end

function plot_chloride(params_object)
    sh = plot_chloride_choice()
    #choose for plotting
    lc_ON = sh["E54"]               #1 for plotting life cycle and 0 or others(except 1) for not plot life cycle
    c_coef_ON = sh["E55"]           #1 for plotting chloride coefficient and 0 or others(except 1) for not plotting
    c_c_ON = sh["E56"]              #1 for plotting chloride ion concentration and 0 or others(except 1) for not plotting
    alternatif_c_c_ON = sh["E57"]   #1 for plotting chloride ion concentration using Dc function and
                                    #2 for plotting chloride ion concentration using mean of Dc (based on data given)
                                    #3 for plotting chloride ion concentration using MCS, MCS analysis should be ON

    @unpack μx, handle, is_real_time_plotting, t_used, t_month, t_year, μCcr, cov_Ccr, h = Params_Data()
    if input_data()["E46"] == 2
        μCcr = exp(μCcr + 0.5 * log(cov_Ccr^2 + 1))    #conversion from μln(x) to μx
    end
    #definition of object for plotting chloride ion concentration
    #typealiases
    crt = Chloride_IC.Range_Time{Array{Float64,1}}
    clc = Chloride_IC.LifeCycle{Float64}
    cdc = Chloride_IC.ChlorideDiffusionCoefficient{Float64}
    cic = Chloride_IC.ChlorideIonConcentration{Float64}

    #Range Time and Life Cycle
    rt = crt(t_month, t_year, t_used)
    lc_m = clc(Chloride_IC.tMonth(rt)) #(rt.tMonth)
    lc_y = clc(Chloride_IC.tYear(rt)) #(rt.tYear)
    lc_u = clc(Chloride_IC.tUsed(rt))

    #Chloride Diffusion Coefficient
    #@btime ChlorideDiffusionCoefficient{Float64}(rt.tYear, lc_y.T, lc_y.RH/100)
    dc = cdc(Chloride_IC.tYear(rt), Chloride_IC.T(lc_y), Chloride_IC.RH(lc_y)/100)
    dc_ = cdc(Chloride_IC.tUsed(rt), Chloride_IC.T(lc_u), Chloride_IC.RH(lc_u)/100)

    #Chloride Ion Concentration
    ic = cic(dc_, μx, Chloride_IC.tYear(rt))
    ic_ = cic(dc_, μx, Chloride_IC.tUsed(rt))

    #plotting
    ################################################ Life Cycle ##################################################
    plt1 = plot([0],[0])
    if lc_ON == 1
        try
            p_lc_t = plot(Chloride_IC.tMonth(rt), Chloride_IC.T(lc_m); color=:red, w=:1, legend= :outertopright,
            xlabel = "t [month]", ylabel = "T [ᵒK]", label = "Fitted; T", xticks = 0:2:12)
            p_lc_rh = plot(Chloride_IC.tMonth(rt), Chloride_IC.RH(lc_m); c=:black, w=:1, legend=:outertopright,
            xlabel = "t [month]", ylabel = "RH [%]", label = "Fitted; RH", xticks = 0:2:12)
            plt1 = plot(p_lc_t, p_lc_rh)
        catch e
            println("Please, do not choose Live Plotting MCS!")
        end
    end
    ################################### Chloride Ion Concentration (Coefficient) #################################
    plt2 = plot([0],[0])
    if c_coef_ON == 1
        try
            p_cc_f1 = plot(t_year, dc.F1, xticks = 0:10:75, yticks = 0.7:0.05:1.0, color=:black, w=:1,
            legend=:outertopright, xlabel = "t [year]", ylabel = "F₂(t)", title = "Function F₁", label = nothing)
            p_cc_f2 = plot(Chloride_IC.T(lc_y), Chloride_IC.F2(dc), xticks = 280:5:300, yticks = 0.4:0.2:1.2,
            color=:black, w=:1, legend=:outertopright, xlabel = "T [ᵒK]", ylabel = "F₂(T)", title = "Function F₂",
            label = nothing)
            p_cc_f3 = plot(Chloride_IC.RH(lc_y), Chloride_IC.F3(dc), xticks = 60:5:80, yticks = 0.1:0.1:1.6,
            c=:black, w=:1, legend=:outertopright, xlabel = "RH [%]", ylabel = "F₃(RH)", title = "Function F₃",
            label = nothing)
            p_cc_dc = plot(Chloride_IC.tUsed(rt)/handle, Chloride_IC.Dc(dc_), c=:black, w=:0.25, legend=:outertopright,
            xlabel = "t [year]", ylabel = "Dc [m²/s]", title = "Diffusion Coefficient", label = nothing)
            plt2 = plot(p_cc_f1, p_cc_f2, p_cc_f3, p_cc_dc)
        catch e
            println("Please, do not choose Live Plotting MCS!")
        end
    end
    ##################################### Chloride Ion Concentration (C) #########################################
    plt3 = plot([0],[0])
    if c_c_ON == 1
        #cycle 75 years
        if alternatif_c_c_ON == 1
            try
                p_cic_atm = plot(t_used/handle, Chloride_IC.C_atm(ic_), c=:red, w=:1, legend=:outertopright,
                label = "Atmosperic")
                p_cic_splash = plot!(p_cic_atm, t_used/handle, Chloride_IC.C_splash(ic_), c=:blue, w=:1,
                legend=:outertopright, label = "Splash/Tidal")
                p_cic_sub = plot!(p_cic_splash, t_used/handle, Chloride_IC.C_sub(ic_), c=:black, w=:1,
                legend=:outertopright, label = "Submersion")
                plt3 = plot(p_cic_sub, title = "Chloride Content of Concrete", xlabel = "t [year]",
                ylabel = "C [kg/m³]")
            catch e
                println("Please, do not choose Live Plotting MCS!")
            end
        elseif alternatif_c_c_ON == 2
            try
                plot_Ccr = plot([0, length(t_used)/(handle*h)], [μCcr, μCcr], c=:magenta, w=:1, s=:dash, 
                label=nothing, annotations=[(length(t_used)/(handle*h*3/2), μCcr+0.15, 
                text("Critical chloride", :left, :magenta, 8))])
                p_cic_atm = plot!(plot_Ccr, t_used/handle, Chloride_IC.C_atm_dc_mean(ic_), c=:red, w=:1,
                legend=:outertopright, label = "Atmosperic")
                p_cic_splash = plot!(p_cic_atm, t_used/handle, Chloride_IC.C_splash_dc_mean(ic_), c=:blue, w=:1,
                legend=:outertopright, label = "Splash/Tidal")
                p_cic_sub = plot!(p_cic_splash, t_used/handle, Chloride_IC.C_sub_dc_mean(ic_), c=:black, w=:1,
                legend=:outertopright, label = "Submersion")
                plt3 = plot(p_cic_sub, title = "Chloride Content of Concrete", xlabel = "t [year]",
                ylabel = "C [kg/m³]")
            catch e
                println("Please, do not choose Live Plotting MCS!")
            end
        elseif alternatif_c_c_ON == 3
            if is_real_time_plotting == 0 
                mcs(params_object)
                obj_C_mcs_cycle_atm, obj_C_mcs_cycle_splash, obj_C_mcs_cycle_sub = params_object
            end
            try
                plot_Ccr = plot([0, length(t_used)/(handle*h)], [μCcr, μCcr], c=:magenta, w=:1, s=:dash, 
                label=nothing, annotations=[(length(t_used)/(handle*h*3/2), μCcr+0.15, 
                text("Critical chloride", :left, :magenta, 8))])                
                plt3 = plot(plot(plot(plot_Ccr, t_used/handle, obj_C_mcs_cycle_atm.C_mcs_cycle[3], 
                c=:red, w=:1, legend=:outertopright, label = "Atmosperic"),
                t_used/handle, obj_C_mcs_cycle_splash.C_mcs_cycle[3], c=:blue, w=:1,
                legend=:outertopright, label = "Splash/Tidal"),
                t_used/handle, obj_C_mcs_cycle_sub.C_mcs_cycle[3], c=:black, w=:1,
                legend=:outertopright, label = "Submersion", title = "Chloride Content of Concrete",
                xlabel = "t [year]", ylabel = "C [kg/m³]") #75 years
            catch e
                println("Please, do not choose Live Plotting MCS!")
            end
        else
            println("Choose the alternatif_c_c_ON for choice 1 or for choice 2!")
        end
    end

    if lc_ON == 1 && c_coef_ON == 0 && c_c_ON == 0
        plot(plt1)
    elseif lc_ON == 0 && c_coef_ON == 1 && c_c_ON == 0
        plot(plt2)
    elseif lc_ON == 0 && c_coef_ON == 0 && c_c_ON == 1
        plot(plt3)
    elseif lc_ON == 1 && c_coef_ON == 1 && c_c_ON == 0
        plot(plt1,plt2) 
    elseif lc_ON == 1 && c_coef_ON == 0 && c_c_ON == 1
        plot(plt1,plt3)
    elseif lc_ON == 0 && c_coef_ON == 1 && c_c_ON == 1
        plot(plt2,plt3)
    elseif lc_ON == 1 && c_coef_ON == 1 && c_c_ON == 1
        plot(plt1,plt2,plt3)
    else
        println("Please choose the choices for plotting-ON chloride ion concentration function!")
    end
end