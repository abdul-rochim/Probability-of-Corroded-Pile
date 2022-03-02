#include("mcs.jl")
include("pdf_cdf.jl")

function plot_pdf_choice()
    xf = XLSX.readxlsx("input/data.xlsx")
    sh = xf["rv_mcs"]
    return sh
end

#choose for plotting
params__ = (
    pdf_C_atm_ON = plot_pdf_choice()["E62"]                #1 for plotting pdf atmospheric and 0 or others(except 1) for not plotting
    ,pdf_C_splash_ON = plot_pdf_choice()["E63"]            #1 for plotting pdf splash/ tidal and 0 or others(except 1) for not plotting
    ,pdf_C_sub_ON = plot_pdf_choice()["E64"]               #1 for plotting pdf submersion and 0 or others(except 1) for not plotting
)

function plot_pdf_(params, params_object)
    pdf_C_atm_ON, pdf_C_splash_ON, pdf_C_sub_ON = params
    obj_C_mcs_cycle_atm, obj_C_mcs_cycle_splash, obj_C_mcs_cycle_sub = params_object
    @unpack handle, h, num_data, μCcr, t_used = Params_Data()
    @unpack σCcr = params_mu_sig

    myLabel = String[]
    push!(myLabel, "25 years")
    push!(myLabel, "50 years")
    push!(myLabel, "75 years")

    # Probability Distributions Function
    μC_mcs_atm = zeros(num_data)         
    σC_mcs_atm = zeros(num_data) 
    
    μC_mcs_splash = zeros(num_data)
    σC_mcs_splash = zeros(num_data)
    
    μC_mcs_sub = zeros(num_data)
    σC_mcs_sub = zeros(num_data)
    
    plot_atm_25 = plot([0],[0])
    plot_atm_25_50 = plot([0],[0])
    plot_atm_25_50_75 = plot([0],[0])
    plot_splash_25 = plot([0],[0])
    plot_splash_25_50 = plot([0],[0])
    plot_splash_25_50_75 = plot([0],[0])
    plot_sub_25 = plot([0],[0])
    plot_sub_25_50 = plot([0],[0])
    plot_sub_25_50_75 = plot([0],[0])

    plot_Ccr_1 = plot(t_used/(handle*h), pdf_used_Ccr(μCcr, σCcr, t_used/(handle*h)), color=:red, w=1, 
    legend= :outertopright, label="Critical")
    plot_Ccr_2 = plot(t_used/(handle*h)*2, pdf_used_Ccr(μCcr, σCcr, t_used/(handle*h)*2), color=:red, w=1, 
    legend= :outertopright, label="Critical")
    plot_Ccr_3 = plot(t_used/(handle*h), pdf_used_Ccr(μCcr, σCcr, t_used/(handle*h)), color=:red, w=1, 
    legend= :outertopright, label="Critical")
    mcs(params_object)
    for i in 1 : num_data
        #atmospheric
        if pdf_C_atm_ON == 1
            if plot_pdf_choice()["E42"] == 1    #1 for pdf Normal and 2 for pdf LogNormal
                μC_mcs_atm[i] = mean(obj_C_mcs_cycle_atm.C_mcs_cycle[i])
                σC_mcs_atm[i] = std(obj_C_mcs_cycle_atm.C_mcs_cycle[i])
            elseif plot_pdf_choice()["E42"] == 2
                Vx = std(obj_C_mcs_cycle_atm.C_mcs_cycle[i])/ mean(obj_C_mcs_cycle_atm.C_mcs_cycle[i])
                μC_mcs_atm[i] = log(mean(obj_C_mcs_cycle_atm.C_mcs_cycle[i])) - 0.5*log(Vx^2 + 1)
                σC_mcs_atm[i] = √(log(Vx^2 + 1))
            end

            if i == 1
                plot_atm_25 = plot!(plot_Ccr_1, t_used/(handle*h), pdf_used_C(μC_mcs_atm[i], σC_mcs_atm[i], t_used/(handle*h)), 
                xticks =0:1:10, color=:black, w=1, legend= :outertopright, label=myLabel[i])
            end
            if i == 2
                plot_atm_25_50 = plot!(plot_atm_25, t_used/(handle*h), pdf_used_C(μC_mcs_atm[i], σC_mcs_atm[i], t_used/(handle*h)), 
                xticks =0:1:10, color=:blue, w=1, legend= :outertopright, label=myLabel[i])
            end
            if i == 3
                plot_atm_25_50_75 = plot!(plot_atm_25_50, t_used/(handle*h), pdf_used_C(μC_mcs_atm[i], σC_mcs_atm[i], t_used/(handle*h)), 
                xticks =0:1:10, xlabel="C [kg/m³]", ylabel="PDF", c=:cyan, w=1, legend= :outertopright, label=myLabel[i])
            end
        end
        #splash/ tidal
        if pdf_C_splash_ON == 1
            if plot_pdf_choice()["E42"] == 1    #1 for pdf Normal and 2 for pdf LogNormal
                μC_mcs_splash[i] = mean(obj_C_mcs_cycle_splash.C_mcs_cycle[i])
                σC_mcs_splash[i] = std(obj_C_mcs_cycle_splash.C_mcs_cycle[i])
            elseif plot_pdf_choice()["E42"] == 2
                Vx = std(obj_C_mcs_cycle_splash.C_mcs_cycle[i])/ mean(obj_C_mcs_cycle_splash.C_mcs_cycle[i])
                μC_mcs_splash[i] = log(mean(obj_C_mcs_cycle_splash.C_mcs_cycle[i])) - 0.5*log(Vx^2 + 1)
                σC_mcs_splash[i] = √(log(Vx^2 + 1))
            end

            if i == 1
                plot_splash_25 = plot!(plot_Ccr_2, t_used/(handle*h)*2, pdf_used_C(μC_mcs_splash[i], σC_mcs_splash[i], t_used/(handle*h)*2), 
                xticks =0:2:15, color=:black, w=1, legend= :outertopright, label=myLabel[i])
            end
            if i == 2
                plot_splash_25_50 = plot!(plot_splash_25, t_used/(handle*h)*2, pdf_used_C(μC_mcs_splash[i], σC_mcs_splash[i], t_used/(handle*h)*2), 
                xticks =0:2:15, color=:blue, w=1, legend= :outertopright, label=myLabel[i])
            end
            if i == 3
                plot_splash_25_50_75 = plot!(plot_splash_25_50, t_used/(handle*h)*2, pdf_used_C(μC_mcs_splash[i], σC_mcs_splash[i], t_used/(handle*h)*2), 
                xticks =0:2:15, xlabel="C [kg/m³]", ylabel="PDF", c=:cyan, w=1, legend= :outertopright, label=myLabel[i])
            end
        end
        #submersion
        if pdf_C_sub_ON == 1
            if plot_pdf_choice()["E42"] == 1    #1 for pdf Normal and 2 for pdf LogNormal
                μC_mcs_sub[i] = mean(obj_C_mcs_cycle_sub.C_mcs_cycle[i])
                σC_mcs_sub[i] = std(obj_C_mcs_cycle_sub.C_mcs_cycle[i])
            elseif plot_pdf_choice()["E42"] == 2
                Vx = std(obj_C_mcs_cycle_sub.C_mcs_cycle[i])/ mean(obj_C_mcs_cycle_sub.C_mcs_cycle[i])
                μC_mcs_sub[i] = log(mean(obj_C_mcs_cycle_sub.C_mcs_cycle[i])) - 0.5*log(Vx^2 + 1)
                σC_mcs_sub[i] = √(log(Vx^2 + 1))
            end

            if i == 1
                plot_sub_25 = plot!(plot_Ccr_3, t_used/(handle*h), pdf_used_C(μC_mcs_sub[i], σC_mcs_sub[i], t_used/(handle*h)), 
                xticks =0:1:10, color=:black, w=1, legend= :outertopright, label=myLabel[i])
            end
            if i == 2
                plot_sub_25_50 = plot!(plot_sub_25, t_used/(handle*h), pdf_used_C(μC_mcs_sub[i], σC_mcs_sub[i], t_used/(handle*h)), 
                xticks =0:1:10, color=:blue, w=1, legend= :outertopright, label=myLabel[i])
            end
            if i == 3
                plot_sub_25_50_75 = plot!(plot_sub_25_50, t_used/(handle*h), pdf_used_C(μC_mcs_sub[i], σC_mcs_sub[i], t_used/(handle*h)), 
                xticks =0:1:10, xlabel="C [kg/m³]", ylabel="PDF", c=:cyan, w=1, legend= :outertopright, label=myLabel[i])
            end
        end
    end

    println("this information displays only for the plotting choice(s):")
    println("μ C(x,t) atmospheric, respectively(25, 50y, 75y): ", round.(μC_mcs_atm, digits=2))
    println("σ C(x,t) atmospheric, respectively(25, 50y, 75y): ", round.(σC_mcs_atm, digits=2))
    println("μ C(x,t) splash/tidal, respectively(25, 50y, 75y): ", round.(μC_mcs_splash, digits=2))
    println("σ C(x,t) splash/tidal, respectively(25, 50y, 75y): ", round.(σC_mcs_splash, digits=2))
    println("μ C(x,t) submersion, respectively(25, 50y, 75y): ", round.(μC_mcs_sub, digits=2))
    println("σ C(x,t) submersion, respectively(25, 50y, 75y): ", round.(σC_mcs_sub, digits=2))

    if pdf_C_atm_ON == 1 && pdf_C_splash_ON == 0 && pdf_C_sub_ON == 0
        plot(plot_atm_25_50_75, title = "Atmosperic Zone")
    elseif pdf_C_atm_ON == 0 && pdf_C_splash_ON == 1 && pdf_C_sub_ON == 0
        plot(plot_splash_25_50_75, title = "Splash/Tidal Zone")
    elseif pdf_C_atm_ON == 0 && pdf_C_splash_ON == 0 && pdf_C_sub_ON == 1
        plot(plot_sub_25_50_75, title = "Submersion Zone")
    elseif pdf_C_atm_ON == 1 && pdf_C_splash_ON == 1 && pdf_C_sub_ON == 0
        plot(plot(plot_atm_25_50_75, title = "Atmosperic Zone"),
        plot(plot_splash_25_50_75, title = "Splash/Tidal Zone"))
    elseif pdf_C_atm_ON == 1 && pdf_C_splash_ON == 0 && pdf_C_sub_ON == 1
        plot(plot(plot_atm_25_50_75, title = "Atmosperic Zone"),
        plot(plot_sub_25_50_75, title = "Submersion Zone"))
    elseif pdf_C_atm_ON == 0 && pdf_C_splash_ON == 1 && pdf_C_sub_ON == 1
        plot(plot(plot_splash_25_50_75, title = "Splash/Tidal Zone"),
        plot(plot_sub_25_50_75, title = "Submersion Zone"))
    elseif pdf_C_atm_ON == 1 && pdf_C_splash_ON == 1 && pdf_C_sub_ON == 1
        plot(plot(plot_atm_25_50_75, title = "Atmosperic Zone"),
        plot(plot_splash_25_50_75, title = "Splash/Tidal Zone"),
        plot(plot_sub_25_50_75, title = "Submersion Zone"))
    else
        println("Please choose the choices for plotting-ON pdf!")
    end
end