using Random
using Distributions
using Statistics
using Parameters

include("errorFunction.jl")
include("gaussian.jl")

###################################################################################################
#------------------------------------- Monte Carlo Simulation ------------------------------------#
###################################################################################################
#author : Abdul Rochim

@unpack is_real_time_plotting = Params_Data()
if is_real_time_plotting == 0
    using Plots
elseif is_real_time_plotting == 1
    using GLMakie
    using Makie
else
    @assert is_real_time_plotting == 0 || is_real_time_plotting == 1 "You must enter 0 or 1!"
end

#definition struct for data(mean and std) => MCS
mutable struct TimeCycle{T<:Real}
    time_cycle::Array{Vector{T}}
    time_cycle_method::Function

    function TimeCycle{T}() where T<:Real
        this = new()
        this.time_cycle = Vector{T}[] #Array{Vector{T},1}()
        this.time_cycle_method = function(v::Vector{T})
            push!(this.time_cycle, v)
            return nothing
        end
        return this
    end
end

@with_kw mutable struct Params_Var_Cov
    var_mean_25::Vector{Float64}
    var_mean_50::Vector{Float64}
    var_mean_75::Vector{Float64}
    var_cov_25::Vector{Float64}
    var_cov_50::Vector{Float64}
    var_cov_75::Vector{Float64}

    function Params_Var_Cov(params)
        @unpack μCs_atm_25, μCs_splash_25, μCs_sub_25, 
                μCs_atm_50, μCs_splash_50, μCs_sub_50,
                μCs_atm_75, μCs_splash_75, μCs_sub_75,
                cov_Cs_atm_25, cov_Cs_splash_25, cov_Cs_sub_25,
                cov_Cs_atm_50, cov_Cs_splash_50, cov_Cs_sub_50,
                cov_Cs_atm_75, cov_Cs_splash_75, cov_Cs_sub_75 = params
        this = new()        
        this.var_mean_25 = [μCs_atm_25, μCs_splash_25, μCs_sub_25]
        this.var_mean_50 = [μCs_atm_50, μCs_splash_50, μCs_sub_50]
        this.var_mean_75 = [μCs_atm_75, μCs_splash_75, μCs_sub_75]

        this.var_cov_25 = [cov_Cs_atm_25, cov_Cs_splash_25, cov_Cs_sub_25]
        this.var_cov_50 = [cov_Cs_atm_50, cov_Cs_splash_50, cov_Cs_sub_50]
        this.var_cov_75 = [cov_Cs_atm_75, cov_Cs_splash_75, cov_Cs_sub_75]
        return this
    end
end
params_var_cov = Params_Var_Cov(Params_Data())

@with_kw mutable struct Params_mu_and_sig
    σx::Float64
    σdc::Float64
    
    μCs_atm::Vector{Float64}
    μCs_splash::Vector{Float64}
    μCs_sub::Vector{Float64}
    σCs_atm::Vector{Float64}
    σCs_splash::Vector{Float64}
    σCs_sub::Vector{Float64}
    
    σCcr::Float64

    function Params_mu_and_sig(params_vc)
        @unpack var_mean_25, var_mean_50, var_mean_75,
                var_cov_25, var_cov_50, var_cov_75 = params_vc
        @unpack cov_x, μx, cov_dc, μdc, cov_Ccr, μCcr,
                num_data = Params_Data()

        #object of TimeCycle
        var_time_cycle_mean = TimeCycle{Float64}()
        var_time_cycle_std = TimeCycle{Float64}()

        var_time_cycle_mean.time_cycle_method(var_mean_25)
        var_time_cycle_mean.time_cycle_method(var_mean_50)
        var_time_cycle_mean.time_cycle_method(var_mean_75)

        var_time_cycle_std.time_cycle_method(var_cov_25)
        var_time_cycle_std.time_cycle_method(var_cov_50)
        var_time_cycle_std.time_cycle_method(var_cov_75)

        sig_x = cov_x * μx

        if input_data()["E45"] == 1
            sig_dc = cov_dc * μdc
        elseif input_data()["E45"] == 2
            sig_dc = √(log(cov_dc^2 + 1))    #conversion from σx to σln(x) (Eq. 2.48 Nowak)
        end

        mu_Cs_atm = Vector{Float64}()
        mu_Cs_splash = Vector{Float64}()
        mu_Cs_sub = Vector{Float64}()

        for i in var_time_cycle_mean.time_cycle
            for j in 1:num_data
                if j == 1
                    push!(mu_Cs_atm,i[j])
                elseif j == 2
                    push!(mu_Cs_splash,i[j])
                elseif j == 3
                    push!(mu_Cs_sub,i[j])
                end
            end
        end

        cov_Cs_atm = Vector{Float64}()
        cov_Cs_splash = Vector{Float64}()
        cov_Cs_sub = Vector{Float64}()

        for i in var_time_cycle_std.time_cycle
            for j in 1:num_data
                if j == 1
                    push!(cov_Cs_atm,i[j])
                elseif j == 2
                    push!(cov_Cs_splash,i[j])
                elseif j == 3
                    push!(cov_Cs_sub,i[j])
                end
            end
        end

        sig_Cs_atm = Float64[]
        sig_Cs_splash = Vector{Float64}()
        sig_Cs_sub = Vector{Float64}()

        for i in 1 : num_data
            push!(sig_Cs_atm, mu_Cs_atm[i] * cov_Cs_atm[i])
            push!(sig_Cs_splash, mu_Cs_splash[i] * cov_Cs_splash[i])
            push!(sig_Cs_sub, mu_Cs_sub[i] * cov_Cs_sub[i])
        end

        if input_data()["E46"] == 1
            sig_Ccr = cov_Ccr * μCcr
        elseif input_data()["E46"] == 2
            sig_Ccr = √(log(cov_Ccr^2 + 1))  #conversion from σx to σln(x) (Eq. 2.48 Nowak)
        end

        this = new()
        this.σx = sig_x
        this.σdc = sig_dc
        this.μCs_atm = mu_Cs_atm
        this.μCs_splash = mu_Cs_splash
        this.μCs_sub = mu_Cs_sub
        this.σCs_atm = sig_Cs_atm
        this.σCs_splash = sig_Cs_splash
        this.σCs_sub = sig_Cs_sub
        this.σCcr = sig_Ccr
        return this
    end
end
params_mu_sig = Params_mu_and_sig(params_var_cov)

#=
#object of TimeCycle
var_time_cycle_mean = TimeCycle{Float64}()
var_time_cycle_std = TimeCycle{Float64}()

var_time_cycle_mean.time_cycle_method(var_mean_25)
var_time_cycle_mean.time_cycle_method(var_mean_50)
var_time_cycle_mean.time_cycle_method(var_mean_75)

var_time_cycle_std.time_cycle_method(var_cov_25)
var_time_cycle_std.time_cycle_method(var_cov_50)
var_time_cycle_std.time_cycle_method(var_cov_75)

μCs_atm = Vector{Float64}()
μCs_splash = Vector{Float64}()
μCs_sub = Vector{Float64}()

for i in var_time_cycle_mean.time_cycle
    for j in 1:num_data
        if j == 1
            push!(μCs_atm,i[j])
        elseif j == 2
            push!(μCs_splash,i[j])
        elseif j == 3
            push!(μCs_sub,i[j])
        end
    end
end

cov_Cs_atm = Vector{Float64}()
cov_Cs_splash = Vector{Float64}()
cov_Cs_sub = Vector{Float64}()

for i in var_time_cycle_std.time_cycle
    for j in 1:num_data
        if j == 1
            push!(cov_Cs_atm,i[j])
        elseif j == 2
            push!(cov_Cs_splash,i[j])
        elseif j == 3
            push!(cov_Cs_sub,i[j])
        end
    end
end

σCs_atm = Float64[]
σCs_splash = Vector{Float64}()
σCs_sub = Vector{Float64}()

#standard deviation
σx = cov_x * μx
σdc = cov_dc * μdc
for i in 1 : num_data
    push!(σCs_atm, μCs_atm[i] * cov_Cs_atm[i])
    push!(σCs_splash, μCs_splash[i] * cov_Cs_splash[i])
    push!(σCs_sub, μCs_sub[i] * cov_Cs_sub[i])
end
σCcr = cov_Ccr * μCcr
=#

#using Parameters (@with_kw)
@with_kw mutable struct Params_Gaussian
    G_x::Gaussian
    G_dc::Gaussian
    G_Cs_atm::Vector{Gaussian}
    G_Cs_splash::Vector{Gaussian}
    G_Cs_sub::Vector{Gaussian}
    G_Ccr::Gaussian

    function Params_Gaussian(params)
        @unpack σx, σdc, μCs_atm, μCs_splash, μCs_sub, σCs_atm, σCs_splash, σCs_sub, σCcr = params
        @unpack num_data, μx, μdc, μCcr = Params_Data()

        var_G_Cs_atm = Array{Gaussian,1}(undef, num_data)           #object Gaussian
        for i in 1 : num_data
            var_G_Cs_atm[i] = Gaussian(μCs_atm[i],σCs_atm[i])
        end

        var_G_Cs_splash = Array{Gaussian,1}(undef, num_data)        #object Gaussian
        for i in 1 : num_data
            var_G_Cs_splash[i] = Gaussian(μCs_splash[i],σCs_splash[i])
        end

        var_G_Cs_sub = Array{Gaussian,1}(undef, num_data)           #object Gaussian
        for i in 1 : num_data
            var_G_Cs_sub[i] = Gaussian(μCs_sub[i],σCs_sub[i])
        end

        this = new()
        this.G_x = Gaussian(μx,σx)
        this.G_dc = Gaussian(μdc,σdc)
        this.G_Cs_atm = var_G_Cs_atm
        this.G_Cs_splash = var_G_Cs_splash
        this.G_Cs_sub = var_G_Cs_sub
        this.G_Ccr = Gaussian(μCcr,σCcr)
        return this
    end
end
params_gaussian = Params_Gaussian(params_mu_sig)
#@unpack G_x, G_dc, G_Cs_atm, G_Cs_splash, G_Cs_sub, G_Ccr = params_gaussian

#=
G_Cs_atm = Array{Gaussian,1}(undef, num_data)          #object Gaussian
G_Cs_splash = Array{Gaussian,1}(undef, num_data)       #object Gaussian
G_Cs_sub = Array{Gaussian,1}(undef, num_data)          #object Gaussian

#using Gaussian function
begin
    G_x = Gaussian(μx,σx)
    G_dc = Gaussian(μdc,σdc)
    for i in 1 : num_data
        G_Cs_atm[i] = Gaussian(μCs_atm[i],σCs_atm[i])
        G_Cs_splash[i] = Gaussian(μCs_splash[i],σCs_splash[i])
        G_Cs_sub[i] = Gaussian(μCs_sub[i],σCs_sub[i])
    end
    G_Ccr = Gaussian(μCcr,σCcr)
end
=#

@with_kw mutable struct Params_Cs
    x_mcs::Vector{Float64}
    Dc_mcs::Vector{Float64}
    Cs_atm_mcs::Array{Vector{Float64}}
    Cs_splash_mcs::Array{Vector{Float64}}
    Cs_sub_mcs::Array{Vector{Float64}}

    #critical chloride value for corrosion initiation as capacity (resistance) [R(U,t)]
    Ccr_mcs::Vector{Float64}

    function Params_Cs(params)
        @unpack G_x, G_dc, G_Cs_atm, G_Cs_splash, G_Cs_sub, G_Ccr = params
        @unpack N, distribution_which_Dc, num_data, distribution_which_Ccr = Params_Data()
        this = new()
        this.x_mcs = rand(Normal(G_x.μ,G_x.σ), N)

        if distribution_which_Dc == 1
            distribution_used_Dc = Normal
        elseif distribution_which_Dc == 2
            distribution_used_Dc = LogNormal
        else
            println("Choose which distribution(Dc) type you want to use!")
        end
        this.Dc_mcs = rand(distribution_used_Dc(G_dc.μ,G_dc.σ),N)

        var_Cs_atm_mcs = Array{Vector{Float64},1}(undef, num_data)
        for i in 1 : num_data
            var_Cs_atm_mcs[i] = rand(Normal(G_Cs_atm[i].μ,G_Cs_atm[i].σ), N)
        end
        this.Cs_atm_mcs = var_Cs_atm_mcs

        var_Cs_splash_mcs = Array{Vector{Float64},1}(undef, num_data)
        for i in 1 : num_data
            var_Cs_splash_mcs[i] = rand(Normal(G_Cs_splash[i].μ,G_Cs_splash[i].σ), N)
        end
        this.Cs_splash_mcs = var_Cs_splash_mcs
        
        var_Cs_sub_mcs = Array{Vector{Float64},1}(undef, num_data)
        for i in 1 : num_data
            var_Cs_sub_mcs[i] = rand(Normal(G_Cs_sub[i].μ,G_Cs_sub[i].σ), N)
        end
        this.Cs_sub_mcs = var_Cs_sub_mcs

        #critical chloride value for corrosion initiation as capacity (resistance) [R(U,t)]
        if distribution_which_Ccr == 1
            distribution_used_Ccr = Normal
        elseif distribution_which_Ccr == 2
            distribution_used_Ccr = LogNormal
        else
            println("Choose which distribution(Ccr) type you want to use!")
        end
        this.Ccr_mcs = rand(distribution_used_Ccr(G_Ccr.μ,G_Ccr.σ), N)
        return this
    end
end
params_cs = Params_Cs(params_gaussian)

#=
function x_mcs_(params)
    @unpack G_x = params
    @unpack N = Params_Data()
    return rand(Normal(G_x.μ,G_x.σ), N)
end

function Dc_mcs_(params)
    @unpack G_dc = params
    @unpack N, distribution_which_Dc = Params_Data()
    if distribution_which_Dc == 1
        distribution_used_Dc = Normal
    elseif distribution_which_Dc == 2
        distribution_used_Dc = LogNormal
    else
        println("Choose which distribution(Dc) type you want to use!")
    end

    return rand(distribution_used_Dc(G_dc.μ,G_dc.σ),N)
end

function Cs_atm_mcs_(params)
    @unpack G_Cs_atm = params
    @unpack num_data, N = Params_Data()
    var = Array{Vector{Float64},1}()
    for i in 1 : num_data
        push!(var, rand(Normal(G_Cs_atm[i].μ,G_Cs_atm[i].σ), N))
    end
    return var
end

function Cs_splash_mcs_(params)
    @unpack G_Cs_splash = params
    @unpack num_data, N = Params_Data()
    var = Array{Vector{Float64},1}()
    for i in 1 : num_data
        push!(var, rand(Normal(G_Cs_splash[i].μ,G_Cs_splash[i].σ), N))
    end
    return var
end

function Cs_sub_mcs_(params)
    @unpack G_Cs_sub = params
    @unpack num_data, N = Params_Data()
    var = Vector{Float64}[]
    for i in 1 : num_data
        push!(var, rand(Normal(G_Cs_sub[i].μ,G_Cs_sub[i].σ), N))
    end
    return var
end

#critical chloride value for corrosion initiation as capacity (resistance) [R(U,t)]
function Ccr_mcs_(params)
    @unpack G_Ccr = params
    @unpack N, distribution_which_Ccr = Params_Data()
    if distribution_which_Ccr == 1
        distribution_used_Ccr = Normal
    elseif distribution_which_Ccr == 2
        distribution_used_Ccr = LogNormal
    else
        println("Choose which distribution(Ccr) type you want to use!")
    end

    return rand(distribution_used_Ccr(G_Ccr.μ,G_Ccr.σ), N)
end

params_cs = Params_Cs(x_mcs_(params_gaussian), Dc_mcs_(params_gaussian), Cs_atm_mcs_(params_gaussian), 
Cs_splash_mcs_(params_gaussian), Cs_sub_mcs_(params_gaussian), Ccr_mcs_(params_gaussian))
=#

#=
Cs_atm_mcs = Array{Vector{Float64},1}()
Cs_splash_mcs = Array{Vector{Float64},1}()
Cs_sub_mcs = Vector{Float64}[]

#data distribution
begin
    x_mcs = rand(Normal(G_x.μ,G_x.σ), N)
    Dc_mcs = rand(Normal(G_dc.μ,G_dc.σ),N)
    for i in 1 : num_data
        push!(Cs_atm_mcs, rand(Normal(G_Cs_atm[i].μ,G_Cs_atm[i].σ), N))
        push!(Cs_splash_mcs, rand(Normal(G_Cs_splash[i].μ,G_Cs_splash[i].σ), N))
        push!(Cs_sub_mcs, rand(Normal(G_Cs_sub[i].μ,G_Cs_sub[i].σ), N))
    end
    #critical chloride value for corrosion initiation as capacity (resistance) [R(U,t)]
    Ccr_mcs = rand(Normal(G_Ccr.μ,G_Ccr.σ), N)
end
=#

#definition C(x,t) Vector => chloride ion concentration (cycle)
mutable struct C_MCS{T<:Real}
    c_cp_mcs::Array{Vector{T}}
    C_mcs_cycle::Array{Vector{T}}
    C_mcs_cycle_method::Function

    function C_MCS{T}() where T<:Real
        this = new()
        this.c_cp_mcs = Array{Vector{T},1}()
        this.C_mcs_cycle = Array{Vector{T},1}()
        
        this.C_mcs_cycle_method = function(u::Vector{T}, v::Vector{T}, w::Vector{T}, x::Vector{T})
            push!(this.c_cp_mcs, u)         #cumulative of corrosion probability (of MCS)
            push!(this.C_mcs_cycle, v)      #C(chloride content cycle) -> 25y
            push!(this.C_mcs_cycle, w)      #C(chloride content cycle) -> 50y
            push!(this.C_mcs_cycle, x)      #C(chloride content cycle) -> 75y
            return nothing
        end
        return this
    end
end

params_object_mcs = (
    obj_C_mcs_cycle_atm = C_MCS{Float64}(),
    obj_C_mcs_cycle_splash = C_MCS{Float64}(),
    obj_C_mcs_cycle_sub = C_MCS{Float64}()
)

#analysis
function mcs(params_object)
    @unpack N, mcs_ON, num_data, handle, t_used, μCcr, cov_Ccr = Params_Data()
    obj_C_mcs_cycle_atm, obj_C_mcs_cycle_splash, obj_C_mcs_cycle_sub = params_object
    μCcr_convert = exp(μCcr + 0.5 * log(cov_Ccr^2 + 1))    #conversion from μln(x) to μx
    myLabel = Vector{String}()
    push!(myLabel, "25 years")
    push!(myLabel, "50 years")
    push!(myLabel, "75 years")

    if mcs_ON == 1
        println("\n********** Analysis starts **********")
        println("Monte Carlo Simulation is running...")
        println("It will take time...\n")
        #G(U,t) = R(U,t) - S(U,t)
        #U is set of random variables(RVs) and t represents the time
        #G > 0 represents the safe region. while G < 0 refers to the failure region
   
        println("prepare for calculation...")       
        #initial
        cp_atm = 0          #corrosion probability (atmospheric)
        cp_splash = 0       #corrosion probability (splash/ tidal)
        cp_sub = 0          #corrosion probability (submersion)
        
        #initialist
        c_cp_atm = zeros(length(t_used))            #cumulative of corrosion probability (atmospheric)
        c_cp_splash = zeros(length(t_used))         #cumulative of corrosion probability (splash/ tidal)
        c_cp_sub = zeros(length(t_used))            #cumulative of corrosion probability (submersion)

        #surface chloride content as demand or stressor [S(U,t)]
        C_mcs_atm_25 = zeros(convert(Int64,round(length(t_used)/3, digits=0)))      #chloride ion concentration (atmospheric)
        C_mcs_splash_25 = zeros(convert(Int64,round(length(t_used)/3, digits=0)))   #chloride ion concentration (splash/ tidal)
        C_mcs_sub_25 = zeros(convert(Int64,round(length(t_used)/3, digits=0)))      #chloride ion concentration (submersion)

        C_mcs_atm_25_calc = fill(zeros(N),convert(Int64,round(length(t_used)/3, digits=0)))
        C_mcs_splash_25_calc = fill(zeros(N),convert(Int64,round(length(t_used)/3, digits=0)))
        C_mcs_sub_25_calc = fill(zeros(N),convert(Int64,round(length(t_used)/3, digits=0)))

        C_mcs_atm_50 = zeros(convert(Int64,round(length(t_used)*2/3, digits=0)))
        C_mcs_splash_50 = zeros(convert(Int64,round(length(t_used)*2/3, digits=0)))
        C_mcs_sub_50 = zeros(convert(Int64,round(length(t_used)*2/3, digits=0)))

        C_mcs_atm_50_calc = fill(zeros(N),convert(Int64,round(length(t_used)*2/3, digits=0)))
        C_mcs_splash_50_calc = fill(zeros(N),convert(Int64,round(length(t_used)*2/3, digits=0)))
        C_mcs_sub_50_calc = fill(zeros(N),convert(Int64,round(length(t_used)*2/3, digits=0)))

        C_mcs_atm_75 = zeros(length(t_used))
        C_mcs_splash_75 = zeros(length(t_used))
        C_mcs_sub_75 = zeros(length(t_used))

        C_mcs_atm_75_calc = fill(zeros(N),length(t_used))
        C_mcs_splash_75_calc = fill(zeros(N),length(t_used))
        C_mcs_sub_75_calc = fill(zeros(N),length(t_used))

        @unpack x_mcs, Dc_mcs, Cs_atm_mcs, Cs_splash_mcs, Cs_sub_mcs, Ccr_mcs = params_cs
        h = 0       #initial
        ax_at_Ccr_atm = length(t_used)
        ax_at_Ccr_splash = length(t_used)
        ax_at_Ccr_sub = length(t_used)
        ax_at_10percent_atm = length(t_used)
        ax_at_10percent_splash = length(t_used)
        ax_at_10percent_sub = length(t_used)
        if is_real_time_plotting == 0
            println("the next plot is not live plotting, so wait 'till the plot shown...\n")
        elseif is_real_time_plotting == 1
            @unpack num_year = Params_Data()
            #set_theme!(theme_black())
            points_atm = Observable(Point2f[(0.0,0.0)])
            points_splash = Observable(Point2f[(0.0,0.0)])
            points_sub = Observable(Point2f[zeros(2)])
            points_10percent = Observable(Point2f[(0,0.1)])
            points_50percent = Observable(Point2f[(0,0.5)])
            fig, ax = lines(points_atm, color= :red, linewidth= :3,
                figure = (resolution=(750,750),), axis= (xlabel = "t [year]",
                ylabel = "Prob. of Corrosion", xticks = 0:25:75, yticks = 0:0.2:1.0,))
            lines!(points_splash, color= :blue, linewidth= :3)
            lines!(points_sub, color= :black, linewidth= :3)
            lines!(points_10percent, color= :magenta, linewidth= :2, linestyle=:dash)
            lines!(points_50percent, color= :magenta, linewidth= :2, linestyle=:dash)
            text!("Atmosperic", color =:red, position = Point2f(3.0, 0.95),
            textsize = 15, align = (:left, :center))
            text!("Splash/Tidal", color =:blue, position = Point2f(3.0, 0.90),
            textsize = 15, align = (:left, :center))
            text!("Submersion", color =:black, position = Point2f(3.0, 0.85),
            textsize = 15, align = (:left, :center))
            t = t_used
            limits!(ax, 0, num_year, 0, 1)
            fps = 120        #30, 60, 120, ...  frame per second
            #nframes = length(t)
            display(fig) #do not display if using record
        end #take this "end" below =>as the end of elseif statement
            #record(fig, joinpath(@__DIR__, "output", "mcs_plot.mp4")) do io
        for i in eachindex(t_used) # = 1 : length(t_used) #nframes
            if i <= convert(Int64,round(length(t_used)/3, digits=0))
                h = 1       #h = 1 for 25 years
                
                if i == 1
                    println("initial calculation for period ", 0, " to ", myLabel[h],
                    " on ", i, "st data")
                end

                for j = 1 : N
                    C_mcs_atm_25_calc[i][j] = Cs_atm_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))
                    C_mcs_splash_25_calc[i][j] = Cs_splash_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))
                    C_mcs_sub_25_calc[i][j] = Cs_sub_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))
                end

                C_mcs_atm_25[i] = mean(C_mcs_atm_25_calc[i])
                C_mcs_splash_25[i] = mean(C_mcs_splash_25_calc[i])
                C_mcs_sub_25[i] = mean(C_mcs_sub_25_calc[i])

                if i == convert(Int64,round(length(t_used)/3, digits=0))
                    println("end of data 25y at: ", i, "\n")
                end
            end

            if i <= convert(Int64,round(length(t_used)*2/3, digits=0))
                h = 2       #h = 2 for 50 years

                if i == convert(Int64,round(length(t_used)/3+1, digits=0))
                    println("initial calculation for period ", myLabel[h-1], " to ", myLabel[h],
                    " on ", i, "th data")
                end

                for j = 1 : N
                    C_mcs_atm_50_calc[i][j] = Cs_atm_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))
                    C_mcs_splash_50_calc[i][j] = Cs_splash_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))
                    C_mcs_sub_50_calc[i][j] = Cs_sub_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))
                end

                C_mcs_atm_50[i] = mean(C_mcs_atm_50_calc[i])
                C_mcs_splash_50[i] = mean(C_mcs_splash_50_calc[i])
                C_mcs_sub_50[i] = mean(C_mcs_sub_50_calc[i])

                if i == convert(Int64,round(length(t_used)*2/3, digits=0))
                    println("end of data 50y at: ", i, "\n")
                end
            end

            if i <= length(t_used)
                h = 3       #h = 3 for 75 years
                
                if i == convert(Int64,round(length(t_used)*2/3+1, digits=0))
                    println("initial calculation for period ", myLabel[h-1], " to ", myLabel[h],
                    " on ", i, "th data")
                end

                for j = 1 : N
                    C_mcs_atm_75_calc[i][j] = Cs_atm_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))   
                    C_mcs_splash_75_calc[i][j] = Cs_splash_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))   
                    C_mcs_sub_75_calc[i][j] = Cs_sub_mcs[h][j]*(1-erf_used(x_mcs[j]/(2*√(abs(Dc_mcs[j])*t_used[i]/handle))))    
                    #S(U,t) > R(U,t)  ->  failure region
                    #75 years
                    if(C_mcs_atm_75_calc[i][j] > Ccr_mcs[j])
                        cp_atm += 1
                    end
                    if(C_mcs_splash_75_calc[i][j] > Ccr_mcs[j])
                        cp_splash += 1
                    end
                    if(C_mcs_sub_75_calc[i][j] > Ccr_mcs[j])
                        cp_sub += 1
                    end
                end
                
                C_mcs_atm_75[i] = mean(C_mcs_atm_75_calc[i])
                C_mcs_splash_75[i] = mean(C_mcs_splash_75_calc[i])
                C_mcs_sub_75[i] = mean(C_mcs_sub_75_calc[i])
                
                c_cp_atm[i] = Float64(cp_atm)/(N*length(t_used))
                c_cp_splash[i] = Float64(cp_splash)/(N*length(t_used))
                c_cp_sub[i] = Float64(cp_sub)/(N*length(t_used))                
                if is_real_time_plotting == 1
                    new_points_atm = [t[i]/handle, c_cp_atm[i]]
                    new_points_splash = [t[i]/handle, c_cp_splash[i]]
                    new_points_sub = [t[i]/handle, c_cp_sub[i]]
                    new_points_10p = [length(t)/(handle*h), 0.1]
                    new_points_50p = [length(t)/(handle*h), 0.5]
                    ax.title ="t = $(round(t[i]/handle,digits=2))"
                    points_atm[] = push!(points_atm[], new_points_atm)
                    points_splash[] = push!(points_splash[], new_points_splash)
                    points_sub[] = push!(points_sub[], new_points_sub)
                    points_10percent[] = push!(points_10percent[], new_points_10p)
                    points_50percent[] = push!(points_50percent[], new_points_50p)
                    sleep(1/fps) #refreshes the display!
                    #recordframe!(io) #for recording frame                    
                end

                if input_data()["E46"] == 1
                    if C_mcs_atm_75[i] >= μCcr && i < ax_at_Ccr_atm
                        ax_at_Ccr_atm = i
                    end

                    if C_mcs_splash_75[i] >= μCcr && i < ax_at_Ccr_splash
                        ax_at_Ccr_splash = i
                    end

                    if C_mcs_sub_75[i] >= μCcr && i < ax_at_Ccr_sub
                        ax_at_Ccr_sub = i
                    end
                elseif input_data()["E46"] == 2
                    if C_mcs_atm_75[i] >= μCcr_convert && i < ax_at_Ccr_atm
                        ax_at_Ccr_atm = i
                    end

                    if C_mcs_splash_75[i] >= μCcr_convert && i < ax_at_Ccr_splash
                        ax_at_Ccr_splash = i
                    end

                    if C_mcs_sub_75[i] >= μCcr_convert && i < ax_at_Ccr_sub
                        ax_at_Ccr_sub = i
                    end
                end

                if c_cp_atm[i] >= 0.1 && i < ax_at_10percent_atm #10%
                    ax_at_10percent_atm = i + 1
                end

                if c_cp_splash[i] >= 0.1 && i < ax_at_10percent_splash #10%
                    ax_at_10percent_splash = i + 1
                end

                if c_cp_sub[i] >= 0.1 && i < ax_at_10percent_sub #10%
                    ax_at_10percent_sub = i + 1
                end

                if i == convert(Int64,round(length(t_used), digits=0))
                    println("end of data 75y at: ", i, "\n")
                end
            end
        end #end of 75 years
            #end #end record
        #end #end of elseif

        println("The chloride ion concentration achieve at critical chloride value:")
        if input_data()["E46"] == 1
            if(C_mcs_atm_75[ax_at_Ccr_atm] >= μCcr)
                println("\tin atmospheric zone was ", round(t_used[ax_at_Ccr_atm]/handle, digits=2), " years")
            end
            if(C_mcs_splash_75[ax_at_Ccr_splash] >= μCcr)
                println("\tin splash/tidal zone was ", round(t_used[ax_at_Ccr_splash]/handle, digits=2), " years")
            end
            if(C_mcs_sub_75[ax_at_Ccr_sub] >= μCcr)
                println("\tin submersion zone was ", round(t_used[ax_at_Ccr_sub]/handle, digits=2), " years")
            end
        elseif input_data()["E46"] == 2
            if(C_mcs_atm_75[ax_at_Ccr_atm] >= μCcr_convert)
                println("\tin atmospheric zone was ", round(t_used[ax_at_Ccr_atm]/handle, digits=2), " years")
            end
            if(C_mcs_splash_75[ax_at_Ccr_splash] >= μCcr_convert)
                println("\tin splash/tidal zone was ", round(t_used[ax_at_Ccr_splash]/handle, digits=2), " years")
            end
            if(C_mcs_sub_75[ax_at_Ccr_sub] >= μCcr_convert)
                println("\tin submersion zone was ", round(t_used[ax_at_Ccr_sub]/handle, digits=2), " years")
            end
        end
        
        println("\nThe corrosion initiation time obtained from deterministic approach, which is close to 10% probability of corrosion occurance criteria:")
        if(c_cp_atm[ax_at_10percent_atm] >= 0.1) #10%
            println("\tin atmospheric zone was ", round(t_used[ax_at_10percent_atm-1]/handle, digits=2), " years ", )
        end
        if(c_cp_splash[ax_at_10percent_splash] >= 0.1) #10%
            println("\tin splash/tidal zone was ", round(t_used[ax_at_10percent_splash-1]/handle, digits=2), " years ", )
        end
        if(c_cp_sub[ax_at_10percent_sub] >= 0.1) #10%
            println("\tin submersion zone was ", round(t_used[ax_at_10percent_sub-1]/handle, digits=2), " years ", )
        end

        println("\nMaximum Probability of Corrosion each zone in 75 years, respectively:")
        println("\tin atmospheric zone, ", round(maximum(c_cp_atm)*100, digits=2), "%")
        println("\tin splash/tidal zone, ", round(maximum(c_cp_splash)*100, digits=2), "%")
        println("\tin submersion zone, ", round(maximum(c_cp_sub)*100, digits=2), "%\n")

        obj_C_mcs_cycle_atm.C_mcs_cycle_method(c_cp_atm, C_mcs_atm_25, C_mcs_atm_50, C_mcs_atm_75)
        obj_C_mcs_cycle_splash.C_mcs_cycle_method(c_cp_splash, C_mcs_splash_25, C_mcs_splash_50, C_mcs_splash_75)
        obj_C_mcs_cycle_sub.C_mcs_cycle_method(c_cp_sub, C_mcs_sub_25, C_mcs_sub_50, C_mcs_sub_75)

        println("The result of analysis:")
        for k in 1:num_data
            println("μ C(x,t) atmospheric ", myLabel[k], " : ", round(mean(obj_C_mcs_cycle_atm.C_mcs_cycle[k]), digits=2), " kg/m³")
            println("σ C(x,t) atmospheric ", myLabel[k], " : ", round(std(obj_C_mcs_cycle_atm.C_mcs_cycle[k]), digits=2))
            println("μ C(x,t) splash/tidal ", myLabel[k], " : ", round(mean(obj_C_mcs_cycle_splash.C_mcs_cycle[k]), digits=2), " kg/m³")
            println("σ C(x,t) splash/tidal ", myLabel[k], " : ", round(std(obj_C_mcs_cycle_splash.C_mcs_cycle[k]), digits=2))
            println("μ C(x,t) submersion ", myLabel[k], " : ", round(mean(obj_C_mcs_cycle_sub.C_mcs_cycle[k]), digits=2), " kg/m³")
            println("σ C(x,t) submersion ", myLabel[k], " : ", round(std(obj_C_mcs_cycle_sub.C_mcs_cycle[k]), digits=2))
            #println()
        end
        println("********** End of Analysis **********\n")
    else
        println("Monte Carlo Simulation is NOT running!\n")
    end
end #end of analysis


### this code below does not a part of the analysis ###
#find the minimum Dc at i-th
#=
n = 12*h
trial = zeros(Int64(n))
for i in 1: Int64(n)
    trial[i] = Choride_IC.Dc(dc_)[i]
end
println(trial)
for i in 1 : length(trial)
    if trial[i] == minimum(trial)
        println("i at minimum : ", i)
    end
end
plot(1:Int64(n), trial)
=#
