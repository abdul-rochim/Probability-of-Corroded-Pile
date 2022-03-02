#include("rv.jl") #cannot double call include("...")

module Chloride_IC
include("errorFunction.jl")
using Statistics

export Range_Time, tMonth, tYear, tUsed
export LifeCycle, T, RH
export ChlorideDiffusionCoefficient, Dref, F1, F2, F3, Dc
export ChlorideIonConcentration, C_atm, C_splash, C_sub

#struct in Julia
struct Range_Time{T<:Array}
    tMonth::T
    tYear::T
    tUsed::T
    Range_Time{T}(tMonth::T, tYear::T, tUsed::T) where {T<:Array} = new{T}(tMonth, tYear, tUsed)
end
#Inherite
tMonth(r::Range_Time) = r.tMonth
tYear(r::Range_Time) = r.tYear
tUsed(r::Range_Time) = r.tUsed

mutable struct LifeCycle{T<:Real}
    T::Array{T}
    RH::Array{T}
    function LifeCycle{T}(j::Array{T}) where T<:Real
        this = new{T}()
        this.T = @. 282.7 + (7.4 * sin(0.000045 * j + 1.6) + 8.0 * sin(0.53 * j - 2.1))
        this.RH = @. 62.5 + (8.0 * sin(0.000045 * j + 1.6) + 6.9 * sin(0.53 * j - 1.9))
        return this;
    end
end

T(l::LifeCycle) = l.T
RH(l::LifeCycle) = l.RH

#=
T(j) = @. 282.7 + (7.4 * sin(0.000045 * j + 1.6) + 8.0 * sin(0.53 * j - 2.1))
RH(j) = @. 62.5 + (8.0 * sin(0.000045 * j + 1.6) + 6.9 * sin(0.53 * j - 1.9))

@btime T(rt.tMonth)
@btime T(rt.tMonth)

@btime T(t_month)
@btime T(t_year)

@btime LifeCycle(rt.tMonth)
@btime LifeCycle(rt.tYear)
=#

mutable struct ChlorideDiffusionCoefficient{S}
    Dref::S
    F1::Vector{S}
    F2::Vector{S}
    F3::Vector{S}
    Dc::Vector{S}

    function ChlorideDiffusionCoefficient{S}(t::Vector{S}, T::Vector{S}, RH::Vector{S}) where {S}
        @unpack wc, m, tref, Uc, R, Tref, RHc = Params_Data()
        Dref_ = 10^(-12.03 + 2.4*wc) #10^(-12.6 + 2.4*wc)
        
        F1_ = @. (tref/t)^m
        
        F2_ = @. exp((Uc/R) * (1/Tref - 1/T))

        F3_ = @. (1 + (1-RH)^4 / (1-RHc)^4)^(-1)
        
        Dc_ = @. Dref_ * F1_ * F2_ * F3_

        return new{S}(Dref_, F1_, F2_, F3_, Dc_)
    end
end

#Outer constructors
Dref(d::ChlorideDiffusionCoefficient) = d.Dref
F1(d::ChlorideDiffusionCoefficient) = d.F1
F2(d::ChlorideDiffusionCoefficient) = d.F2
F3(d::ChlorideDiffusionCoefficient) = d.F3
Dc(d::ChlorideDiffusionCoefficient) = d.Dc

#Chloride Ion Concentration
mutable struct ChlorideIonConcentration{T}
    C::Function
    C_dc_mean::Function
    function ChlorideIonConcentration{T}(d::ChlorideDiffusionCoefficient,x::T,t::Vector{T}) where T
        this = new()
        this.C = function(Cs::T)
            println("Mean of Dc :", mean(d.Dc), " m2/s")
            println("Standard Deviation of Dc :", std(d.Dc))
            println("CoV of Dc :", round(std(d.Dc)/mean(d.Dc), digits=2))
            @unpack handle = Params_Data()
            return @. Cs * (1 - erf_used( x / (2 * √(abs(d.Dc/ (1/(60*60*24*365))) * t/handle))))
        end
        this.C_dc_mean = function(Cs::T, μdc::T)
            @unpack handle, cov_dc = Params_Data()
            if input_data()["E45"] == 2
                μdc = exp(μdc + 0.5 * log(cov_dc^2 + 1))    #conversion from μln(x) to μx
            end
            return @. Cs * (1 - erf_used( x / (2 * √(μdc * t/handle))))
        end
        return this
    end
end

#Outer constructors
C_atm(c::ChlorideIonConcentration) = (@unpack μCs_atm_75 = Params_Data(); c.C(μCs_atm_75))               #75 years
C_splash(c::ChlorideIonConcentration) = (@unpack μCs_splash_75 = Params_Data(); c.C(μCs_splash_75))      #75 years
C_sub(c::ChlorideIonConcentration) = (@unpack μCs_sub_75 = Params_Data(); c.C(μCs_sub_75))               #75 years
#using mean of Dc
C_atm_dc_mean(c::ChlorideIonConcentration) = (@unpack μCs_atm_75, μdc = Params_Data(); c.C_dc_mean(μCs_atm_75, μdc))           #75 years
C_splash_dc_mean(c::ChlorideIonConcentration) = (@unpack μCs_splash_75, μdc = Params_Data(); c.C_dc_mean(μCs_splash_75, μdc))  #75 years
C_sub_dc_mean(c::ChlorideIonConcentration) = (@unpack μCs_sub_75, μdc = Params_Data(); c.C_dc_mean(μCs_sub_75, μdc))           #75 years

end #end module