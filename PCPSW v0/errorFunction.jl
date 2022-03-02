include("rv.jl")
using SpecialFunctions
using Parameters

#Abramowitz and Stegun give several approximations
#maximum error: 5 x 10^-4
function erf_1(x::Float64)
    a1 = 0.278393
    a2 = 0.230389
    a3 = 0.000972
    a4 = 0.078108
    val =  (1 - (1 /(1 + a1*x + a2*x^2 + a3*x^3 + a4*x^4)^4))
    if(x >= 0)
        return val
    else
        return -val
    end
end
#println(erf_1(0.1))

#maximum error: 2.5 x 10^-5
function erf_2(x::Float64)
    p = 0.47047
    a1 = 0.3480242
    a2 = -0.0958798
    a3 = 0.7478556
    t = 1/(1 + p*x)
    val = 1 - (a1*t + a2*t^2 + a3*t^3) * exp(-x^2)
    if(x >= 0)
        return val
    else
        return -val
    end
end
#println(erf_2(0.1))

#maximum error: 3 x 10^-7
function erf_3(x::Float64)
    a1 = 0.0705230784
    a2 = 0.0422820123
    a3 = 0.0092705272
    a4 = 0.0001520143
    a5 = 0.0002765672
    a6 = 0.0000430638
    val =  (1 - (1 /(1 + a1*x + a2*x^2 + a3*x^3 + a4*x^4 + a5*x^5 + a6*x^6)^16))
    if(x >= 0)
        return val
    else
        return -val
    end
end
#println(erf_3(0.1))

#maximum error: 1.5 x 10^-7
function erf_4(x::Float64)
    p = 0.3275911
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    t = 1/(1 + p*x)
    val = 1 - (a1*t + a2*t^2 + a3*t^3 + a4*t^4 + a5*t^5) * exp(-x^2)
    if(x >= 0)
        return val
    else
        return -val
    end
end
#println(erf_4(0.1))

#println(erf(0.1))

#typeAliases
#erf_used = erf
@unpack erf_which = Params_Data()
if erf_which == 1
    erf_used = erf_1
elseif  erf_which == 2
    erf_used = erf_2
elseif  erf_which == 3
    erf_used = erf_3
elseif  erf_which == 4
    erf_used = erf_4
else
    erf_used = erf
end