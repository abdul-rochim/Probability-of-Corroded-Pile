using Parameters
import XLSX
#reference : Andrzej S. Nowak & Kevin R. Collins, Reliability of Structures
function input_data()
    xf = XLSX.readxlsx("input/data.xlsx")
    sh = xf["rv_mcs"]
    return sh
end

function μdc_f()
    sh = input_data()
    var = 0
    if sh["E45"] == 1
        var = sh["E25"]/(1/(60*60*24*365))
    elseif sh["E45"] == 2
        var = log(sh["E25"]/(1/(60*60*24*365))) - 0.5 * log(sh["G25"]^2 + 1) #conversion from μx to μln(x) (Eq. 2.49 Nowak)
    else
        println("Enter 1 for Normal Distribution and 2 for LogNormal Distribution (μDc)")
    end
    return var
end

function μCcr_f()
    sh = input_data()
    var = 0
    if sh["E46"] == 1
        var = sh["E29"]
    elseif sh["E46"] == 2
        var = log(sh["E29"]) - 0.5 * log(sh["G29"]^2 + 1) #conversion from μx to μln(x) (Eq. 2.49 Nowak)
    else
        println("Enter 1 for Normal Distribution and 2 for LogNormal Distribution (μCcr)")
    end
    return var
end

@with_kw mutable struct Params_Data
    sh = input_data()
    #period cycle
    num_month = Int64(sh["E4"])        #months (number of months)
    num_year = Int64(sh["E5"])         #years (number of years)
    scale = Int64(sh["E49"])           #for t_month and t_year

    handle = Int64(sh["E50"])          #1.0, 10.0
    h = Int64(sh["E51"])               #1.0. 10.0, 100.0

    t_month = collect(range(0, num_month, length = num_month*scale))
    
    t_year = collect(range(0, num_year, length = num_year*scale))
    
    t_used = collect(range(1.0/h, num_year*handle, step = 1.0/h))

    #Data given
    wc = sh["E8"]               #water cement ratio
    m = sh["E9"]
    tref = sh["E10"]/365.0      #year (28 days in a year)
    Uc = sh["E11"]              #kJ/mol
    R = sh["E12"]/1000          #kJ/(mol.K)
    Tref = sh["E13"]            #K
    RHc = sh["E14"]             #in procentage, ex 75%
    
    #Data -> Monte Carlo Simulation
    #number of iteration
    N = Int64(sh["E17"])

    #choose for turn-on or turn-off -> the analysis of Monte Carlo
    mcs_ON = Int64(sh["E18"])   # 1 for run monte carlo analysis and 0 or others(except 1) for do not run monte carlo analysis
                                # if monte carlo simulation(MCS) is running will take "time", input 0 for not running MCS.
    
    #real time plotting or not?
    is_real_time_plotting = Int64(sh["E19"])   #1 for real time plotting and 0 for not real time plotting

    #choose erf function
    erf_which = Int64(sh["E20"])   # 1 for erf_1, 2 for erf_2, 3 for erf_3, 4 for erf_4  which (1,2,3,4) is based on aproximation numerics
                                   # and 0 or others(except 1,2,3,4) for using "erf function" from SpecialFunctions Package

    #choose time of cycle
    num_data = 3                #3 for 25, 50 and 75 years

    #Mean and CoV
    #Mean
    μx = sh["E24"]                              #meter
    μdc = μdc_f()                               #m²/y
    μCs_atm_75 = sh["E26"]                      #kg/m³      #75 years
    μCs_splash_75 = sh["E27"]                   #kg/m³      #75 years
    μCs_sub_75 = sh["E28"]                      #kg/m³      #75 years
    μCcr = μCcr_f()                             #kg/m³

    #Coefficient of variance
    cov_x = sh["G24"]
    cov_dc = sh["G25"]
    cov_Cs_atm_75 = sh["G26"]                   #75 years
    cov_Cs_splash_75 = sh["G27"]                #75 years
    cov_Cs_sub_75 = sh["G28"]                   #75 year
    cov_Ccr = sh["G29"]

    #50 years
    #Mean
    μCs_atm_50 = sh["E31"]
    μCs_splash_50 = sh["E32"]
    μCs_sub_50 = sh["E33"]
    #CoV
    cov_Cs_atm_50 = sh["G31"]
    cov_Cs_splash_50 = sh["G32"]
    cov_Cs_sub_50 = sh["G33"]

    #25 years
    #Mean
    μCs_atm_25 = sh["E35"]
    μCs_splash_25 = sh["E36"]
    μCs_sub_25 = sh["E37"]
    #CoV
    cov_Cs_atm_25 = sh["G35"]
    cov_Cs_splash_25 = sh["G36"]
    cov_Cs_sub_25 = sh["G37"]

    #choose for pdf type
    pdf_which_Ccr = Int64(sh["E41"])               #1 for pdf Normal and 2 for pdf LogNormal => Ccr
    pdf_which_C = Int64(sh["E42"])                 #1 for pdf Normal and 2 for pdf LogNormal => C

    #choose for the distribution type
    distribution_which_Dc = Int64(sh["E45"])       #1 for Normal distribution and 2 LogNormal distribution => Dc
    distribution_which_Ccr = Int64(sh["E46"])      #1 for Normal distribution and 2 LogNormal distribution => Ccr

    #plotting choices
    plotting_choice = Int64(sh["L40"])
end

#Params_Data()
#=
@unpack t_used = Params_Data()
for i in 1 : length(t_used)
    println(i, " - ", t_used[i]/12)
end=#