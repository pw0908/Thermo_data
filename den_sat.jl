using Clapeyron, CoolProp, NPZ, PyCall, JSON

include("TVTPR.jl")
include("iPCSAFT.jl")
CoolProp_=pyimport("CoolProp")
np=pyimport("numpy")


function vapor_pressure(temperature, EOS)
    EOS.update(CoolProp_.QT_INPUTS, 1, temperature)
    return EOS.p()
end

function error_den_vp(Name,CES)
    EXP = CoolProp_.AbstractState("HEOS", Name)
    
    temperatures=LinRange(CoolProp_.AbstractState.Tmin(EXP), 0.9999*CoolProp_.AbstractState.T_critical(EXP), 70)
    Error_liq_sat = zeros(length(temperatures))
    Error_gas_sat = zeros(length(temperatures))

    Pvap = zeros(length(temperatures))
    handle = eval(Meta.parse(CES*"([\""*Name*"\"])"))

    for i in 1:length(temperatures)
        Pvap[i] = vapor_pressure(temperatures[i], EXP)

        EXP.update(CoolProp_.QT_INPUTS , 0, temperatures[i])
        Liq_den_exp=EXP.rhomolar()
        
        EXP.update(CoolProp_.QT_INPUTS , 1, temperatures[i])
        Gas_den_exp=EXP.rhomolar()

        Psat,Vl,Vg = saturation_pressure(handle, temperatures[i])
        
        Liq_den_eos = 1/Vl
        Gas_den_eos = 1/Vg

        Error_liq_sat[i] = abs( Liq_den_exp - Liq_den_eos ) * 100 / Liq_den_exp
        Error_gas_sat[i] = abs( Gas_den_exp - Gas_den_eos ) * 100 / Gas_den_exp

    end
    
    println("Errors_L: ",Error_liq_sat)
    println("Errors_G: ",Error_gas_sat)
    npzwrite(CES*"_"*Name*"_density_sat.npz", Dict("Temperature" => temperatures, "VaporPressure" => Pvap, "Error_liq_sat" => Error_liq_sat,"Error_gas_sat"=>Error_gas_sat))
    
    return Error_liq_sat,Error_gas_sat
end

error_den_vp(ARGS[1],ARGS[2])

