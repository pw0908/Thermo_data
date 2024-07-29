using Clapeyron, CoolProp, NPZ, PyCall

include("TVTPR.jl")
include("iPCSAFT.jl")
CoolProp_=pyimport("CoolProp")
np=pyimport("numpy")
pcp=pyimport("pubchempy")


function vapor_pressure(temperature, EOS)
    EOS.update(CoolProp_.QT_INPUTS, 1, temperature)
    return EOS.p()
end

EOS=["Berthelot",
    "Clausius",
    "KU",
    "PR",
    "PTV",
    "PatelTeja",
    "RK",
    "RKPR",
    "vdW",
    "EPPR78",
    "PR78",
    "PSRK",
    "QCPR",
    "SRK",
    "UMRPR",
    "cPR","tcRK",]

function error_fi_vp(Name,CES)
    EXP = CoolProp_.AbstractState("HEOS", Name)
    
    temperatures=LinRange(CoolProp_.AbstractState.Tmin(EXP), 0.99*CoolProp_.AbstractState.T_critical(EXP), 30)
    errors = zeros(length(temperatures))
    Pvap = zeros(length(temperatures))
    
    if CES == "PCSAFT"
        handle = PCSAFT([Name])
    elseif CES == "ADPCSAFT"
        handle = ADPCSAFT([Name])
    elseif CES == "BACKSAFT"
        handle = BACKSAFT([Name])
    elseif CES == "CKSAFT"
        handle = CKSAFT([Name])
    elseif CES == "CPA"
        handle = CPA([Name])
    elseif CES == "CPPCSAFT"
        handle = CPPCSAFT([Name])
    elseif CES == "DAPT"
        handle = DAPT([Name])
    elseif CES == "GEPCSAFT"
        handle = GEPCSAFT([Name])
    elseif CES == "HeterogcPCPSAFT"
        handle = HeterogcPCPSAFT([Name])
    elseif CES == "HomogcPCPSAFT"
        handle = HomogcPCPSAFT([Name])
    elseif CES == "LJSAFT"
        handle = LJSAFT([Name])
    elseif CES == "PCPSAFT"
        handle = PCPSAFT([Name])
    elseif CES == "QPCPSAFT"
        handle = QPCPSAFT([Name])
    elseif CES == "SAFTVRMie"
        handle = SAFTVRMie([Name])
    elseif CES == "SAFTVRMie15"
        handle = SAFTVRMie15([Name])
    elseif CES == "SAFTVRQMie"
        handle = SAFTVRQMie([Name])
    elseif CES == "SAFTVRSMie"
        handle = SAFTVRSMie([Name])
    elseif CES == "SAFTVRSW"
        handle = SAFTVRSW([Name])
    elseif CES == "SAFTgammaMie"
        handle = SAFTgammaMie([Name])
    elseif CES == "ogSAFT"
        handle = ogSAFT([Name])
    elseif CES == "pharmaPCSAFT"
        handle = pharmaPCSAFT([Name])
    elseif CES == "sCKSAFT"
        handle = sCKSAFT([Name])
    elseif CES == "sPCSAFT"
        handle = sPCSAFT([Name])
    elseif CES == "softSAFT2016"
        handle = softSAFT2016([Name])
    elseif CES == "solidsoftSAFT"
        handle = solidsoftSAFT([Name])
    elseif CES == "structSAFTgammaMie"
        handle = structSAFTgammaMie([Name])
    elseif CES == "sCPA"
        handle = sCPA([Name])
    elseif CES == "tcPR"
        handle = tcPR([Name])
    elseif CES == "iPCSAFT"
        handle = iPCSAFT([Name]; userlocations=["iPCSAFT_like.csv"])
    elseif CES == "TVTPR"
        handle = TVTPR([Name])
    elseif CES == "Berthelot"
        handle = Berthelot([Name])
    elseif CES == "Clausius"
        handle = Clausius([Name])
    elseif CES == "KU"
        handle = KU([Name])
    elseif CES == "PR"
        handle = PR([Name])
    elseif CES == "PTV"
        handle = PTV([Name])
    elseif CES == "PatelTeja"
        handle = PatelTeja([Name])
    elseif CES == "RK"
        handle = RK([Name])
    elseif CES == "RKPR"
        handle = RKPR([Name])
    elseif CES == "vdW"
        handle = vdW([Name])
    elseif CES == "EPPR78"
        handle = EPPR78([Name])
    elseif CES == "PR78"
        handle = PR78([Name])
    elseif CES == "PSRK"
        handle = PSRK([Name])
    elseif CES == "QCPR"
        handle = QCPR([Name])
    elseif CES == "SRK"
        handle = SRK([Name])
    elseif CES == "UMRPR"
        handle = UMRPR([Name])
    elseif CES == "cPR"
        handle = cPR([Name])
    elseif CES == "tcRK"
        handle = tcRK([Name])
    end

    for i in 1:length(temperatures)
        Pvap[i] = vapor_pressure(temperatures[i], EXP)
        EXP.update(CoolProp_.PQ_INPUTS, Pvap[i], 1)
        Ori_val=EXP.fugacity_coefficient(0)
        cal_val=fugacity_coefficient(handle, Pvap[i], temperatures[i])[1]
        errors[i] = abs( Ori_val-cal_val ) * 100 / Ori_val
    end


    npzwrite(CES*"_"*Name*"_VP.npz", Dict("Temperature" => temperatures, "VaporPressure" => Pvap, "Errors" => errors))
    
    return errors
end

error_fi_vp(ARGS[1],ARGS[2])






