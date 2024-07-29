using Clapeyron, CoolProp, NPZ, PyCall

CoolProp_=pyimport("CoolProp")
include("iPCSAFT.jl")
include("TVTPR.jl")

N = 6000

function phi_PR(model, t, p)
    return fugacity_coefficient(model, p, t)
end

Phi_PR = zeros(N, N)
CES=ARGS[2]
for Name in [ARGS[1]]
    if CES == "PCSAFT"
        model1 = PCSAFT([Name])
    elseif CES == "ADPCSAFT"
        model1 = ADPCSAFT([Name])
    elseif CES == "BACKSAFT"
        model1 = BACKSAFT([Name])
    elseif CES == "CKSAFT"
        model1 = CKSAFT([Name])
    elseif CES == "CPA"
        model1 = CPA([Name])
    elseif CES == "CPPCSAFT"
        model1 = CPPCSAFT([Name])
    elseif CES == "DAPT"
        model1 = DAPT([Name])
    elseif CES == "GEPCSAFT"
        model1 = GEPCSAFT([Name])
    elseif CES == "HeterogcPCPSAFT"
        model1 = HeterogcPCPSAFT([Name])
    elseif CES == "HomogcPCPSAFT"
        model1 = HomogcPCPSAFT([Name])
    elseif CES == "LJSAFT"
        model1 = LJSAFT([Name])
    elseif CES == "PCPSAFT"
        model1 = PCPSAFT([Name])
    elseif CES == "QPCPSAFT"
        model1 = QPCPSAFT([Name])
    elseif CES == "SAFTVRMie"
        model1 = SAFTVRMie([Name])
    elseif CES == "SAFTVRMie15"
        model1 = SAFTVRMie15([Name])
    elseif CES == "SAFTVRQMie"
        model1 = SAFTVRQMie([Name])
    elseif CES == "SAFTVRSMie"
        model1 = SAFTVRSMie([Name])
    elseif CES == "SAFTVRSW"
        model1 = SAFTVRSW([Name])
    elseif CES == "SAFTgammaMie"
        model1 = SAFTgammaMie([Name])
    elseif CES == "ogSAFT"
        model1 = ogSAFT([Name])
    elseif CES == "pharmaPCSAFT"
        model1 = pharmaPCSAFT([Name])
    elseif CES == "sCKSAFT"
        model1 = sCKSAFT([Name])
    elseif CES == "sPCSAFT"
        model1 = sPCSAFT([Name])
    elseif CES == "softSAFT2016"
        model1 = softSAFT2016([Name])
    elseif CES == "solidsoftSAFT"
        model1 = solidsoftSAFT([Name])
    elseif CES == "structSAFTgammaMie"
        model1 = structSAFTgammaMie([Name])
    elseif CES == "sCPA"
        model1 = sCPA([Name])
    elseif CES == "tcPR"
        model1 = tcPR([Name])
    elseif CES == "iPCSAFT"
        model1 = iPCSAFT([Name]; userlocations=["iPCSAFT_like.csv"])
    elseif CES == "TVTPR"
        model1 = TVTPR([Name])
    elseif CES == "Berthelot"
        model1 = Berthelot([Name])
    elseif CES == "Clausius"
        model1 = Clausius([Name])
    elseif CES == "KU"
        model1 = KU([Name])
    elseif CES == "PR"
        model1 = PR([Name])
    elseif CES == "PTV"
        model1 = PTV([Name])
    elseif CES == "PatelTeja"
        model1 = PatelTeja([Name])
    elseif CES == "RK"
        model1 = RK([Name])
    elseif CES == "RKPR"
        model1 = RKPR([Name])
    elseif CES == "vdW"
        model1 = vdW([Name])
    elseif CES == "EPPR78"
        model1 = EPPR78([Name])
    elseif CES == "PR78"
        model1 = PR78([Name])
    elseif CES == "PSRK"
        model1 = PSRK([Name])
    elseif CES == "QCPR"
        model1 = QCPR([Name])
    elseif CES == "SRK"
        model1 = SRK([Name])
    elseif CES == "UMRPR"
        model1 = UMRPR([Name])
    elseif CES == "cPR"
        model1 = cPR([Name])
    elseif CES == "tcRK"
        model1 = tcRK([Name])
    end
    handle = CoolProp_.AbstractState("HEOS", Name)

    pc = CoolProp_.AbstractState.p_critical(handle)
    Tmin = CoolProp_.AbstractState.Tmin(handle)
    Tmax = CoolProp_.AbstractState.Tmax(handle)
    pmax = CoolProp_.AbstractState.pmax(handle)
    
    T = LinRange(Tmin, Tmax, N)
    P = LinRange(0.01 * pc, pmax, N)
    for t in T
        # println("Temperature: ", t)
        for pr in P
            Phi_PR[T .== t, P .== pr] = phi_PR(model1, t, pr)
        end
    end
    
    filename = CES *" " * Name * " 1.npz"
    NPZ.npzwrite(filename, Dict("Phi_PR" => Phi_PR))
end

