using Clapeyron, CoolProp, NPZ, PyCall

CoolProp_=pyimport("CoolProp")
include("iPCSAFT.jl")
include("TVTPR.jl")

# N = 200
N = 6000 # I still think this is overkill

function phi_PR(model, t, p)
    return fugacity_coefficient(model, p, t)
end

Phi_PR = zeros(N, N)
CES=ARGS[2]
for Name in [ARGS[1]]
    model1 = eval(Meta.parse(CES*"([\""*Name*"\"])"))
    
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

