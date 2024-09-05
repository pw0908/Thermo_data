using Clapeyron, NPZ, PyCall

CoolProp_ = pyimport("CoolProp")
include("iPCSAFT.jl")
include("TVTPR.jl")

handle = CoolProp_.AbstractState("HEOS", "R40")
Tmin = CoolProp_.AbstractState.Tmin(handle)
Tc = CoolProp_.AbstractState.T_critical(handle)

# Function to calculate vapor pressure across a temperature range
function calculate_vp(compound, CES)
    N = 200
    vp_values = zeros(N)

    model = eval(Meta.parse(CES * "([" * "\"" * compound * "\"])"))

    Tmax = 0.9999 * Tc
    T = LinRange(Tmin, Tmax, N)

    for i in 1:length(T)
        vp_values[i] = saturation_pressure(model, T[i])[1]
    end

    filename_vp = "vp " * CES * " " * compound * ".npz"
    NPZ.npzwrite(filename_vp, Dict("vp" => vp_values))
end

# Function to calculate enthalpy of vaporization across a temperature range
function calculate_Hv(compound, CES)
    N = 200

    model = eval(Meta.parse(CES * "([" * "\"" * compound * "\"])"))

    Tmax = 0.9999 * Tc

    T = LinRange(Tmin, Tmax, N)
    Hv_values = enthalpy_vap.(model,T)

    filename_Hv = "Hv " * CES * " " * compound * ".npz"
    NPZ.npzwrite(filename_Hv, Dict("Hv" => Hv_values))
    return Hv_values
end

# Function to calculate density matrix
function density(compound, CES)
    N = 6000
    
    model = eval(Meta.parse(CES * "([" * "\"" * compound * "\"])"))

    pc = CoolProp_.AbstractState.p_critical(handle)
    Tmax = CoolProp_.AbstractState.Tmax(handle)

    pmin = 0.001 * pc
    pmax = CoolProp_.AbstractState.pmax(handle)

    T = LinRange(Tmin, Tmax, N)
    P = LinRange(pmin, pmax, N)

    # Preallocate the density matrix
    Phi_PR = zeros(Float64, length(T), length(P))

    for t in T
        for pr in P

            density_value = 1/volume(model, pr, t)
            if isnan(density_value)
                handle.update(CoolProp_.PT_INPUTS, pr, t)
                density_value = 1/volume(model, pr, t, vol0=1/handle.rhomolar())
            end
            Phi_PR[T .== t, P .== pr] .= density_value

        end
    end

    filename =  CES * " " * compound * " density.npz"
    NPZ.npzwrite(filename, Dict("Phi_PR" => Phi_PR))
    
    return Phi_PR
end


# Example usage:
# To calculate and save the vapor pressure
calculate_vp(ARGS[1], ARGS[2])

# To calculate and save the enthalpy of vaporization
calculate_Hv(ARGS[1], ARGS[2])

# To calculate and print the density matrix
density(ARGS[1], ARGS[2])
