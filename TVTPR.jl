using Clapeyron
import Clapeyron: TranslationModel, CubicModel
import Clapeyron: critical_data, default_locations, default_references, transform_params
import Clapeyron: @comps, translation, translation!, recombine_translation!, α_function, R̄
abstract type TVTPRTranslationModel <: TranslationModel end

struct TVTPRTranslationParam <: EoSParam
    Vc::SingleParam{Float64}
    v_shift::SingleParam{Float64}
end

@newmodelsimple TVTPRTranslation TVTPRTranslationModel TVTPRTranslationParam

TVTPRTranslation

export TVTPRTranslation

default_locations(::Type{TVTPRTranslation}) = critical_data()
default_references(::Type{TVTPRTranslation}) = ["10.1016/0378-3812(82)80002-2"]

function transform_params(::Type{TVTPRTranslation},params,components)
    # println("transforming params")
    v_shift = SingleParam("Volume shift",components,zeros(length(components)))
    v_shift.ismissingvalues .= true
    params["v_shift"] = v_shift
    return params
end

function translation(model::CubicModel,V,T,z,translation_model::TVTPRTranslation)
    res = zeros(eltype(V+T+first(z)),length(z))
    translation!(model,V,T,z,translation_model,res)
    return res
end

function translation!(model::CubicModel,V,T,z,translation_model::TVTPRTranslation,c)
    Tc = model.params.Tc.values
    Pc = model.params.Pc.values
    Vc = translation_model.params.Vc.values
    α = α_function(model,V,T,z,model.alpha)
    
    for i ∈ @comps
        Tci = Tc[i]
        Pci = Pc[i]
        RT = Tci*R̄
        Zc = Pci*Vc[i]/RT
        Tr = T/Tci

        η = -74.458*Zc+26.966
        γ = 246.78*Zc^2-107.21*Zc+12.67
        β = 0.35/(0.35+(η*abs(Tr-α[i]))^γ)

        cc = (0.3074-Zc)*RT/Pci
        c[i] = cc*β
    end
    return c
end

function Clapeyron.recombine_translation!(model::Clapeyron.CubicModel,translation_model::TVTPRTranslation)
    println("recombining translation")
    c = translation_model.params.v_shift
    translation!(model,0.0,0.0,0.0,translation_model,c.values)
    c.ismissingvalues .= false
    return translation_model
end

function TVTPR(components;
    idealmodel = BasicIdeal,
    alpha = TwuAlpha, #here just for compatibility with the notebooks.
    translation = TVTPRTranslation,
    userlocations = String[], 
    group_userlocations = String[],
    ideal_userlocations = String[],
    alpha_userlocations = String[],
    mixing_userlocations = String[],
    activity_userlocations = String[],
    translation_userlocations = String[],
    reference_state = nothing,
    verbose = false)

    activity = VTPRUNIFAC(components,
            userlocations = activity_userlocations,
            group_userlocations = group_userlocations,
            verbose = verbose)

    _components = activity.groups.components #extract pure component list

    mixing = VTPRRule

    return PR(_components;
    idealmodel = idealmodel,
    alpha = alpha,
    mixing = mixing,
    activity = activity,
    translation = translation,
    userlocations = userlocations,
    ideal_userlocations = ideal_userlocations,
    alpha_userlocations = alpha_userlocations,
    mixing_userlocations = mixing_userlocations,
    translation_userlocations = translation_userlocations,
    reference_state = reference_state,
    verbose = verbose)
end
export TVTPR