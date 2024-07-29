using Clapeyron
import Clapeyron: PCSAFTModel, PCSAFTParam
import Clapeyron: @f, data, eos, eos_impl, a_res, a_hc, a_disp, a_assoc, transform_params, saft_lorentz_berthelot

abstract type iPCSAFTModel <: PCSAFTModel end

struct iPCSAFTParam{T} <: EoSParam
    Mw::SingleParam{T}
    segment::SingleParam{T}
    c::SingleParam{T}
    sigma::PairParam{T}
    epsilon::PairParam{T}
    epsilon_assoc::AssocParam{T}
    bondvol::AssocParam{T}
end

function iPCSAFTParam(Mw,segment,c,sigma,epsilon,epsilon_assoc,bondvol)
    el(x) = eltype(x.values)
    el(x::AssocParam) = eltype(x.values.values)
    T = mapreduce(el,promote_type,(Mw,segment,sigma,epsilon,epsilon_assoc,bondvol))
    Mw = convert(SingleParam{T},Mw)
    segment = convert(SingleParam{T},segment)
    c = convert(SingleParam{T},c)
    sigma = convert(PairParam{T},sigma)
    epsilon = convert(PairParam{T},epsilon)
    epsilon_assoc = convert(AssocParam{T},epsilon_assoc)
    bondvol = convert(AssocParam{T},bondvol)
    return iPCSAFTParam{T}(Mw,segment,c,sigma,epsilon,epsilon_assoc,bondvol) 
end

Base.eltype(p::iPCSAFTParam{T}) where T = T

@newmodel iPCSAFT iPCSAFTModel iPCSAFTParam{T}

default_references(::Type{iPCSAFT}) = ["10.1021/acs.iecr.9b04660"]
default_locations(::Type{iPCSAFT}) = []

function transform_params(::Type{iPCSAFT},params)
    sigma = params["sigma"]
    sigma.values .*= 1E-10
    return saft_lorentz_berthelot(params)
end

export iPCSAFT

function a_res(model::iPCSAFTModel, V, T, z)
    _data = @f(data)
    return @f(a_hc,_data) + @f(a_disp,_data) + @f(a_assoc,_data)
end

function eos(model::iPCSAFTModel, V, T, z = SA[1.0])
    c = model.params.c[1]
    V = V+c
    return eos_impl(model,V,T,z)
end