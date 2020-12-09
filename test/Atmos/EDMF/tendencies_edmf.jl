#####
##### Tendency types
#####

##### First order fluxes
struct SGSFlux{PV <: Union{Momentum, Energy, TotalMoisture}} <:
       TendencyDef{Flux{SecondOrder}, PV} end

##### Second order fluxes

##### Sources

const EnvVars =
    Union{en_ρatke, en_ρaθ_liq_cv, en_ρaq_tot_cv, en_ρaθ_liq_q_tot_cv}

const CVEnvVars = Union{en_ρaθ_liq_cv, en_ρaq_tot_cv, en_ρaθ_liq_q_tot_cv}

const UpVars = Union{UP_ρa, UP_ρaw, UP_ρaθ_liq, UP_ρaq_tot}

const EntrDetrVars = Union{UpVars, EnvVars}

struct EntrSource{PV <: EntrDetrVars} <: TendencyDef{Source, PV} end
struct DetrSource{PV <: EntrDetrVars} <: TendencyDef{Source, PV} end
struct BuoySource{PV <: Union{en_ρatke, UP_ρaw}} <: TendencyDef{Source, PV} end
struct BuoyPressSource{PV <: UP_ρaw} <: TendencyDef{Source, PV} end
struct TurbEntrSource{PV <: EnvVars} <: TendencyDef{Source, PV} end
struct PressSource{PV <: en_ρatke} <: TendencyDef{Source, PV} end
struct ShearSource{PV <: en_ρatke} <: TendencyDef{Source, PV} end
struct DissSource{PV <: EnvVars} <: TendencyDef{Source, PV} end
struct GradProdSource{PV <: CVEnvVars} <: TendencyDef{Source, PV} end

#####
##### Tendency definitions
#####


##### First order fluxes
##### Second order fluxes

##### Sources

function precomute(::EDMF, args, ::Source)
    @unpack atmos, state, aux = args

    N_up = n_updrafts(atmos.turbconv)
    EΔ_up = ntuple(N_up) do i
        entr_detr(atmos, atmos.turbconv.entr_detr, args, i)
    end
    up = state.turbconv.updraft
    ρ_inv = 1 / state.ρ
    a_min = turbconv.subdomains.a_min
    a_max = turbconv.subdomains.a_max
    ρa_up = vuntuple(N_up) do i
        state.ρ * enforce_unit_bounds(up[i].ρa * ρ_inv, a_min, a_max)
    end

    E_dyn, Δ_dyn, E_trb = ntuple(i -> map(x -> x[i], EΔ_up), 3)
    return (
        ts_all = recover_thermo_state_all(atmos, state, aux),
        env = environment_vars(state, aux, N_up),
        E_dyn = E_dyn,
        Δ_dyn = Δ_dyn,
        E_trb = E_trb,
        ρa_up = ρa_up,
    )
end

source(s::EntrSource{up_ρa{i}}, args) where {i} = args.precomputed.turbconv.E_dyn[i]
source(s::DetrSource{up_ρa{i}}, args) where {i} = -args.precomputed.turbconv.Δ_dyn[i]

function source(s::EntrSource{up_ρaw{i}}, args) where {i}
    @unpack env, E_dyn, E_trb = args.precomputed
    return (E_dyn[i] + E_trb[i]) * env.w
end
function source(s::DetrSource{up_ρaw{i}}, args) where {i}
    @unpack state = args
    @unpack Δ_dyn, E_trb, ρa_up = args.precomputed.turbconv
    up = state.turbconv.updraft
    return - (Δ_dyn[i] + E_trb[i]) * (up[i].ρaw / ρa_up[i])
end

function source(s::EntrSource{up_ρaθ_liq{i}}, args) where {i}
    @unpack E_trb, E_dyn, ts_all = args.precomputed.turbconv
    θ_liq_en = liquid_ice_pottemp(ts_all.en)
    return  (E_dyn[i] + E_trb[i]) * θ_liq_en
end
function source(s::DetrSource{up_ρaθ_liq{i}}, args) where {i}
    @unpack state = args
    @unpack E_trb, Δ_dyn, ρa_up = args.precomputed.turbconv
    up = state.turbconv.updraft
    return -(Δ_dyn[i] + E_trb[i]) * up[i].ρaθ_liq * 1/ρa_up[i]
end

function source(s::EntrSource{up_ρaq_tot{i}}, args) where {i}
    @unpack E_trb, E_dyn, ts_all = args.precomputed.turbconv
    q_tot_en = total_specific_humidity(ts_all.en)
    return  (E_dyn[i] + E_trb[i]) * q_tot_en
end
function source(s::DetrSource{up_ρaq_tot{i}}, args) where {i}
    @unpack state = args
    @unpack E_trb, Δ_dyn = args.precomputed.turbconv
    up = state.turbconv.updraft
    return -(Δ_dyn[i] + E_trb[i]) * up[i].ρaq_tot * 1/ρa_up[i]
end

function source(s::EntrSource{en_ρatke}, args)
    @unpack state = args
    @unpack E_trb, E_dyn, env = args.precomputed.turbconv
    gm = state
    ρ_inv = 1/gm.ρ
    en = state.turbconv.environment
    tke_en = enforce_positivity(en.ρatke) * ρ_inv / env.a
    return - (E_dyn[i] + E_trb[i]) * tke_en
end
function source(s::DetrSource{en_ρatke}, args)
    @unpack state = args
    @unpack Δ_dyn, env, ρa_up = args.precomputed.turbconv
    up = state.turbconv.updraft
    w_up_i = up[i].ρaw/ρa_up[i]
    return FT(0.5)*Δ_dyn[i] * (w_up_i - env.w) * (w_up_i - env.w)
end

source(s::EntrSource{en_ρaθ_liq_cv}, args)
source(s::DetrSource{en_ρaθ_liq_cv}, args)
source(s::EntrSource{en_ρaq_tot_cv}, args)
source(s::DetrSource{en_ρaq_tot_cv}, args)
source(s::EntrSource{en_ρaθ_liq_q_tot_cv}, args)
source(s::DetrSource{en_ρaθ_liq_q_tot_cv}, args)



source(s::TurbEntrSource{en_ρatke}, args)
source(s::TurbEntrSource{en_ρaθ_liq_cv}, args)
source(s::TurbEntrSource{en_ρaq_tot_cv}, args)
source(s::TurbEntrSource{en_ρaθ_liq_q_tot_cv}, args)

source(s::DissSource{en_ρatke}, args)
source(s::DissSource{en_ρaθ_liq_cv}, args)
source(s::DissSource{en_ρaq_tot_cv}, args)
source(s::DissSource{en_ρaθ_liq_q_tot_cv}, args)

source(s::GradProdSource{en_ρaθ_liq_cv}, args)
source(s::GradProdSource{en_ρaq_tot_cv}, args)
source(s::GradProdSource{en_ρaθ_liq_q_tot_cv}, args)

function source(s::BuoySource{en_ρatke}, args)
    @unpack state, atmos, diffusive, aux, t = args
    @unpack env, ts_all, Δ_dyn, E_trb = args.precomputed.turbconv

    l_mix, ∂b∂z_env, _ = mixing_length(
        atmos,
        atmos.turbconv.mix_len,
        state,
        diffusive,
        aux,
        t,
        Δ_dyn,
        E_trb,
        ts_all,
        env,
    )

    gm = state
    ρa₀ = gm.ρ * env.a
    K_m = atmos.turbconv.mix_len.c_m * l_mix * sqrt(tke_en)
    return -ρa₀ * K_m * ∂b∂z_env
end

function source(s::BuoySource{up_ρaw{i}}, args) where {i}
    # TODO: store buoyancy in precomputed instead of aux
    @unpack state, aux = args
    up = state.turbconv.updraft
    up_aux = aux.turbconv.updraft
    return up[i].ρa * up_aux[i].buoyancy
end

function source(s::BuoyPressSource{up_ρaw{i}}, args) where {i}
    @unpack state, atmos, diffusive, aux, t, env = args
    up = state.turbconv.updraft
    dpdz = perturbation_pressure(
        atmos,
        atmos.turbconv.pressure,
        state,
        diffusive,
        aux,
        t,
        env,
        i,
    )
    return -up[i].ρa * dpdz
end

source(s::PressSource{en_ρatke}, args)
source(s::ShearSource{en_ρatke}, args)
