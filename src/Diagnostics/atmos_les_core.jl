using ..Atmos
using ..Atmos: recover_thermo_state
using ..Mesh.Topologies
using ..Mesh.Grids
using ..Thermodynamics
using LinearAlgebra

"""
    setup_atmos_core_diagnostics(
        ::AtmosLESConfigType,
        interval::String,
        out_prefix::String;
        writer::AbstractWriter,
        interpol = nothing,
    )

Create the "AtmosLESCore" `DiagnosticsGroup` which contains the following
diagnostic variables, all of which are density-averaged horizontal averages
conditional upon `q_liq > 0 && w > 0` except for `core_frac`:

- core_frac: cloud core fraction
- u_core: cloud core x-velocity
- v_core: cloud core y-velocity
- w_core: cloud core z-velocity
- avg_rho_core: cloud core air density (_not_ density averaged)
- rho_core: cloud core air density
- qt_core: cloud core total specific humidity
- ql_core: cloud core liquid water specific humidity
- thv_core: cloud core virtual potential temperature
- thl_core: cloud core liquid-ice potential temperature
- ei_core: cloud core specific internal energy
- var_u_core: cloud core variance of x-velocity
- var_v_core: cloud core variance of y-velocity
- var_w_core: cloud core variance of z-velocity
- var_qt_core: cloud core variance of total specific humidity
- var_thl_core: cloud core variance of liquid-ice potential temperature
- var_ei_core: cloud core variance of specific internal energy
- cov_w_rho_core: cloud core vertical eddy flux of density
- cov_w_qt_core: cloud core vertical eddy flux of specific humidity
- cov_w_thl_core: cloud core vertical eddy flux of liquid-ice potential temperature
- cov_w_ei_core: cloud core vertical eddy flux of specific internal energy
- cov_qt_thl_core: cloud core covariance of total specific humidity and liquid-ice potential temperature
- cov_qt_ei_core: cloud core covariance of total specific humidity and specific internal energy

All these variables are output with the `z` dimension (`x3id`) on the DG grid
(`interpol` may _not_ be specified) as well as a (unlimited) `time` dimension
at the specified `interval`.
"""
function setup_atmos_core_diagnostics(
    ::AtmosLESConfigType,
    interval::String,
    out_prefix::String;
    writer = NetCDFWriter(),
    interpol = nothing,
)
    # TODO: remove this
    @assert isnothing(interpol)

    return DiagnosticsGroup(
        "AtmosLESCore",
        Diagnostics.atmos_les_core_init,
        Diagnostics.atmos_les_core_fini,
        Diagnostics.atmos_les_core_collect,
        interval,
        out_prefix,
        writer,
        interpol,
    )
end

# Simple horizontal averages
function vars_atmos_les_core_simple(m::AtmosModel, FT)
    @vars begin
        u_core::FT
        v_core::FT
        w_core::FT
        avg_rho_core::FT        # ρ
        rho_core::FT            # ρρ
        qt_core::FT             # q_tot
        ql_core::FT             # q_liq
        thv_core::FT            # θ_vir
        thl_core::FT            # θ_liq
        ei_core::FT             # e_int
    end
end
num_atmos_les_core_simple_vars(m, FT) =
    varsize(vars_atmos_les_core_simple(m, FT))
atmos_les_core_simple_vars(m, array) =
    Vars{vars_atmos_les_core_simple(m, eltype(array))}(array)

function atmos_les_core_simple_sums!(atmos::AtmosModel, state, thermo, MH, sums)
    sums.u_core += MH * state.ρu[1]
    sums.v_core += MH * state.ρu[2]
    sums.w_core += MH * state.ρu[3]
    sums.avg_rho_core += MH * state.ρ
    sums.rho_core += MH * state.ρ * state.ρ
    sums.qt_core += MH * state.moisture.ρq_tot
    sums.ql_core += MH * thermo.moisture.q_liq * state.ρ
    sums.thv_core += MH * thermo.moisture.θ_vir * state.ρ
    sums.thl_core += MH * thermo.moisture.θ_liq_ice * state.ρ
    sums.ei_core += MH * thermo.e_int * state.ρ

    return nothing
end

# Variances and covariances
function vars_atmos_les_core_ho(m::AtmosModel, FT)
    @vars begin
        var_u_core::FT          # u′u′
        var_v_core::FT          # v′v′
        var_w_core::FT          # w′w′
        var_qt_core::FT         # q_tot′q_tot′
        var_thl_core::FT        # θ_liq_ice′θ_liq_ice′
        var_ei_core::FT         # e_int′e_int′

        cov_w_rho_core::FT      # w′ρ′
        cov_w_qt_core::FT       # w′q_tot′
        cov_w_thl_core::FT      # w′θ_liq_ice′
        cov_w_ei_core::FT       # w′e_int′
        cov_qt_thl_core::FT     # q_tot′θ_liq_ice′
        cov_qt_ei_core::FT      # q_tot′e_int′
    end
end
num_atmos_les_core_ho_vars(m, FT) = varsize(vars_atmos_les_core_ho(m, FT))
atmos_les_core_ho_vars(m, array) =
    Vars{vars_atmos_les_core_ho(m, eltype(array))}(array)

function atmos_les_core_ho_sums!(atmos::AtmosModel, state, thermo, MH, ha, sums)
    u = state.ρu[1] / state.ρ
    u′ = u - ha.u_core
    v = state.ρu[2] / state.ρ
    v′ = v - ha.v_core
    w = state.ρu[3] / state.ρ
    w′ = w - ha.w_core
    q_tot = state.moisture.ρq_tot / state.ρ
    q_tot′ = q_tot - ha.qt_core
    θ_liq_ice′ = thermo.moisture.θ_liq_ice - ha.thl_core
    e_int′ = thermo.e_int - ha.ei_core

    sums.var_u_core += MH * u′^2 * state.ρ
    sums.var_v_core += MH * v′^2 * state.ρ
    sums.var_w_core += MH * w′^2 * state.ρ
    sums.var_qt_core += MH * q_tot′^2 * state.ρ
    sums.var_thl_core += MH * θ_liq_ice′^2 * state.ρ
    sums.var_ei_core += MH * e_int′^2 * state.ρ

    sums.cov_w_rho_core += MH * w′ * (state.ρ - ha.avg_rho_core) * state.ρ
    sums.cov_w_qt_core += MH * w′ * q_tot′ * state.ρ
    sums.cov_w_thl_core += MH * w′ * θ_liq_ice′ * state.ρ
    sums.cov_qt_thl_core += MH * q_tot′ * θ_liq_ice′ * state.ρ
    sums.cov_qt_ei_core += MH * q_tot′ * e_int′ * state.ρ
    sums.cov_w_ei_core += MH * w′ * e_int′ * state.ρ

    return nothing
end

"""
    atmos_les_core_init(dgngrp, currtime)

Initialize the 'AtmosLESCore' diagnostics group.
"""
function atmos_les_core_init(dgngrp::DiagnosticsGroup, currtime)
    atmos = Settings.dg.balance_law
    FT = eltype(Settings.Q)
    mpicomm = Settings.mpicomm
    mpirank = MPI.Comm_rank(mpicomm)

    # FIXME properly
    if !isa(atmos.moisture, EquilMoist)
        @warn """
            Diagnostics $(dgngrp.name): can only be used with the `EquilMoist` moisture model
            """
        return nothing
    end

    atmos_collect_onetime(Settings.mpicomm, Settings.dg, Settings.Q)

    if mpirank == 0
        dims = OrderedDict("z" => (AtmosCollected.zvals, Dict()))

        # set up the variables we're going to be writing
        vars = OrderedDict()
        vars["core_frac"] = (("z",), FT, Dict())

        varnames = map(
            s -> startswith(s, "moisture.") ? s[10:end] : s,
            flattenednames(vars_atmos_les_core_simple(atmos, FT)),
        )
        ho_varnames = map(
            s -> startswith(s, "moisture.") ? s[10:end] : s,
            flattenednames(vars_atmos_les_core_ho(atmos, FT)),
        )
        append!(varnames, ho_varnames)
        for varname in varnames
            var = Variables[varname]
            vars[varname] = (("z",), FT, var.attrib)
        end

        # create the output file
        dprefix = @sprintf(
            "%s_%s_%s",
            dgngrp.out_prefix,
            dgngrp.name,
            Settings.starttime,
        )
        dfilename = joinpath(Settings.output_dir, dprefix)
        init_data(dgngrp.writer, dfilename, dims, vars)
    end

    return nothing
end

"""
    atmos_les_core_collect(dgngrp, currtime)

Perform a global grid traversal to compute various diagnostics.
"""
function atmos_les_core_collect(dgngrp::DiagnosticsGroup, currtime)
    mpicomm = Settings.mpicomm
    dg = Settings.dg
    Q = Settings.Q
    mpirank = MPI.Comm_rank(mpicomm)
    atmos = dg.balance_law
    if !isa(atmos.moisture, EquilMoist)
        @warn """
            Diagnostics $(dgngrp.name): can only be used with the `EquilMoist` moisture model
            """
        return nothing
    end
    grid = dg.grid
    grid_info = basic_grid_info(dg)
    topl_info = basic_topology_info(grid.topology)
    Nqk = grid_info.Nqk
    Nqh = grid_info.Nqh
    npoints = prod(grid_info.Nq)
    nrealelem = topl_info.nrealelem
    nvertelem = topl_info.nvertelem
    nhorzelem = topl_info.nhorzrealelem

    # get needed arrays onto the CPU
    if array_device(Q) isa CPU
        state_data = Q.realdata
        aux_data = dg.state_auxiliary.realdata
        vgeo = grid.vgeo
    else
        state_data = Array(Q.realdata)
        aux_data = Array(dg.state_auxiliary.realdata)
        vgeo = Array(grid.vgeo)
    end
    FT = eltype(state_data)

    zvals = AtmosCollected.zvals

    # Visit each node of the state variables array and:
    # - generate and store the thermo variables,
    # - if core condition holds (q_liq > 0 && w > 0)
    #   - count that point in the core fraction for that z
    #   - count the point's weighting towards averaging for that z, and
    #   - accumulate the simple horizontal sums
    #
    core_MH_z = zeros(FT, Nqk * nvertelem)
    thermo_array =
        [zeros(FT, num_thermo(atmos, FT)) for _ in 1:npoints, _ in 1:nrealelem]
    simple_sums = [
        zeros(FT, num_atmos_les_core_simple_vars(atmos, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    ql_w_gt_0 = [zeros(FT, (Nqh * nhorzelem)) for _ in 1:(Nqk * nvertelem)]
    @traverse_dg_grid grid_info topl_info begin
        state = extract_state(dg, state_data, ijk, e, Prognostic())
        aux = extract_state(dg, aux_data, ijk, e, Auxiliary())
        MH = vgeo[ijk, grid.MHid, e]

        thermo = thermo_vars(atmos, thermo_array[ijk, e])
        compute_thermo!(atmos, state, aux, thermo)

        if thermo.moisture.q_liq > 0 && state.ρu[3] > 0
            idx = (Nqh * (eh - 1)) + (grid_info.Nq[2] * (j - 1)) + i
            ql_w_gt_0[evk][idx] = one(FT)
            core_MH_z[evk] += MH

            simple = atmos_les_core_simple_vars(atmos, simple_sums[evk])
            atmos_les_core_simple_sums!(atmos, state, thermo, MH, simple)
        end
    end

    # reduce horizontal sums and core fraction across ranks and compute averages
    simple_avgs = [
        zeros(FT, num_atmos_les_core_simple_vars(atmos, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    core_frac = zeros(FT, Nqk * nvertelem)
    MPI.Allreduce!(core_MH_z, +, mpicomm)
    for evk in 1:(Nqk * nvertelem)
        tot_ql_w_gt_0 = MPI.Reduce(sum(ql_w_gt_0[evk]), +, 0, mpicomm)
        tot_horz = MPI.Reduce(length(ql_w_gt_0[evk]), +, 0, mpicomm)

        MPI.Allreduce!(simple_sums[evk], simple_avgs[evk], +, mpicomm)
        simple_avgs[evk] .= simple_avgs[evk] ./ core_MH_z[evk]

        if mpirank == 0
            core_frac[evk] = tot_ql_w_gt_0 / tot_horz
        end
    end

    # complete density averaging
    simple_varnames = flattenednames(vars_atmos_les_core_simple(atmos, FT))
    for vari in 1:length(simple_varnames)
        for evk in 1:(Nqk * nvertelem)
            simple_ha = atmos_les_core_simple_vars(atmos, simple_avgs[evk])
            avg_rho = simple_ha.avg_rho_core
            if simple_varnames[vari] != "avg_rho_core"
                simple_avgs[evk][vari] /= avg_rho
            end
        end
    end

    # compute the variances and covariances
    ho_sums = [
        zeros(FT, num_atmos_les_core_ho_vars(atmos, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    @traverse_dg_grid grid_info topl_info begin
        state = extract_state(dg, state_data, ijk, e, Prognostic())
        thermo = thermo_vars(atmos, thermo_array[ijk, e])
        MH = vgeo[ijk, grid.MHid, e]

        if thermo.moisture.q_liq > 0 && state.ρu[3] > 0
            simple_ha = atmos_les_core_simple_vars(atmos, simple_avgs[evk])
            ho = atmos_les_core_ho_vars(atmos, ho_sums[evk])
            atmos_les_core_ho_sums!(atmos, state, thermo, MH, simple_ha, ho)
        end
    end

    # reduce across ranks and compute averages
    ho_avgs = [
        zeros(FT, num_atmos_les_core_ho_vars(atmos, FT))
        for _ in 1:(Nqk * nvertelem)
    ]
    for evk in 1:(Nqk * nvertelem)
        MPI.Reduce!(ho_sums[evk], ho_avgs[evk], +, 0, mpicomm)
        if mpirank == 0
            ho_avgs[evk] .= ho_avgs[evk] ./ core_MH_z[evk]
        end
    end

    # complete density averaging and prepare output
    if mpirank == 0
        varvals = OrderedDict()
        varvals["core_frac"] = core_frac

        for (vari, varname) in enumerate(simple_varnames)
            davg = zeros(FT, Nqk * nvertelem)
            for evk in 1:(Nqk * nvertelem)
                davg[evk] = simple_avgs[evk][vari]
            end
            varvals[varname] = davg
        end

        ho_varnames = flattenednames(vars_atmos_les_core_ho(atmos, FT))
        for (vari, varname) in enumerate(ho_varnames)
            davg = zeros(FT, Nqk * nvertelem)
            for evk in 1:(Nqk * nvertelem)
                simple_ha = atmos_les_core_simple_vars(atmos, simple_avgs[evk])
                avg_rho = simple_ha.avg_rho_core
                davg[evk] = ho_avgs[evk][vari] / avg_rho
            end
            varvals[varname] = davg
        end

        # write output
        append_data(dgngrp.writer, varvals, currtime)
    end

    MPI.Barrier(mpicomm)
    return nothing
end # function collect

function atmos_les_core_fini(dgngrp::DiagnosticsGroup, currtime) end
