# Copy the specified kind of state's variables into an `MArray`.
# This is GPU-safe.
function extract_state(bl, state, ijk, e, st::AbstractStateType)
    FT = eltype(state)
    num_state = number_states(bl, st)
    local_state = MArray{Tuple{num_state}, FT}(undef)
    for s in 1:num_state
        local_state[s] = state[ijk, s, e]
    end
    return Vars{vars_state(bl, st, FT)}(local_state)
end

# Return `true` if the specified symbol is a type name that is a subtype
# of `BalanceLaw` and `false` otherwise.
isa_bl(sym::Symbol) = any(
    bl -> endswith(bl, "." * String(sym)),
    map(bl -> String(Symbol(bl)), subtypes(BalanceLaw)),
)
isa_bl(ex) = false

uppers_in(s) = foldl((f, c) -> isuppercase(c) ? f * c : f, s, init = "")

# Generate the common definitions used in many places.
function generate_common_defs()
    quote
        mpicomm = Settings.mpicomm
        dg = Settings.dg
        Q = Settings.Q
        mpirank = MPI.Comm_rank(mpicomm)
        bl = dg.balance_law
        grid = dg.grid
        grid_info = basic_grid_info(dg)
        topl_info = basic_topology_info(grid.topology)
        Nq1 = grid_info.Nq[1]
        Nq2 = grid_info.Nq[2]
        Nqh = grid_info.Nqh
        Nqk = grid_info.Nqk
        npoints = prod(grid_info.Nq)
        nrealelem = topl_info.nrealelem
        nvertelem = topl_info.nvertelem
        nhorzelem = topl_info.nhorzrealelem
        FT = eltype(Q)
        interpol = dgngrp.interpol
        params = dgngrp.params
    end
end

# Generate the `dims` dictionary for `Writers.init_data`.
function generate_init_dims(::NoInterpolation, cfg, dvtype_dvars_map)
    dimslst = Any[]
    for dvtype in keys(dvtype_dvars_map)
        dimnames = DiagnosticsMachine.dv_dg_dimnames(cfg, dvtype)
        dimranges = DiagnosticsMachine.dv_dg_dimranges(cfg, dvtype)
        for (dimname, dimrange) in zip(dimnames, dimranges)
            lhs = :($dimname)
            rhs = :(collect($dimrange), Dict())
            push!(dimslst, :($lhs => $rhs))
        end
    end

    quote
        OrderedDict($(Expr(:tuple, dimslst...))...)
    end
end
function generate_init_dims(::InterpolationType, cfg, dvtype_dvars_map)
    quote
        dims = dimensions(interpol)
        if interpol isa InterpolationCubedSphere
            # Adjust `level` on the sphere.
            level_val = dims["level"]
            dims["level"] = (
                level_val[1] .- FT(planet_radius(Settings.param_set)),
                level_val[2],
            )
        end
        dims
    end
end

get_dimnames(::NoInterpolation, cfg, dvtype) =
    DiagnosticsMachine.dv_dg_dimnames(cfg, dvtype)
get_dimnames(::InterpolationType, cfg, dvtype) = :(tuple(collect(keys(dims))...))

# Generate the `vars` dictionary for `Writers.init_data`.
function generate_init_vars(intrp, cfg, dvtype_dvars_map)
    varslst = Any[]
    for (dvtype, dvlst) in dvtype_dvars_map
        for dvar in dvlst
            lhs = :($(DiagnosticsMachine.dv_name(cfg, dvar)))
            dimnames = get_dimnames(intrp, cfg, dvtype)
            rhs = :($dimnames, FT, $(DiagnosticsMachine.dv_attrib(cfg, dvar)))
            push!(varslst, :($lhs => $rhs))
        end
    end

    quote
        # TODO: add code to filter this based on what's actually in `bl`.
        OrderedDict($(Expr(:tuple, varslst...))...)
    end
end

# Generate `Diagnostics.$(name)_init(...)` which will initialize the
# `DiagnosticsGroup` when called.
function generate_init_fun(
    name,
    config_type,
    params_type,
    init_fun,
    interpolate,
    dvtype_dvars_map,
)
    init_name = Symbol(name, "_init")
    cfg = getfield(ConfigTypes, config_type)
    intrp = getfield(@__MODULE__, interpolate)
    quote
        function $init_name(dgngrp, curr_time)
            $(generate_common_defs())

            $(init_fun)(dgngrp, curr_time)

            # temporarily if dgngrp.onetime
            DiagnosticsMachine.collect_onetime(mpicomm, dg, Q)
            # temporarily end

            if mpirank == 0
                dims = $(generate_init_dims(intrp(), cfg(), dvtype_dvars_map))
                vars = $(generate_init_vars(intrp(), cfg(), dvtype_dvars_map))

                # create the output file
                dprefix = @sprintf("%s_%s", dgngrp.out_prefix, dgngrp.name,)
                dfilename = joinpath(Settings.output_dir, dprefix)
                init_data(dgngrp.writer, dfilename, dims, vars)
            end

            return nothing
        end
    end
end

# Generate code snippet for copying arrays to the CPU if needed. Ideally,
# this will be removed when diagnostics are made to run on GPU.
function generate_array_copies()
    quote
        # get needed arrays onto the CPU
        if array_device(Q) isa CPU
            state_data = Q.realdata
            gradflux_data = dg.state_gradient_flux.realdata
            aux_data = dg.state_auxiliary.realdata
            vgeo = grid.vgeo
        else
            state_data = Array(Q.realdata)
            gradflux_data = Array(dg.state_gradient_flux.realdata)
            aux_data = Array(dg.state_auxiliary.realdata)
            vgeo = Array(grid.vgeo)
        end
    end
end

# Generate code to create the necessary arrays for the diagnostics
# variables.
function generate_create_vars_arrays(
    ::CollectOnInterpolatedGrid,
    cfg,
    dvtype_dvars_map,
)
    quote end
end
function generate_create_vars_arrays(::InterpolationType, cfg, dvtype_dvars_map)
    cva_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        dvt_short = uppers_in(split(String(Symbol(dvtype)), ".")[end])
        arr_name = Symbol("vars_", dvt_short, "_array")
        npoints = DiagnosticsMachine.dv_dg_points_length(cfg, dvtype)
        nvars = length(dvlst)
        nelems = DiagnosticsMachine.dv_dg_elems_length(cfg, dvtype)
        cva_ex = quote
            $arr_name = Array{FT}(undef, $npoints, $nvars, $nelems)
            fill!($arr_name, 0)
        end
        push!(cva_exs, cva_ex)
    end
    return Expr(:block, (cva_exs...))
end

# Generate calls to the implementations for the `DiagnosticVar`s in this
# group and store the results.
function generate_collect_calls(::InterpolationType, cfg, dvtype_dvars_map)
    cc_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        dvt = split(String(Symbol(dvtype)), ".")[end]
        dvt_short = uppers_in(dvt)
        arr_name = Symbol("vars_", dvt_short, "_array")
        pt = DiagnosticsMachine.dv_dg_points_index(cfg, dvtype)
        elem = DiagnosticsMachine.dv_dg_elems_index(cfg, dvtype)
        var_impl = Symbol("dv_", dvt)

        for (v, dvar) in enumerate(dvlst)
            impl_args = DiagnosticsMachine.dv_args(cfg, dvar)
            AT1 = impl_args[1][2] # the type of the first argument
            if isa_bl(AT1)
                impl_extra_params = ()
            else
                AT2 = impl_args[2][2] # the type of the second argument
                @assert isa_bl(AT2)
                AN1 = impl_args[1][1] # the name of the first argument
                impl_extra_params = (:(getproperty(bl, $(QuoteNode(AN1)))),)
            end
            cc_ex = DiagnosticsMachine.dv_op(
                cfg,
                dvtype,
                :($arr_name[$pt, $v, $elem]),
                :(DiagnosticsMachine.$(var_impl)(
                    $cfg,
                    $dvar,
                    $(impl_extra_params...),
                    bl,
                    states,
                    curr_time,
                    cache,
                )),
            )
            push!(cc_exs, cc_ex)
        end
    end

    return Expr(:block, (cc_exs...))
end

# Generate the nested loops to traverse the DG grid within which we extract
# the various states and then generate the individual collection calls.
function generate_dg_collections(
    ::CollectOnInterpolatedGrid,
    cfg,
    dvtype_dvars_map,
)
    quote end
end
function generate_dg_collections(
    intrp::InterpolationType,
    cfg,
    dvtype_dvars_map,
)
    quote
        cache = Dict{Symbol, Any}()
        for eh in 1:nhorzelem, ev in 1:nvertelem
            e = ev + (eh - 1) * nvertelem
            for k in 1:Nqk, j in 1:Nq2, i in 1:Nq1
                ijk = i + Nq1 * ((j - 1) + Nq2 * (k - 1))
                MH = vgeo[ijk, grid.MHid, e]
                state = DiagnosticsMachine.extract_state(
                    bl,
                    state_data,
                    ijk,
                    e,
                    Prognostic(),
                )
                gradflux = DiagnosticsMachine.extract_state(
                    bl,
                    gradflux_data,
                    ijk,
                    e,
                    GradientFlux(),
                )
                aux = DiagnosticsMachine.extract_state(
                    bl,
                    aux_data,
                    ijk,
                    e,
                    Auxiliary(),
                )
                states = States(state, gradflux, aux)
                $(generate_collect_calls(intrp, cfg, dvtype_dvars_map))
            end
        end
        empty!(cache)
    end
end

# Generate any reductions needed for the data collected thus far.
function generate_dg_reductions(::NoInterpolation, cfg, dvtype_dvars_map)
    red_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        dvt_short = uppers_in(split(String(Symbol(dvtype)), ".")[end])
        arr_name = Symbol("vars_", dvt_short, "_array")
        red_ex = DiagnosticsMachine.dv_reduce(cfg, dvtype, arr_name)
        push!(red_exs, red_ex)
    end

    return Expr(:block, (red_exs...))
end
function generate_dg_reductions(::InterpolationType, cfg, dvtype_dvars_map)
    quote end
end

# Generate interpolation calls as needed. None for `NoInterpolation`.
function generate_interpolations(::NoInterpolation, cfg, dvtype_dvars_map)
    quote end
end
# Interpolate only the diagnostic variables arrays.
function generate_interpolations(
    ::InterpolateAfterCollection,
    cfg,
    dvtype_dvars_map,
)
    ic_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        nvars = length(dvlst)
        dvt = split(String(Symbol(dvtype)), ".")[end]
        dvt_short = uppers_in(dvt)
        arr_name = Symbol("vars_", dvt_short, "_array")
        iarr_name = Symbol("i", arr_name)
        acc_arr_name = Symbol("acc_", arr_name)
        ic_ex = quote
            # TODO: `interpolate_local!` expects arrays of the same type
            # as `grid`. This should be solved when we make this GPU-capable.
            if !(array_device(Q) isa CPU)
                ArrayType = arraytype(grid)
                $arr_name = ArrayType($arr_name)
            end
            $iarr_name = similar($arr_name, interpol.Npl, $nvars)
            interpolate_local!(interpol, $arr_name, $iarr_name)
            # TODO: projection
            $acc_arr_name =
                accumulate_interpolated_data(mpicomm, interpol, $iarr_name)
            if !(array_device(Q) isa CPU)
                $acc_arr_name = Array($acc_arr_name)
            end
        end
        push!(ic_exs, ic_ex)
    end

    return Expr(:block, (ic_exs...))
end
# Interpolate all the arrays needed for `States`.
function generate_interpolations(
    ::CollectOnInterpolatedGrid,
    cfg,
    dvtype_dvars_map,
)
    quote
        istate_array =
            similar(Q.realdata, interpol.Npl, number_states(bl, Prognostic()))
        interpolate_local!(interpol, Q.realdata, istate_array)
        igradflux_array =
            similar(Q.realdata, interpol.Npl, number_states(bl, GradientFlux()))
        interpolate_local!(
            interpol,
            dg.state_gradient_flux.realdata,
            igradflux_array,
        )
        iaux_array =
            similar(Q.realdata, interpol.Npl, number_states(bl, Auxiliary()))
        interpolate_local!(interpol, dg.state_auxiliary.realdata, iaux_array)

        i_ρu = varsindex(vars_state(bl, Prognostic(), FT), :ρu)
        project_cubed_sphere!(interpol, istate_array, tuple(collect(i_ρu)...))

        # FIXME: accumulating to rank 0 is not scalable
        all_state_data =
            accumulate_interpolated_data(mpicomm, interpol, istate_array)
        all_gradflux_data =
            accumulate_interpolated_data(mpicomm, interpol, igradflux_array)
        all_aux_data =
            accumulate_interpolated_data(mpicomm, interpol, iaux_array)
    end
end

# Generate code to create the necessary arrays to collect the diagnostics
# variables on the interpolated grid.
function generate_create_i_vars_arrays(
    ::CollectOnInterpolatedGrid,
    cfg,
    dvtype_dvars_map,
)
    cva_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        dvt_short = uppers_in(split(String(Symbol(dvtype)), ".")[end])
        acc_arr_name = Symbol("acc_vars_", dvt_short, "_array")
        nvars = length(dvlst)
        cva_ex = quote
            $acc_arr_name = Array{FT}(undef, nx, ny, nz, $nvars)
            fill!($acc_arr_name, 0)
        end
        push!(cva_exs, cva_ex)
    end

    return Expr(:block, (cva_exs...))
end
function generate_create_i_vars_arrays(
    ::InterpolationType,
    cfg,
    dvtype_dvars_map,
)
    quote end
end

# Generate the nested loops to traverse the interpolated grid within
# which we extract the various (interpolated) states and then generate
# the individual collection calls.
function generate_i_collections(
    intrp::CollectOnInterpolatedGrid,
    cfg,
    dvtype_dvars_map,
)
    quote
        (x1, x2, x3) = map(k -> length(dims[k][1]), collect(keys(dims)))
        for x in 1:x1, y in 1:x2, z in 1:x3
            istate = Vars{vars_state(bl, Prognostic(), FT)}(view(
                all_state_data,
                x,
                y,
                z,
                :,
            ))
            igradflux = Vars{vars_state(bl, GradientFlux(), FT)}(view(
                all_gradflux_data,
                x,
                y,
                z,
                :,
            ))
            iaux = Vars{vars_state(bl, Auxiliary(), FT)}(view(
                all_aux_data,
                x,
                y,
                z,
                :,
            ))
            states = States(istate, igradflux, iaux)
            $(generate_collect_calls(intrp, cfg, dvtype_dvars_map))
        end
    end
end
function generate_i_collections(::InterpolationType, cfg, dvtype_dvars_map)
    quote end
end

# Generate assignments into `varvals` for writing.
function generate_varvals(::NoInterpolation, cfg, dvtype_dvars_map)
    vv_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        dvt = split(String(Symbol(dvtype)), ".")[end]
        dvt_short = uppers_in(dvt)
        arr_name = Symbol("vars_", dvt_short, "_array")
        for (v, dvar) in enumerate(dvlst)
            vv_ex = quote
                varvals[$(DiagnosticsMachine.dv_name(cfg, dvar))] =
                    reshape(view($(arr_name), :, $v, :), :)
            end
            push!(vv_exs, vv_ex)
        end
    end

    return Expr(:block, (vv_exs...))
end
function generate_varvals(::InterpolationType, cfg, dvtype_dvars_map)
    vv_exs = []
    for (dvtype, dvlst) in dvtype_dvars_map
        dvt = split(String(Symbol(dvtype)), ".")[end]
        dvt_short = uppers_in(dvt)
        acc_arr_name = Symbol("acc_vars_", dvt_short, "_array")
        for (v, dvar) in enumerate(dvlst)
            vv_ex = quote
                varvals[$(DiagnosticsMachine.dv_name(cfg, dvar))] =
                    $(acc_arr_name)[:, :, :, $v]
            end
            push!(vv_exs, vv_ex)
        end
    end

    return Expr(:block, (vv_exs...))
end

# Generate `Diagnostics.$(name)_collect(...)` which when called,
# performs a collection of all the diagnostic variables in the group
# and writes them out.
function generate_collect_fun(
    name,
    config_type,
    params_type,
    init_fun,
    interpolate,
    dvtype_dvars_map,
)
    collect_name = Symbol(name, "_collect")
    cfg = getfield(ConfigTypes, config_type)
    intrp = getfield(@__MODULE__, interpolate)
    quote
        function $collect_name(dgngrp, curr_time)
            $(generate_common_defs())
            $(generate_array_copies())
            $(generate_create_vars_arrays(intrp(), cfg(), dvtype_dvars_map))

            # Traverse the DG grid and collect diagnostics as needed.
            $(generate_dg_collections(intrp(), cfg(), dvtype_dvars_map))

            # Perform any reductions necessary.
            $(generate_dg_reductions(intrp(), cfg(), dvtype_dvars_map))

            # TODO: density averaging.

            # Interpolate and accumulate if needed.
            $(generate_interpolations(intrp(), cfg(), dvtype_dvars_map))

            if mpirank == 0
                dims = dimensions(interpol)

                $(generate_create_i_vars_arrays(
                    intrp(),
                    cfg(),
                    dvtype_dvars_map,
                ))

                # Traverse the interpolated grid and collect diagnostics if needed.
                $(generate_i_collections(intrp(), cfg(), dvtype_dvars_map))

                # Assemble the diagnostic variables and write them.
                varvals = OrderedDict()
                $(generate_varvals(intrp(), cfg(), dvtype_dvars_map))
                append_data(dgngrp.writer, varvals, curr_time)
            end

            MPI.Barrier(mpicomm)
            return nothing
        end
    end
end

# Generate `Diagnostics.$(name)_fini(...)`, which does nothing right now.
function generate_fini_fun(name, args...)
    fini_name = Symbol(name, "_fini")
    quote
        function $fini_name(dgngrp, curr_time) end
    end
end

# Generate `$(name)(...)` which will create the `DiagnosticsGroup` for 
# $name when called.
function generate_setup_fun(
    name,
    config_type,
    params_type,
    init_fun,
    interpolate,
    dvtype_dvars_map,
)
    init_name = Symbol(name, "_init")
    collect_name = Symbol(name, "_collect")
    fini_name = Symbol(name, "_fini")

    setup_name = Symbol(name)
    intrp = getfield(@__MODULE__, interpolate)

    no_intrp_err = quote end
    some_intrp_err = quote end
    if intrp() isa NoInterpolation
        no_intrp_err = quote
            @warn "$($name) does not specify interpolation, but an " *
                  "`InterpolationTopology` has been provided; ignoring."
            interpol = nothing
        end
    else
        some_intrp_err = quote
            throw(ArgumentError(
                "$($name) specifies interpolation, but no " *
                "`InterpolationTopology` has been provided.",
            ))
        end
    end
    quote
        function $setup_name(
            interval::String,
            out_prefix::String,
            params::$params_type = nothing;
            writer = NetCDFWriter(),
            interpol = nothing,
        ) where {
            $config_type <: ClimateMachineConfigType,
            $params_type <: Union{Nothing, DiagnosticsGroupParams},
        }
            if isnothing(interpol)
                $(some_intrp_err)
            else
                $(no_intrp_err)
            end

            return DiagnosticsMachine.DiagnosticsGroup(
                $name,
                $init_name,
                $fini_name,
                $collect_name,
                interval,
                out_prefix,
                writer,
                interpol,
                # temporarily $(intrp() isa NoInterpolation),
                params,
            )
        end
    end
end
