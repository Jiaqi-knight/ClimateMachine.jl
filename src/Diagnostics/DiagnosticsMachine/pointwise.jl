"""
    PointwiseDiagnostic

A diagnostic with the same dimensions as the original grid (DG or interpolated).
"""
abstract type PointwiseDiagnostic <: DiagnosticVar end
function dv_PointwiseDiagnostic end

function dv_dg_points_length(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
)
    :(npoints)
end
function dv_dg_points_index(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
)
    :(ijk)
end

function dv_dg_elems_length(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
)
    :(nrealelem)
end
function dv_dg_elems_index(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
)
    :(e)
end

function dv_dg_dimnames(::ClimateMachineConfigType, ::Type{PointwiseDiagnostic})
    ("nodes", "elements")
end
function dv_dg_dimranges(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
)
    (:(1:npoints), :(1:nrealelem))
end

function dv_op(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
    x,
    y,
    scale_with = 1,
)
    x = y
end

# Reduction for point-wise diagnostics would be a gather, but that will probably
# blow up memory. TODO.
function dv_reduce(
    ::ClimateMachineConfigType,
    ::Type{PointwiseDiagnostic},
    array_name,
)
    quote end
end

macro pointwise_diagnostic(impl, config_type, name)
    iex = quote
        $(generate_dv_interface(:PointwiseDiagnostic, config_type, name))
        $(generate_dv_function(:PointwiseDiagnostic, config_type, name, impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end

"""
    @pointwise_diagnostic(
        impl,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )

Define a point-wise diagnostic variable.
"""
macro pointwise_diagnostic(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :PointwiseDiagnostic,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(:PointwiseDiagnostic, config_type, name, impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end
