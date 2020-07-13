"""
    ScalarDiagnostic

A reduction into a scalar value.
"""
abstract type ScalarDiagnostic <: DiagnosticVar end
dv_ScalarDiagnostic(
    ::ClimateMachineConfigType,
    ::ScalarDiagnostic,
    ::BalanceLaw,
    ::States,
    ::AbstractFloat,
) = nothing

function dv_dg_points_length(
    ::ClimateMachineConfigType,
    ::Type{ScalarDiagnostic},
)
    1
end
function dv_dg_points_index(
    ::ClimateMachineConfigType,
    ::Type{ScalarDiagnostic},
)
    1
end

function dv_dg_elems_length(
    ::ClimateMachineConfigType,
    ::Type{ScalarDiagnostic},
)
    1
end
function dv_dg_elems_index(::ClimateMachineConfigType, ::Type{ScalarDiagnostic})
    1
end

function dv_dg_dimnames(::ClimateMachineConfigType, ::Type{ScalarDiagnostic})
    ()
end
function dv_dg_dimranges(::ClimateMachineConfigType, ::Type{ScalarDiagnostic})
    (1, 1)
end

dv_dimnames(::ClimateMachineConfigType, ::Type{ScalarDiagnostic}, ::Any) = ()

macro scalar_diagnostic(impl, config_type, name)
    iex = quote
        $(generate_dv_interface(:ScalarDiagnostic, config_type, name))
        $(generate_dv_function(:ScalarDiagnostic, config_type, name, impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end

macro scalar_diagnostic(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :ScalarDiagnostic,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(:ScalarDiagnostic, config_type, [name], impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end
