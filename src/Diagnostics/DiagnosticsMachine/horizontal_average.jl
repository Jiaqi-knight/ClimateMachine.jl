"""
    HorizontalAverage

A horizontal reduction into a single vertical dimension.
"""
abstract type HorizontalAverage <: DiagnosticVar end
function dv_HorizontalAverage end

# replace these with a `dv_array_dims` that takes `nvars` and returns the dims for the array
# or create the array? Use `Array`? `similar`?
function dv_dg_points_length(
    ::ClimateMachineConfigType,
    ::Type{HorizontalAverage},
)
    :(Nqk)
end
function dv_dg_points_index(
    ::ClimateMachineConfigType,
    ::Type{HorizontalAverage},
)
    :(k)
end

function dv_dg_elems_length(
    ::ClimateMachineConfigType,
    ::Type{HorizontalAverage},
)
    :(nvertelem)
end
function dv_dg_elems_index(
    ::ClimateMachineConfigType,
    ::Type{HorizontalAverage},
)
    :(ev)
end

function dv_dg_dimnames(::ClimateMachineConfigType, ::Type{HorizontalAverage})
    ("z",)
end
function dv_dg_dimranges(::ClimateMachineConfigType, ::Type{HorizontalAverage})
    z = quote
        ijk_range = 1:Nqh:npoints
        e_range = 1:nhorzelem:nrealelem
        reshape(grid.vgeo[ijk_range, grid.x3id, e_range], :)
    end
    (z,)
end

function dv_op(::ClimateMachineConfigType, ::Type{HorizontalAverage}, lhs, rhs)
    :($lhs += MH * $rhs)
end
function dv_reduce(
    ::ClimateMachineConfigType,
    ::Type{HorizontalAverage},
    array_name,
)
    quote
        MPI.Reduce!($array_name, +, 0, mpicomm)
        if mpirank == 0
            for v in 1:size($array_name, 2)
                $(array_name)[:, v, :] ./= DiagnosticsMachine.Collected.ΣMH_z
            end
        end
    end
end

macro horizontal_average(impl, config_type, name)
    iex = quote
        $(generate_dv_interface(:HorizontalAverage, config_type, name))
        $(generate_dv_function(:HorizontalAverage, config_type, name, impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end

"""
    @horizontal_average(
        impl,
        config_type,
        name,
        units,
        long_name,
        standard_name,
    )

Define a horizontal average diagnostic variable.
"""
macro horizontal_average(
    impl,
    config_type,
    name,
    units,
    long_name,
    standard_name,
)
    iex = quote
        $(generate_dv_interface(
            :HorizontalAverage,
            config_type,
            name,
            units,
            long_name,
            standard_name,
        ))
        $(generate_dv_function(:HorizontalAverage, config_type, name, impl))
    end
    esc(MacroTools.prewalk(unblock, iex))
end
