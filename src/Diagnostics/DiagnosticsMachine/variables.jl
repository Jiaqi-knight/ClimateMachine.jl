"""
    DiagnosticVar

The base type for diagnostic variables.
"""
abstract type DiagnosticVar end

"""
    dv_name(::CT, ::DVT) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }
"""
function dv_name end

"""
    dv_attrib(::CT, ::DVT) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }
"""
function dv_attrib end

# Default method for variable attributes.
dv_attrib(::ClimateMachineConfigType, ::DiagnosticVar) = Dict()

"""
    dv_args(::CT, ::DVT) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns a tuple of the arguments specified by the implementation of the
diagnostic variables.
"""
function dv_args end

"""
    dv_dg_points_length(::CT, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns an expression that evaluates to the length of the points
dimension of a DG array of the diagnostic variable type.
"""
function dv_dg_points_length end

"""
    dv_dg_points_index(::CT, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns an expression that, within a `@traverse_dg_grid`, evaluates
to an index into the points dimension of a DG array of the diagnostic
variable type.
"""
function dv_dg_points_index end

"""
    dv_dg_elems_length(::CT, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns an expression that evaluates to the length of the elements
dimension of a DG array of the diagnostic variable type.
"""
function dv_dg_elems_length end

"""
    dv_dg_elems_index(::CT, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns an expression that, within a `@traverse_dg_grid`, evaluates
to an index into the elements dimension of a DG array of the diagnostic
variable type.
"""
function dv_dg_elems_index end

"""
    dv_dg_dimnames(::CT, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns a tuple of the names of the dimensions of the diagnostic
variable type.
"""
function dv_dg_dimnames end

"""
    dv_dg_dimranges(::CT, ::Type{DVT}) where {
        CT <: ClimateMachineConfigType,
        DVT <: DiagnosticVar,
    }

Returns a tuple of expressions that evaluate to the range of the
dimensions for the diagnostic variable type.
"""
function dv_dg_dimranges end

# Default method for variable dimension ranges.
function dv_dimranges(
    ::ClimateMachineConfigType,
    ::Type{DiagnosticVar},
    ::Nothing,
)
    (:(1:npoints), :(1:nrealelem))
end

"""
"""
function dv_op end

"""
"""
function dv_reduce end

####
#
# Generate a standardized type name from the diagnostic variable name.
function dv_type_name(dvtype, config_type, name)
    let uppers_in(s) =
            foldl((f, c) -> isuppercase(c) ? f * c : f, String(s), init = "")
        return uppers_in(config_type) *
               "_" *
               uppers_in(dvtype) *
               "_" *
               String(name)
    end
end

# Generate the type and interface functions for a diagnostic variable.
function generate_dv_interface(
    dvtype,
    config_type,
    name,
    units = "",
    long_name = "",
    standard_name = "",
)
    dvtypname = Symbol(dv_type_name(dvtype, config_type, name))
    attrib_ex = quote end
    if any(a -> a != "", [units, long_name, standard_name])
        attrib_ex = quote
            dv_attrib(::$config_type, ::$dvtypname) = OrderedDict(
                "units" => $units,
                "long_name" => $long_name,
                "standard_name" => $standard_name,
            )
        end
    end
    quote
        struct $dvtypname <: $dvtype end
        DiagnosticsMachine.AllDiagnosticVars[$config_type][$(String(name))] =
            $dvtypname()
        dv_name(::$config_type, ::$dvtypname) = $(String(name))
        $(attrib_ex)
    end
end

# Helper to generate the implementation function for a diagnostic variable.
function generate_dv_function(dvtype, config_type, name, impl)
    dvfun = Symbol("dv_", dvtype)
    dvtypname = Symbol(dv_type_name(dvtype, config_type, name))
    @capture(impl, ((args__,),) -> (body_)) ||
        @capture(impl, (args_) -> (body_)) ||
        error("Bad implementation for $(esc(names[1]))")
    split_fun_args = map(splitarg, args)
    fun_args = map(a -> :($(a[1])::$(a[2])), split_fun_args)
    quote
        function dv_args(::$config_type, ::$dvtypname)
            $split_fun_args
        end
        function $dvfun(::$config_type, ::$dvtypname, $(fun_args...))
            $body
        end
    end
end

"""
    States

Composite of the various states, used as a parameter to diagnostic
collection functions.
"""
struct States{PS, GFS, AS}
    prognostic::PS
    gradient_flux::GFS
    auxiliary::AS
end

# Diagnostic variable types and interfaces to create diagnostic variables
# of these types.
include("pointwise.jl")
include("horizontal_average.jl")
include("scalar.jl")
