"""
    StdDiagnostics

This module defines many standard diagnostic variables and groups that may
be used directly by experiments.
"""
module StdDiagnostics

using KernelAbstractions
using MPI
using OrderedCollections
using Printf

using ..Diagnostics # temporarily
using ..Atmos
using ..ConfigTypes
using ..DGMethods
using ..DiagnosticsMachine
import ..DiagnosticsMachine: Settings, dv_name, dv_attrib, dv_args
using ..Mesh.Interpolation
using ..Mesh.Topologies
using ..VariableTemplates
using ..Writers

export setup_atmos_default_diagnostics


# Pre-defined diagnostic variables

# Atmos
include("atmos_les_diagnostic_vars.jl")
include("atmos_gcm_diagnostic_vars.jl")


# Pre-defined diagnostics groups

# Atmos
include("atmos_les_default.jl")
include("atmos_gcm_default.jl")

end # module StdDiagnostics
