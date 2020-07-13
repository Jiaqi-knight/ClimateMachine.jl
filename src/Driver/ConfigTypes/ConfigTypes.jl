"""
    ConfigTypes

Module containing ClimateMachine configuration types.
"""
module ConfigTypes

export ClimateMachineConfigType,
    AtmosConfigType,
    AtmosLESConfigType,
    AtmosGCMConfigType,
    OceanConfigType,
    OceanBoxGCMConfigType,
    OceanSplitExplicitConfigType,
    SingleStackConfigType

abstract type ClimateMachineConfigType end
abstract type AtmosConfigType <: ClimateMachineConfigType end
struct AtmosLESConfigType <: AtmosConfigType end
struct AtmosGCMConfigType <: AtmosConfigType end
abstract type OceanConfigType <: ClimateMachineConfigType end
struct OceanBoxGCMConfigType <: OceanConfigType end
struct OceanSplitExplicitConfigType <: OceanConfigType end
struct SingleStackConfigType <: ClimateMachineConfigType end

end
