using ..Atmos
using ..TurbulenceClosures
using ..TurbulenceConvection

# Method definitions for diagnostics collection for all the components
# of `AtmosModel`.

dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::PointwiseDiagnostic,
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::HorizontalAverage,
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::ScalarDiagnostic,
    ::MoistureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0

dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::PointwiseDiagnostic,
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::HorizontalAverage,
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::ScalarDiagnostic,
    ::PrecipitationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0

dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::PointwiseDiagnostic,
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::HorizontalAverage,
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::ScalarDiagnostic,
    ::RadiationModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0

dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::PointwiseDiagnostic,
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::HorizontalAverage,
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::ScalarDiagnostic,
    ::TracerModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0

dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::PointwiseDiagnostic,
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::HorizontalAverage,
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::ScalarDiagnostic,
    ::TurbulenceClosureModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0

dv_PointwiseDiagnostic(
    ::AtmosConfigType,
    ::PointwiseDiagnostic,
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_HorizontalAverage(
    ::AtmosConfigType,
    ::HorizontalAverage,
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
dv_ScalarDiagnostic(
    ::AtmosConfigType,
    ::ScalarDiagnostic,
    ::TurbulenceConvectionModel,
    ::AtmosModel,
    ::States,
    ::AbstractFloat,
    ::Dict{Symbol, Any},
) = 0
