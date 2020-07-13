using ..Atmos
using ..ConfigTypes
using ..DiagnosticsMachine

@diagnostics_group(
    "AtmosGCMDefault",
    AtmosGCMConfigType,
    Nothing,
    (_...) -> nothing,
    # CollectOnInterpolatedGrid | InterpolateAfterCollection | NoInterpolation
    #NoInterpolation,
    InterpolateAfterCollection,
    #CollectOnInterpolatedGrid,
    # simple horizontal averages
    u,
    v,
    w,
    rho,
    temp,
    pres,
    thd,
    et,
    ei,
    ht,
    hi,
    #vort, TODO
    # moisture related
    qt,
    ql,
    qv,
    qi,
    thv,
    thl,
)
