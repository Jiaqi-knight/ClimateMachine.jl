using ..Atmos
using ..ConfigTypes
using ..DiagnosticsMachine

@diagnostics_group(
    # name
    "AtmosLESDefault",
    # configuration type
    AtmosLESConfigType,
    # params type
    Nothing,
    # user initialization function
    (_...) -> nothing,
    # CollectOnInterpolatedGrid | InterpolateAfterCollection | NoInterpolation
    NoInterpolation,
    #InterpolateAfterCollection,
    #CollectOnInterpolatedGrid,
    # simple horizontal averages
    u,
    v,
    w,
    avg_rho,
    rho,
    temp,
    pres,
    thd,
    et,
    ei,
    ht,
    hi,
    w_ht_sgs,
    # moisture related
    qt,
    ql,
    qi,
    qv,
    thv,
    thl,
    w_qt_sgs,
    # for variances and co-variances
    uu,
    vv,
    ww,
    www,
    eiei,
    wu,
    wv,
    wrho,
    wthd,
    wei,
    qtqt,
    thlthl,
    wqt,
    wql,
    wqi,
    wqv,
    wthv,
    wthl,
    qtthl,
    qtei,
    #cld_top, TODO
    #cld_base, TODO
    #cld_cover, TODO
    #lwp, TODO
    #rwp, TODO
)
