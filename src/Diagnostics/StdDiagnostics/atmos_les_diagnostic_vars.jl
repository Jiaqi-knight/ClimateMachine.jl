@horizontal_average(
    AtmosLESConfigType,
    u,
    "m s^-1",
    "x-velocity",
    "",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[1]
end

@horizontal_average(
    AtmosLESConfigType,
    v,
    "m s^-1",
    "y-velocity",
    "",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[2]
end

@horizontal_average(
    AtmosLESConfigType,
    w,
    "m s^-1",
    "z-velocity",
    "",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[3]
end

@horizontal_average(
    AtmosLESConfigType,
    uu,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[1]^2 / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    vv,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[2]^2 / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    ww,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[3]^2 / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    www,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[3]^3 / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    wu,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[3] * states.prognostic.ρu[1] / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    wv,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[3] * states.prognostic.ρu[2] / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    avg_rho,
    "kg m^-3",
    "air density",
    "air_density",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    rho,
    "kg m^-3",
    "air density",
    "air_density",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρ^2
end

@horizontal_average(
    AtmosLESConfigType,
    "wrho",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρu[3] * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    temp,
    "K",
    "air temperature",
    "air_temperature",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    air_temperature(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    pres,
    "Pa",
    "air pressure",
    "air_pressure",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    air_pressure(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    thd,
    "K",
    "dry potential temperature",
    "air_potential_temperature",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    dry_pottemp(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    et,
    "J kg^-1",
    "total specific energy",
    "specific_dry_energy_of_air",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    states.prognostic.ρe
end

@horizontal_average(
    AtmosLESConfigType,
    ei,
    "J kg^-1",
    "specific internal energy",
    "internal_energy",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    internal_energy(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    ht,
    "J kg^-1",
    "specific enthalpy based on total energy",
    "",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    total_specific_enthalpy(ts, states.prognostic.ρe) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    hi,
    "J kg^-1",
    "specific enthalpy based on internal energy",
    "atmosphere_enthalpy_content",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    specific_enthalpy(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    w_ht_sgs,
    "kg kg^-1 m s^-1",
    "vertical sgs flux of total specific enthalpy",
    "",
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    D_t = get!(cache, :D_t) do
        _, D_t, _ = turbulence_tensors(
            atmos,
            states.prognostic,
            states.gradient_flux,
            states.auxiliary,
            curr_time,
        )
        D_t
    end
    d_h_tot = -D_t .* states.gradient_flux.∇h_tot
    d_h_tot[end] * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    eiei,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    internal_energy(ts)^2 * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    wthd,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * dry_pottemp(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    wei,
) do (atmos::AtmosModel, states::States, curr_time, cache)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * internal_energy(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    qt,
    "kg kg^-1",
    "mass fraction of total water in air (qv+ql+qi)",
    "mass_fraction_of_water_in_air",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.moisture.ρq_tot
end

@horizontal_average(
    AtmosLESConfigType,
    ql,
    "kg kg^-1",
    "mass fraction of liquid water in air",
    "mass_fraction_of_cloud_liquid_water_in_air",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    liquid_specific_humidity(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    qi,
    "kg kg^-1",
    "mass fraction of ice in air",
    "mass_fraction_of_cloud_ice_in_air",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    ice_specific_humidity(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    qv,
    "kg kg^-1",
    "mass fraction of water vapor in air",
    "specific_humidity",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    vapor_specific_humidity(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    thv,
    "K",
    "virtual potential temperature",
    "virtual_potential_temperature",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    virtual_pottemp(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    thl,
    "K",
    "liquid-ice potential temperature",
    "",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    liquid_ice_pottemp(ts) * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    w_qt_sgs,
    "kg kg^-1 m s^-1",
    "vertical sgs flux of total specific humidity",
    "",
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    D_t = get!(cache, :D_t) do
        _, D_t, _ = turbulence_tensors(
            atmos,
            states.prognostic,
            states.gradient_flux,
            states.auxiliary,
            curr_time,
        )
        D_t
    end
    d_q_tot = -D_t .* states.gradient_flux.moisture.∇q_tot
    d_q_tot[end] * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    qtqt,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.moisture.ρq_tot^2 / states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    thlthl,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    liquid_ice_pottemp(ts)^2 * states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    wqt,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * states.prognostic.moisture.ρq_tot /
    states.prognostic.ρ
end

@horizontal_average(
    AtmosLESConfigType,
    wql,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * liquid_specific_humidity(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    wqi,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * ice_specific_humidity(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    wqv,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * vapor_specific_humidity(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    wthv,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * virtual_pottemp(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    wthl,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.ρu[3] * liquid_ice_pottemp(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    qtthl,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.moisture.ρq_tot * liquid_ice_pottemp(ts)
end

@horizontal_average(
    AtmosLESConfigType,
    qtei,
) do (
    moisture::Union{EquilMoist, NonEquilMoist},
    atmos::AtmosModel,
    states,
    curr_time,
    cache,
)
    ts = get!(cache, :ts) do
        recover_thermo_state(atmos, states.prognostic, states.auxiliary)
    end
    states.prognostic.moisture.ρq_tot * internal_energy(ts)
end

#= TODO
    Variables["cld_frac"] = DiagnosticVariable(
        "cld_frac",
        diagnostic_var_attrib(
            "",
            "cloud fraction",
            "cloud_area_fraction_in_atmosphere_layer",
        ),
    )
    Variables["cld_cover"] = DiagnosticVariable(
        "cld_cover",
        diagnostic_var_attrib("", "cloud cover", "cloud_area_fraction"),
    )
    Variables["cld_top"] = DiagnosticVariable(
        "cld_top",
        diagnostic_var_attrib("m", "cloud top", "cloud_top_altitude"),
    )
    Variables["cld_base"] = DiagnosticVariable(
        "cld_base",
        diagnostic_var_attrib("m", "cloud base", "cloud_base_altitude"),
    )
    Variables["lwp"] = DiagnosticVariable(
        "lwp",
        diagnostic_var_attrib(
            "kg m^-2",
            "liquid water path",
            "atmosphere_mass_content_of_cloud_condensed_water",
        ),
    )
=#
