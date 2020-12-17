using .FVReconstructions: FVConstant, FVLinear
using ..BalanceLaws: prognostic_to_primitive!, primitive_to_prognostic!

@doc """
    function vert_fvm_interface_tendency!(
        balance_law::BalanceLaw,
        ::Val{info},
        ::Val{nvertelem},
        ::Val{periodicstack},
        ::VerticalDirection,
        numerical_flux_first_order,
        tendency,
        state_prognostic,
        state_auxiliary,
        vgeo,
        sgeo,
        t,
        vmap⁻,
        vmap⁺,
        elemtobndy,
        elems,
        α,
    )

Compute kernel for evaluating the interface tendencies using vertical FVM
reconstructions with a DG method of the form:

∫ₑ ψ⋅ ∂q/∂t dx - ∫ₑ ∇ψ⋅(Fⁱⁿᵛ + Fᵛⁱˢᶜ) dx + ∮ₑ n̂ ψ⋅(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆) dS,

or equivalently in matrix form:

dQ/dt = M⁻¹(MS + DᵀM(Fⁱⁿᵛ + Fᵛⁱˢᶜ) + ∑ᶠ LᵀMf(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆)).

This kernel computes the surface terms: M⁻¹ ∑ᶠ LᵀMf(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆)), where M
is the mass matrix, Mf is the face mass matrix, L is an interpolator from
volume to face, and Fⁱⁿᵛ⋆, Fᵛⁱˢᶜ⋆ are the numerical fluxes for the inviscid
and viscous fluxes, respectively.

A finite volume reconstruction is used to construction `Fⁱⁿᵛ⋆`

!!! note
    Currently `Fᵛⁱˢᶜ⋆` is not supported
""" vert_fvm_interface_tendency!
@kernel function vert_fvm_interface_tendency!(
    balance_law::BalanceLaw,
    ::Val{info},
    ::Val{nvertelem},
    ::Val{periodicstack},
    ::VerticalDirection,
    ::FVConstant,
    numerical_flux_first_order,
    numerical_flux_second_order,
    tendency,
    state_prognostic,
    state_gradient_flux,
    state_auxiliary,
    _,
    sgeo,
    t,
    elemtobndy,
    elems,
    α,
) where {info, nvertelem, periodicstack}
    @uniform begin
        dim = info.dim
        FT = eltype(state_prognostic)
        num_state_prognostic = number_states(balance_law, Prognostic())
        num_state_gradient_flux = number_states(balance_law, GradientFlux())
        num_state_auxiliary = number_states(balance_law, Auxiliary())
        num_state_hyperdiffusion = number_states(balance_law, Hyperdiffusive())
        @assert num_state_hyperdiffusion == 0

        nface = info.nface
        Np = info.Np
        Nqk = info.Nqk # can only be 1 for the FVM method!
        @assert Nqk == 1

        # We only have the vertical faces
        faces = (nface - 1):nface

        local_state_prognostic⁻ = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_gradient_flux⁻ =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)
        local_state_auxiliary⁻ = MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_hyperdiffusion⁻ =
            MArray{Tuple{num_state_hyperdiffusion}, FT}(undef)

        local_state_prognostic⁺ = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_gradient_flux⁺ =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)
        local_state_auxiliary⁺ = MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_hyperdiffusion⁺ =
            MArray{Tuple{num_state_hyperdiffusion}, FT}(undef)

        local_state_prognostic_bottom1 =
            MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_gradient_flux_bottom1 =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)
        local_state_auxiliary_bottom1 =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)

        # XXX: will revisit this later for FVM
        fill!(local_state_prognostic_bottom1, NaN)
        fill!(local_state_gradient_flux_bottom1, NaN)
        fill!(local_state_auxiliary_bottom1, NaN)

        local_flux_top = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_flux_bottom = MArray{Tuple{num_state_prognostic}, FT}(undef)

        sM = MArray{Tuple{2}, FT}(undef)

        # The remainder model needs to know which direction of face the model is
        # being evaluated for. In this case we only have `VerticalDirection()`
        # faces
        face_direction = (EveryDirection(), VerticalDirection())
    end

    # Get the horizontal group IDs
    grp_H = @index(Group, Linear)

    # Determine the index for the element at the bottom of the stack
    eHI = (grp_H - 1) * nvertelem + 1

    # Compute bottom stack element index minus one (so we can add vert element
    # number directly)
    eH = elems[eHI] - 1

    # Which degree of freedom do we handle in the element
    n = @index(Local, Linear)

    # Loads the data for a given element
    function load_data!(
        local_state_prognostic,
        local_state_auxiliary,
        local_state_gradient_flux,
        e,
    )
        @unroll for s in 1:num_state_prognostic
            local_state_prognostic[s] = state_prognostic[n, s, e]
        end

        @unroll for s in 1:num_state_auxiliary
            local_state_auxiliary[s] = state_auxiliary[n, s, e]
        end

        @unroll for s in 1:num_state_gradient_flux
            local_state_gradient_flux[s] = state_gradient_flux[n, s, e]
        end
    end

    # We need to compute the first element we handles bottom flux (future
    # elements will just copied from the prior element)
    @inbounds begin
        eV = 1

        # Minus is the top
        e⁻ = eH + eV

        # bottom face
        f⁻ = faces[1]

        # surface mass
        sM[1] = sgeo[_sM, n, f⁻, e⁻]

        # outward normal for this face
        normal_vector = SVector(
            sgeo[_n1, n, f⁻, e⁻],
            sgeo[_n2, n, f⁻, e⁻],
            sgeo[_n3, n, f⁻, e⁻],
        )

        # determine the plus side element (bottom)
        e⁺, bctag = if periodicstack
            eH + nvertelem, 0
        else
            e⁻, elemtobndy[f⁻, e⁻]
        end

        # Load minus and plus side data
        load_data!(
            local_state_prognostic⁻,
            local_state_auxiliary⁻,
            local_state_gradient_flux⁻,
            e⁻,
        )
        load_data!(
            local_state_prognostic⁺,
            local_state_auxiliary⁺,
            local_state_gradient_flux⁺,
            e⁺,
        )

        # compute the flux
        fill!(local_flux_bottom, -zero(eltype(local_flux_bottom)))
        if bctag == 0
            numerical_flux_first_order!(
                numerical_flux_first_order,
                balance_law,
                local_flux_bottom,
                normal_vector,
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                t,
                face_direction,
            )
            numerical_flux_second_order!(
                numerical_flux_second_order,
                balance_law,
                local_flux_bottom,
                normal_vector,
                local_state_prognostic⁻,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_gradient_flux⁺,
                local_state_hyperdiffusion⁺,
                local_state_auxiliary⁺,
                t,
            )
        else
            numerical_boundary_flux_first_order!(
                numerical_flux_first_order,
                bctag,
                balance_law,
                local_flux_bottom,
                normal_vector,
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                t,
                face_direction,
                local_state_prognostic_bottom1,
                local_state_auxiliary_bottom1,
            )
            numerical_boundary_flux_second_order!(
                numerical_flux_second_order,
                bctag,
                balance_law,
                local_flux_bottom,
                normal_vector,
                local_state_prognostic⁻,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_gradient_flux⁺,
                local_state_hyperdiffusion⁺,
                local_state_auxiliary⁺,
                t,
                local_state_prognostic_bottom1,
                local_state_gradient_flux_bottom1,
                local_state_auxiliary_bottom1,
            )
        end
    end

    # Loop up the vertical stack to update the minus side element (we have the
    # bottom flux from the previous element, so only need to calculate the top
    # flux)
    @inbounds for eV in 1:nvertelem
        e⁻ = eH + eV

        # volume mass inverse
        vMI = sgeo[_vMI, n, faces[1], e⁻]

        # Compute the top face numerical flux
        # The minus side is the bottom element
        # The plus side is the top element
        f⁻ = faces[2]

        # surface mass
        sM[2] = sgeo[_sM, n, f⁻, e⁻]

        # normal with respect to the minus side
        normal_vector = SVector(
            sgeo[_n1, n, f⁻, e⁻],
            sgeo[_n2, n, f⁻, e⁻],
            sgeo[_n3, n, f⁻, e⁻],
        )

        # determine the plus side element (top)
        e⁺, bctag = if eV != nvertelem
            e⁻ + 1, 0
        elseif periodicstack
            eH + 1, 0
        else
            e⁻, elemtobndy[f⁻, e⁻]
        end

        # Load plus side data (minus data is already set)
        load_data!(
            local_state_prognostic⁺,
            local_state_auxiliary⁺,
            local_state_gradient_flux⁺,
            e⁺,
        )

        # compute the flux
        fill!(local_flux_top, -zero(eltype(local_flux_top)))
        if bctag == 0
            numerical_flux_first_order!(
                numerical_flux_first_order,
                balance_law,
                local_flux_top,
                normal_vector,
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                t,
                face_direction,
            )
            numerical_flux_second_order!(
                numerical_flux_second_order,
                balance_law,
                local_flux_top,
                normal_vector,
                local_state_prognostic⁻,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_gradient_flux⁺,
                local_state_hyperdiffusion⁺,
                local_state_auxiliary⁺,
                t,
            )
        else
            numerical_boundary_flux_first_order!(
                numerical_flux_first_order,
                bctag,
                balance_law,
                local_flux_top,
                normal_vector,
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                t,
                face_direction,
                local_state_prognostic_bottom1,
                local_state_auxiliary_bottom1,
            )
            numerical_boundary_flux_second_order!(
                numerical_flux_second_order,
                bctag,
                balance_law,
                local_flux_top,
                normal_vector,
                local_state_prognostic⁻,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                local_state_auxiliary⁻,
                local_state_prognostic⁺,
                local_state_gradient_flux⁺,
                local_state_hyperdiffusion⁺,
                local_state_auxiliary⁺,
                t,
                local_state_prognostic_bottom1,
                local_state_gradient_flux_bottom1,
                local_state_auxiliary_bottom1,
            )
        end

        # Update RHS M⁻¹ Mfᵀ(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆))
        @unroll for s in 1:num_state_prognostic
            # FIXME: Should we pretch these?
            tendency[n, s, e⁻] -=
                α *
                vMI *
                (sM[2] * (local_flux_top[s]) + sM[1] * (local_flux_bottom[s]))
        end

        # Set the flux bottom flux for the next element
        local_flux_bottom .= -local_flux_top

        # set the surface mass matrix
        sM[1] = sM[2]

        # the current plus side is the next minus side
        local_state_prognostic⁻ .= local_state_prognostic⁺
        local_state_auxiliary⁻ .= local_state_auxiliary⁺
        local_state_gradient_flux⁻ .= local_state_gradient_flux⁺
    end
end

@kernel function vert_fvm_interface_tendency!(
    balance_law::BalanceLaw,
    ::Val{info},
    ::Val{nvertelem},
    ::Val{periodicstack},
    ::VerticalDirection,
    reconstruction!::FVLinear,
    numerical_flux_first_order,
    numerical_flux_second_order,
    tendency,
    state_prognostic,
    state_gradient_flux,
    state_auxiliary,
    vgeo,
    sgeo,
    t,
    elemtobndy,
    elems,
    α,
) where {info, nvertelem, periodicstack}
    @uniform begin
        dim = info.dim
        FT = eltype(state_prognostic)
        num_state_prognostic = number_states(balance_law, Prognostic())
        num_state_auxiliary = number_states(balance_law, Auxiliary())
        num_state_gradient_flux = number_states(balance_law, GradientFlux())
        num_state_hyperdiffusion = number_states(balance_law, Hyperdiffusive())
        @assert num_state_hyperdiffusion == 0

        vsp = Vars{vars_state(balance_law, Prognostic(), FT)}
        vsa = Vars{vars_state(balance_law, Auxiliary(), FT)}

        nface = info.nface
        Np = info.Np
        Nqk = info.Nqk # can only be 1 for the FVM method!
        @assert Nqk == 1

        # We only have the vertical faces
        faces = (nface - 1):nface

        # +/- indicate top/bottom elements
        # top and bottom indicate face states

        local_flux = MArray{Tuple{num_state_prognostic}, FT}(undef)

        local_state_prognostic = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_prognostic⁺ = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_prognostic⁻ = MArray{Tuple{num_state_prognostic}, FT}(undef)

        local_state_primitive = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_primitive⁺ = MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_primitive⁻ = MArray{Tuple{num_state_prognostic}, FT}(undef)

        local_state_auxiliary = MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_auxiliary⁺ = MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_auxiliary⁻ = MArray{Tuple{num_state_auxiliary}, FT}(undef)

        local_state_gradient_flux⁻ =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)
        local_state_gradient_flux =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)
        local_state_gradient_flux⁺ =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)

        local_state_hyperdiffusion⁻ =
            MArray{Tuple{num_state_hyperdiffusion}, FT}(undef)
        local_state_hyperdiffusion =
            MArray{Tuple{num_state_hyperdiffusion}, FT}(undef)
        local_state_hyperdiffusion⁺ =
            MArray{Tuple{num_state_hyperdiffusion}, FT}(undef)

        local_state_primitive_bottom =
            MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_primitive_top =
            MArray{Tuple{num_state_prognostic}, FT}(undef)

        local_state_prognostic_bottom =
            MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_prognostic_top =
            MArray{Tuple{num_state_prognostic}, FT}(undef)

        local_state_auxiliary_bottom =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_auxiliary_top =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)

        local_state_prognostic⁻_top =
            MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_prognostic⁺_bottom =
            MArray{Tuple{num_state_prognostic}, FT}(undef)

        local_state_auxiliary⁻_top =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_auxiliary⁺_bottom =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)


        local_state_prognostic_boundary =
            MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_auxiliary_boundary =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)
        local_state_gradient_flux_boundary =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)

        # XXX: will revisit this later for FVM
        fill!(local_state_prognostic_boundary, NaN)
        fill!(local_state_auxiliary_boundary, NaN)
        fill!(local_state_gradient_flux_boundary, NaN)

        # The remainder model needs to know which direction of face the model is
        # being evaluated for. In this case we only have `VerticalDirection()`
        # faces
        face_direction = (VerticalDirection())
    end

    # Get the horizontal group IDs
    grp_H = @index(Group, Linear)

    # Determine the index for the element at the bottom of the stack
    eHI = (grp_H - 1) * nvertelem + 1

    # Compute bottom stack element index minus one (so we can add vert element
    # number directly)
    eH = elems[eHI] - 1

    # Which degree of freedom do we handle in the element
    n = @index(Local, Linear)

    # Loads the data for a given element
    function load_data!(
        local_state_prognostic,
        local_state_auxiliary,
        local_state_gradient_flux,
        e,
    )
        @unroll for s in 1:num_state_prognostic
            local_state_prognostic[s] = state_prognostic[n, s, e]
        end

        @unroll for s in 1:num_state_auxiliary
            local_state_auxiliary[s] = state_auxiliary[n, s, e]
        end

        @unroll for s in 1:num_state_gradient_flux
            local_state_gradient_flux[s] = state_gradient_flux[n, s, e]
        end
    end

    # We need to compute the first element we handles bottom flux (only for nonperiodic boundary condition)
    # elements will just copied from the prior element)
    @inbounds begin
        eV = 1
        # the first element
        e = eH + eV

        # bottom face
        f_bottom = faces[1]
        # surface mass
        sM_bottom = sgeo[_sM, n, f_bottom, e]
        cw = vgeo[n, _M, e]
        vMI = vgeo[n, _MI, e]

        # outward normal for this the bottom face
        normal_vector_bottom = SVector(
            sgeo[_n1, n, f_bottom, e],
            sgeo[_n2, n, f_bottom, e],
            sgeo[_n3, n, f_bottom, e],
        )

        if periodicstack
            e⁻ = eH + nvertelem
            e⁺ = e + 1
            cw⁻ = vgeo[n, _M, e⁻]
            cw⁺ = vgeo[n, _M, e⁺]

            load_data!(
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
                local_state_gradient_flux⁻,
                e⁻,
            )
            load_data!(
                local_state_prognostic,
                local_state_auxiliary,
                local_state_gradient_flux,
                e,
            )
            load_data!(
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                local_state_gradient_flux⁺,
                e⁺,
            )

            prognostic_to_primitive!(
                balance_law,
                local_state_primitive⁻,
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
            )
            prognostic_to_primitive!(
                balance_law,
                local_state_primitive,
                local_state_prognostic,
                local_state_auxiliary,
            )
            prognostic_to_primitive!(
                balance_law,
                local_state_primitive⁺,
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
            )

            cell_states_primitive = (
                local_state_primitive⁻,
                local_state_primitive,
                local_state_primitive⁺,
            )

            cell_weights = SVector(cw⁻, cw, cw⁺)

            # Linear Reconstuction
            reconstruction!(
                local_state_primitive_bottom,
                local_state_primitive_top,
                cell_states_primitive,
                cell_weights,
            )

            # TODO
            local_state_auxiliary_bottom .= local_state_auxiliary
            local_state_auxiliary_top .= local_state_auxiliary

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_bottom,
                local_state_primitive_bottom,
                local_state_auxiliary_bottom,
            )

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_top,
                local_state_primitive_top,
                local_state_auxiliary_top,
            )

            # this is used for the stack top element
            local_state_prognostic⁺_bottom .= local_state_prognostic_bottom
            local_state_auxiliary⁺_bottom .= local_state_auxiliary_bottom
            local_state_prognostic⁻_top .= local_state_prognostic_top
            local_state_auxiliary⁻_top .= local_state_auxiliary_top

        else

            e⁺ = e + 1
            load_data!(
                local_state_prognostic,
                local_state_auxiliary,
                local_state_gradient_flux,
                e,
            )
            load_data!(
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                local_state_gradient_flux⁺,
                e⁺,
            )

            prognostic_to_primitive!(
                balance_law,
                local_state_primitive,
                local_state_prognostic,
                local_state_auxiliary,
            )
            cell_states_primitive = (local_state_primitive,)
            cell_weights = SVector(cw)

            # Constant reconstruction
            reconstruction!(
                local_state_primitive_bottom,
                local_state_primitive_top,
                cell_states_primitive,
                cell_weights,
            )

            # TODO
            local_state_auxiliary_bottom .= local_state_auxiliary
            local_state_auxiliary_top .= local_state_auxiliary

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_bottom,
                local_state_primitive_bottom,
                local_state_auxiliary_bottom,
            )

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_top,
                local_state_primitive_top,
                local_state_auxiliary_top,
            )

            fill!(local_flux, -zero(eltype(local_flux)))
            # TODO
            bctag = elemtobndy[f_bottom, e]
            local_state_prognostic⁻ .= local_state_prognostic_bottom
            local_state_auxiliary⁻ .= local_state_auxiliary_bottom
            numerical_boundary_flux_first_order!(
                numerical_flux_first_order,
                bctag,
                balance_law,
                local_flux,
                normal_vector_bottom,
                local_state_prognostic_bottom,
                local_state_auxiliary_bottom,
                # TODO 
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
                t,
                face_direction,
                local_state_prognostic_boundary,
                local_state_auxiliary_boundary,
            )

            local_state_prognostic⁻ .= local_state_prognostic_bottom
            local_state_gradient_flux⁻ .= local_state_gradient_flux
            local_state_hyperdiffusion⁻ .= local_state_hyperdiffusion
            local_state_auxiliary⁻ .= local_state_auxiliary_bottom
            numerical_boundary_flux_second_order!(
                numerical_flux_second_order,
                bctag,
                balance_law,
                local_flux,
                normal_vector_bottom,
                #
                local_state_prognostic_bottom,
                local_state_gradient_flux,
                local_state_hyperdiffusion,
                local_state_auxiliary_bottom,
                #
                local_state_prognostic⁻,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                local_state_auxiliary⁻,
                t,
                local_state_prognostic_boundary,
                local_state_gradient_flux_boundary,
                local_state_auxiliary_boundary,
            )

            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e] -= α * vMI * sM_bottom * (local_flux[s])
            end

            #TODO
            local_state_prognostic⁻_top .= local_state_prognostic_top
            local_state_auxiliary⁻_top .= local_state_auxiliary_top

        end
    end

    # Loop up the vertical stack to update the minus side element (we have the
    # bottom flux from the previous element, so only need to calculate the top
    # flux)
    @inbounds for eV in 2:(nvertelem - 1)
        e = eH + eV
        e⁻ = e - 1
        e⁺ = e + 1

        # volume mass inverse
        vMI = vgeo[n, _MI, e]
        cw⁻, cw, cw⁺ = vgeo[n, _M, e⁻], vgeo[n, _M, e], vgeo[n, _M, e⁺]

        # Compute the top face numerical flux
        # The minus side is the bottom element
        # The plus side is the top element
        f_bottom = faces[1]
        # surface mass
        sM_bottom = sgeo[_sM, n, f_bottom, e]

        # outward normal with respect to the element
        normal_vector_bottom = SVector(
            sgeo[_n1, n, f_bottom, e],
            sgeo[_n2, n, f_bottom, e],
            sgeo[_n3, n, f_bottom, e],
        )


        # Load plus side data (minus data is already set)
        local_state_prognostic⁻ .= local_state_prognostic
        local_state_auxiliary⁻ .= local_state_auxiliary
        local_state_gradient_flux⁻ .= local_state_gradient_flux
        local_state_prognostic .= local_state_prognostic⁺
        local_state_auxiliary .= local_state_auxiliary⁺
        local_state_gradient_flux .= local_state_gradient_flux⁺
        load_data!(
            local_state_prognostic⁺,
            local_state_auxiliary⁺,
            local_state_gradient_flux⁺,
            e⁺,
        )


        prognostic_to_primitive!(
            balance_law,
            local_state_primitive⁻,
            local_state_prognostic⁻,
            local_state_auxiliary⁻,
        )

        prognostic_to_primitive!(
            balance_law,
            local_state_primitive,
            local_state_prognostic,
            local_state_auxiliary,
        )

        prognostic_to_primitive!(
            balance_law,
            local_state_primitive⁺,
            local_state_prognostic⁺,
            local_state_auxiliary⁺,
        )


        cell_states_primitive = (
            local_state_primitive⁻,
            local_state_primitive,
            local_state_primitive⁺,
        )

        cell_weights = SVector(cw⁻, cw, cw⁺)

        # Linear Reconstuction
        reconstruction!(
            local_state_primitive_bottom,
            local_state_primitive_top,
            cell_states_primitive,
            cell_weights,
        )

        # TODO
        local_state_auxiliary_bottom .= local_state_auxiliary
        local_state_auxiliary_top .= local_state_auxiliary

        primitive_to_prognostic!(
            balance_law,
            local_state_prognostic_bottom,
            local_state_primitive_bottom,
            local_state_auxiliary_bottom,
        )

        primitive_to_prognostic!(
            balance_law,
            local_state_prognostic_top,
            local_state_primitive_top,
            local_state_auxiliary_top,
        )

        ###  TODO HYDROSTATIC BALANCE RECONSTRUCTION

        # compute the flux
        fill!(local_flux, -zero(eltype(local_flux)))

        numerical_flux_first_order!(
            numerical_flux_first_order,
            balance_law,
            local_flux,
            normal_vector_bottom,
            local_state_prognostic_bottom,
            local_state_auxiliary_bottom,
            local_state_prognostic⁻_top,
            local_state_auxiliary⁻_top,
            t,
            face_direction,
        )

        numerical_flux_second_order!(
            numerical_flux_second_order,
            balance_law,
            local_flux,
            normal_vector_bottom,
            #
            # local_state_prognostic_bottom,
            local_state_prognostic,
            local_state_gradient_flux,
            local_state_hyperdiffusion,
            # local_state_auxiliary_bottom,
            local_state_auxiliary,
            #
            # local_state_prognostic⁻_top,
            local_state_prognostic⁻,
            local_state_gradient_flux⁻,
            local_state_hyperdiffusion⁻,
            # local_state_auxiliary⁻_top,
            local_state_auxiliary⁻,
            t,
        )

        # Update RHS (in outer face loop): M⁻¹ Mfᵀ(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆))
        # TODO: This isn't correct:
        # FIXME: Should we pretch these?
        @unroll for s in 1:num_state_prognostic
            tendency[n, s, e⁻] += α * vMI * sM_bottom * (local_flux[s])
        end
        @unroll for s in 1:num_state_prognostic
            tendency[n, s, e] -= α * vMI * sM_bottom * (local_flux[s])
        end

        local_state_prognostic⁻_top .= local_state_prognostic_top
        local_state_auxiliary⁻_top .= local_state_auxiliary_top
    end

    # The top element
    @inbounds begin
        eV = nvertelem
        # the first element
        e = eH + eV
        # bottom face
        f_bottom, f_top = faces[1], faces[2]

        # surface mass
        sM_bottom, sM_top = sgeo[_sM, n, f_bottom, e], sgeo[_sM, n, f_top, e]
        cw = vgeo[n, _M, e]
        vMI = vgeo[n, _MI, e]
        # outward normal for this face
        # outward normal for this the bottom face
        normal_vector_bottom = SVector(
            sgeo[_n1, n, f_bottom, e],
            sgeo[_n2, n, f_bottom, e],
            sgeo[_n3, n, f_bottom, e],
        )

        # outward normal for this the bottom face
        normal_vector_top = SVector(
            sgeo[_n1, n, f_top, e],
            sgeo[_n2, n, f_top, e],
            sgeo[_n3, n, f_top, e],
        )

        if periodicstack

            e⁻ = e - 1
            e⁺ = eH + 1
            cw⁻ = vgeo[n, _M, e⁻]
            cw⁺ = vgeo[n, _M, e⁺]

            local_state_prognostic⁻ .= local_state_prognostic
            local_state_auxiliary⁻ .= local_state_auxiliary
            local_state_gradient_flux⁻ .= local_state_gradient_flux
            local_state_prognostic .= local_state_prognostic⁺
            local_state_auxiliary .= local_state_auxiliary⁺
            local_state_gradient_flux .= local_state_gradient_flux⁺
            load_data!(
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                local_state_gradient_flux⁺,
                e⁺,
            )

            prognostic_to_primitive!(
                balance_law,
                local_state_primitive⁻,
                local_state_prognostic⁻,
                local_state_auxiliary⁻,
            )
            prognostic_to_primitive!(
                balance_law,
                local_state_primitive,
                local_state_prognostic,
                local_state_auxiliary,
            )
            prognostic_to_primitive!(
                balance_law,
                local_state_primitive⁺,
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
            )

            cell_states_primitive = (
                local_state_primitive⁻,
                local_state_primitive,
                local_state_primitive⁺,
            )

            cell_weights = SVector(cw⁻, cw, cw⁺)

            # Linear Reconstuction
            reconstruction!(
                local_state_primitive_bottom,
                local_state_primitive_top,
                cell_states_primitive,
                cell_weights,
            )

            # TODO
            local_state_auxiliary_bottom .= local_state_auxiliary
            local_state_auxiliary_top .= local_state_auxiliary

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_bottom,
                local_state_primitive_bottom,
                local_state_auxiliary_bottom,
            )

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_top,
                local_state_primitive_top,
                local_state_auxiliary_top,
            )

            # compute the bottom flux
            fill!(local_flux, -zero(eltype(local_flux)))

            numerical_flux_first_order!(
                numerical_flux_first_order,
                balance_law,
                local_flux,
                normal_vector_bottom,
                local_state_prognostic_bottom,
                local_state_auxiliary_bottom,
                local_state_prognostic⁻_top,
                local_state_auxiliary⁻_top,
                t,
                face_direction,
            )

            numerical_flux_second_order!(
                numerical_flux_second_order,
                balance_law,
                local_flux,
                normal_vector_bottom,
                #
                # local_state_prognostic_bottom,
                local_state_prognostic,
                local_state_gradient_flux,
                local_state_hyperdiffusion,
                # local_state_auxiliary_bottom,
                local_state_auxiliary,
                #
                # local_state_prognostic⁻_top,
                local_state_prognostic⁻,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                # local_state_auxiliary⁻_top,
                local_state_auxiliary⁻,
                t,
            )

            # Update RHS (in outer face loop): M⁻¹ Mfᵀ(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆))
            # TODO: This isn't correct:
            # FIXME: Should we pretch these?
            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e⁻] += α * vMI * sM_bottom * (local_flux[s])
            end
            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e] -= α * vMI * sM_bottom * (local_flux[s])
            end

            # compute the top flux
            fill!(local_flux, -zero(eltype(local_flux)))

            numerical_flux_first_order!(
                numerical_flux_first_order,
                balance_law,
                local_flux,
                normal_vector_top,
                local_state_prognostic_top,
                local_state_auxiliary_top,
                local_state_prognostic⁺_bottom,
                local_state_auxiliary⁺_bottom,
                t,
                face_direction,
            )

            numerical_flux_second_order!(
                numerical_flux_second_order,
                balance_law,
                local_flux,
                normal_vector_top,
                #
                # local_state_prognostic_top,
                local_state_prognostic,
                local_state_gradient_flux,
                local_state_hyperdiffusion,
                # local_state_auxiliary_top,
                local_state_auxiliary,
                #
                # local_state_prognostic⁺_bottom,
                local_state_prognostic⁺,
                local_state_gradient_flux⁺,
                local_state_hyperdiffusion⁺,
                # local_state_auxiliary⁺_bottom,
                local_state_auxiliary⁺,
                t,
            )

            # Update RHS (in outer face loop): M⁻¹ Mfᵀ(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆))

            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e] -= α * vMI * sM_top * (local_flux[s])
            end
            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e⁺] += α * vMI * sM_top * (local_flux[s])
            end


        else
            e⁻ = e - 1

            # Reset arrays to the top element value
            local_state_gradient_flux⁻ .= local_state_gradient_flux
            local_state_prognostic .= local_state_prognostic⁺
            local_state_auxiliary .= local_state_auxiliary⁺
            local_state_gradient_flux .= local_state_gradient_flux⁺
            local_state_auxiliary_bottom .= local_state_auxiliary
            local_state_auxiliary_top .= local_state_auxiliary

            prognostic_to_primitive!(
                balance_law,
                local_state_primitive,
                local_state_prognostic,
                local_state_auxiliary,
            )
            cell_states_primitive = (local_state_primitive,)
            cell_weights = SVector(cw)

            # Constant reconstruction
            reconstruction!(
                local_state_primitive_bottom,
                local_state_primitive_top,
                cell_states_primitive,
                cell_weights,
            )

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_bottom,
                local_state_primitive_bottom,
                local_state_auxiliary_bottom,
            )

            primitive_to_prognostic!(
                balance_law,
                local_state_prognostic_top,
                local_state_primitive_top,
                local_state_auxiliary_top,
            )

            # bottom flux
            fill!(local_flux, -zero(eltype(local_flux)))

            numerical_flux_first_order!(
                numerical_flux_first_order,
                balance_law,
                local_flux,
                normal_vector_bottom,
                local_state_prognostic_bottom,
                local_state_auxiliary_bottom,
                local_state_prognostic⁻_top,
                local_state_auxiliary⁻_top,
                t,
                face_direction,
            )

            numerical_flux_second_order!(
                numerical_flux_second_order,
                balance_law,
                local_flux,
                normal_vector_bottom,
                #
                local_state_prognostic_bottom,
                local_state_gradient_flux,
                local_state_hyperdiffusion,
                local_state_auxiliary_bottom,
                #
                local_state_prognostic⁻_top,
                local_state_gradient_flux⁻,
                local_state_hyperdiffusion⁻,
                local_state_auxiliary⁻_top,
                t,
            )




            # Update RHS (in outer face loop): M⁻¹ Mfᵀ(Fⁱⁿᵛ⋆ + Fᵛⁱˢᶜ⋆))
            # TODO: This isn't correct:
            # FIXME: Should we pretch these?
            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e] -= α * vMI * sM_bottom * (local_flux[s])
            end
            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e⁻] += α * vMI * sM_bottom * (local_flux[s])
            end


            bctag = elemtobndy[f_top, e]

            fill!(local_flux, -zero(eltype(local_flux)))

            local_state_prognostic⁺ .= local_state_prognostic_top
            local_state_auxiliary⁺ .= local_state_auxiliary_top

            numerical_boundary_flux_first_order!(
                numerical_flux_first_order,
                bctag,
                balance_law,
                local_flux,
                normal_vector_top,
                local_state_prognostic_top,
                local_state_auxiliary_top,
                #TODO
                local_state_prognostic⁺,
                local_state_auxiliary⁺,
                t,
                face_direction,
                local_state_prognostic_boundary,
                local_state_auxiliary_boundary,
            )

            local_state_prognostic⁺ .= local_state_prognostic_top
            local_state_gradient_flux⁺ .= local_state_gradient_flux
            local_state_hyperdiffusion⁺ .= local_state_hyperdiffusion
            local_state_auxiliary⁺ .= local_state_auxiliary_top

            numerical_boundary_flux_second_order!(
                numerical_flux_second_order,
                bctag,
                balance_law,
                local_flux,
                normal_vector_top,
                #
                local_state_prognostic_top,
                local_state_gradient_flux,
                local_state_hyperdiffusion,
                local_state_auxiliary_top,
                #
                local_state_prognostic⁺,
                local_state_gradient_flux⁺,
                local_state_hyperdiffusion⁺,
                local_state_auxiliary⁺,
                #
                t,
                local_state_prognostic_boundary,
                local_state_gradient_flux_boundary,
                local_state_auxiliary_boundary,
            )

            @unroll for s in 1:num_state_prognostic
                tendency[n, s, e] -= α * vMI * sM_top * (local_flux[s])
            end
        end
    end
end

@kernel function vert_fvm_interface_gradients!(
    balance_law::BalanceLaw,
    ::Val{info},
    ::Val{nvertelem},
    ::Val{periodicstack},
    ::VerticalDirection,
    state_prognostic,
    state_gradient_flux,
    state_auxiliary,
    vgeo,
    sgeo,
    t,
    elemtobndy,
    elems,
) where {info, nvertelem, periodicstack}
    @uniform begin

        dim = info.dim
        FT = eltype(state_prognostic)
        num_state_prognostic = number_states(balance_law, Prognostic())
        ngradstate = number_states(balance_law, Gradient())
        num_state_gradient_flux = number_states(balance_law, GradientFlux())
        num_state_auxiliary = number_states(balance_law, Auxiliary())
        nface = info.nface
        Np = info.Np
        faces = (nface - 1, nface)

        # Storage for the prognostic state for e-1, e, e+1
        local_state_prognostic =
            ntuple(_ -> MArray{Tuple{num_state_prognostic}, FT}(undef), Val(3))

        # Storage for the auxiliary state for e-1, e, e+1
        local_state_auxiliary =
            ntuple(_ -> MArray{Tuple{num_state_auxiliary}, FT}(undef), Val(3))

        # Storage for the transform state for e-1, e, e+1
        # e.g., the state we take the gradient of
        l_grad_arg = ntuple(_ -> MArray{Tuple{ngradstate}, FT}(undef), Val(3))

        # Storage for the contribution to the gradient from these faces
        l_nG = MArray{Tuple{3, ngradstate}, FT}(undef)

        l_nG_bc = MArray{Tuple{3, ngradstate}, FT}(undef)

        # Storage for the state gradient flux locally
        local_state_gradient_flux =
            MArray{Tuple{num_state_gradient_flux}, FT}(undef)

        local_state_prognostic_bottom1 =
            MArray{Tuple{num_state_prognostic}, FT}(undef)
        local_state_auxiliary_bottom1 =
            MArray{Tuple{num_state_auxiliary}, FT}(undef)

        # XXX: will revisit this later for FVM
        fill!(local_state_prognostic_bottom1, NaN)
        fill!(local_state_auxiliary_bottom1, NaN)
    end

    # Element index
    eI = @index(Group, Linear)
    # Index of a quadrature point on a face
    n = @index(Local, Linear)

    @inbounds begin
        e = elems[eI]
        eV = mod1(e, nvertelem)

        # Figure out the element above and below e
        e_dn, bc_dn = if eV > 1
            e - 1, 0
        elseif periodicstack
            e + nvertelem - 1, 0
        else
            e, elemtobndy[faces[1], e]
        end

        e_up, bc_up = if eV < nvertelem
            e + 1, 0
        elseif periodicstack
            e - nvertelem + 1, 0
        else
            e, elemtobndy[faces[2], e]
        end

        bctag = (bc_dn, bc_up)

        els = (e_dn, e, e_up)

        # Load the normal vectors and surface mass matrix on the faces
        normal_vector = ntuple(
            k -> SVector(
                sgeo[_n1, n, faces[k], e],
                sgeo[_n2, n, faces[k], e],
                sgeo[_n3, n, faces[k], e],
            ),
            Val(2),
        )
        sM = ntuple(k -> sgeo[_sM, n, faces[k], e], Val(2))

        # volume mass same on both faces
        vMI = sgeo[_vMI, n, faces[1], e]

        # Get the mass matrix for each of the elements
        M = ntuple(k -> vgeo[n, _M, els[k]], Val(3))

        # load prognostic data
        @unroll for k in 1:3
            @unroll for s in 1:num_state_prognostic
                local_state_prognostic[k][s] = state_prognostic[n, s, els[k]]
            end
        end

        # load auxiliary data
        @unroll for k in 1:3
            @unroll for s in 1:num_state_auxiliary
                local_state_auxiliary[k][s] = state_auxiliary[n, s, els[k]]
            end
        end

        # transform to the gradient argument (i.e., the values we take the
        # gradient of)
        @unroll for k in 1:3
            fill!(l_grad_arg[k], -zero(eltype(l_grad_arg[k])))
            compute_gradient_argument!(
                balance_law,
                l_grad_arg[k],
                local_state_prognostic[k],
                local_state_auxiliary[k],
                t,
            )
        end

        # Compute the surface integral contribution from these two faces:
        #   M⁻¹ Lᵀ Mf n̂ G*
        # Since this is FVM we do not subtract the interior state
        fill!(l_nG, -zero(eltype(l_nG)))
        @unroll for f in 1:2
            if bctag[f] == 0
                @unroll for s in 1:ngradstate
                    # Interpolate to the face (using the mass matrix for the
                    # interpolation weights) -- "the gradient numerical flux"
                    G =
                        (
                            M[f] * l_grad_arg[f + 1][s] +
                            M[f + 1] * l_grad_arg[f][s]
                        ) / (M[f] + M[f + 1])

                    # Compute the surface integral for this component and face
                    # multiplied by the normal to get the rotation to the
                    # physical space
                    @unroll for i in 1:3
                        l_nG[i, s] += vMI * sM[f] * normal_vector[f][i] * G
                    end
                end
            else
                # Computes G* incorporating boundary conditions
                numerical_boundary_flux_gradient!(
                    CentralNumericalFluxGradient(),
                    bctag[f],
                    balance_law,
                    l_nG_bc,
                    normal_vector[f],
                    l_grad_arg[2],
                    local_state_prognostic[2],
                    local_state_auxiliary[2],
                    l_grad_arg[2f - 1],
                    local_state_prognostic[2f - 1],
                    local_state_auxiliary[2f - 1],
                    t,
                    local_state_prognostic_bottom1,
                    local_state_auxiliary_bottom1,
                )

                @unroll for s in 1:ngradstate
                    @unroll for i in 1:3
                        l_nG[i, s] += vMI * sM[f] * l_nG_bc[i, s]
                    end
                end
            end
        end

        # Applies linear transformation of gradients to the diffusive variables
        # for storage
        compute_gradient_flux!(
            balance_law,
            local_state_gradient_flux,
            l_nG,
            local_state_prognostic[2],
            local_state_auxiliary[2],
            t,
        )

        # This is the surface integral evaluated discretely
        # M^(-1) Mf G*
        @unroll for s in 1:num_state_gradient_flux
            state_gradient_flux[n, s, e] += local_state_gradient_flux[s]
        end
    end
end
