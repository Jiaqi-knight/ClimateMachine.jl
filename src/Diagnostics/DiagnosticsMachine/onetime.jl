Base.@kwdef mutable struct CollectedDiagnostics
    onetime_done::Bool = false
    zvals::Union{Nothing, Array} = nothing
    ΣMH_z::Union{Nothing, Array} = nothing
end
const Collected = CollectedDiagnostics()

function collect_onetime(mpicomm, dg, Q)
    if !Collected.onetime_done
        FT = eltype(Q)
        grid = dg.grid
        topology = grid.topology
        N = polynomialorder(grid)
        Nq = N + 1
        Nqk = dimensionality(grid) == 2 ? 1 : Nq
        nrealelem = length(topology.realelems)
        nvertelem = topology.stacksize
        nhorzelem = div(nrealelem, nvertelem)

        vgeo = array_device(Q) isa CPU ? grid.vgeo : Array(grid.vgeo)

        Collected.ΣMH_z = zeros(FT, Nqk, nvertelem)

        for eh in 1:nhorzelem, ev in 1:nvertelem
            e = ev + (eh - 1) * nvertelem
            for k in 1:Nqk, j in 1:Nq, i in 1:Nq
                ijk = i + Nq * ((j - 1) + Nq * (k - 1))
                evk = Nqk * (ev - 1) + k
                MH = vgeo[ijk, grid.MHid, e]
                Collected.ΣMH_z[k, ev] += MH
            end
        end

        # compute the full number of points on a slab on rank 0
        MPI.Reduce!(Collected.ΣMH_z, +, 0, mpicomm)

        Collected.onetime_done = true
    end

    return nothing
end
