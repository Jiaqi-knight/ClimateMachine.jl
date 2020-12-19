using MPI
using ClimateMachine
using Logging
using ClimateMachine.Mesh.Topologies
using ClimateMachine.MPIStateArrays

MPI.Initialized() || MPI.Init()
mpicomm = MPI.COMM_WORLD

topl = StackedBrickTopology(
    mpicomm,
    (0:10, 0:10, 0:3);
    periodicity = (false, false, false),
    boundary = ((1, 1), (1, 2), (1, 2)),
)

@show MPI.Comm_rank(mpicomm) length(topl.realelems)

# This file was generated using Literate.jl, https://github.com/fredrikekre/Literate.jl

