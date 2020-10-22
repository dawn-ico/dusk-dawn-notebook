from dusk.script import *


@stencil
def diffusion(
    TV: Field[Vertex],
    TE: Field[Edge], TEinit: Field[Edge], TE_t: Field[Edge], TEnabla2: Field[Edge],
    inv_primal_edge_length: Field[Edge], inv_vert_vert_length: Field[Edge], nnbhV: Field[Vertex],
    boundary_edge: Field[Edge],
    kappa: Field[Edge], dt: Field[Edge]
) -> None:

    with levels_upward:
        # initialize
        TEinit = TE
    with levels_upward:
        # predict
        TE = TEinit + 0.5*dt*TE_t

        # interpolate temperature from edges to vertices
        TV = sum_over(Vertex > Edge, TE) / nnbhV

        # compute nabla2 using the finite differences
        TEnabla2 = sum_over(
            Edge > Cell > Vertex,
            4.0 * TV,
            weights=[
                inv_primal_edge_length ** 2.,
                inv_primal_edge_length ** 2.,
                inv_vert_vert_length ** 2.,
                inv_vert_vert_length ** 2.,
            ],
        )
        TEnabla2 = TEnabla2 - (
            (8.0 * TE * inv_primal_edge_length ** 2.)
            + (8.0 * TE * inv_vert_vert_length ** 2.)
        )

    with levels_upward:
        # build ODEs
        if (boundary_edge):
            TE_t = 0.
        else:
            TE_t = kappa*TEnabla2

    with levels_upward:
        # correct
        TE = TEinit + dt*TE_t
