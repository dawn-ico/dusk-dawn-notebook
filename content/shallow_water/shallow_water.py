from dusk.script import *


@stencil
def ws_shallow_water(hC: Field[Cell], hC_t: Field[Cell],
                     vC: Field[Cell], vC_t: Field[Cell],
                     uC: Field[Cell], uC_t: Field[Cell],
                     hC_x: Field[Cell], hC_y: Field[Cell], uvC_div: Field[Cell],
                     hE: Field[Edge], vE: Field[Edge], uE: Field[Edge],
                     nx: Field[Edge], ny: Field[Edge], L: Field[Edge], alpha: Field[Edge],
                     boundary_edges: Field[Edge], boundary_cells: Field[Cell],
                     A: Field[Cell], edge_orientation: Field[Cell > Edge],
                     Grav: Field[Cell]):

    with levels_downward:

        # lerp cell quantities to edges
        hE = sum_over(Edge > Cell, hC, weights=[1-alpha, alpha])
        uE = sum_over(Edge > Cell, uC, weights=[1-alpha, alpha])
        vE = sum_over(Edge > Cell, vC, weights=[1-alpha, alpha])
    
        # boundary conditions on cells
        if (boundary_edges):
            uE = 0.
            vE = 0.

        # height field gradient
        hC_x = sum_over(Cell > Edge, hE * nx * L * edge_orientation)
        hC_y = sum_over(Cell > Edge, hE * ny * L * edge_orientation)

        # height field gradient is zero on the boundaries
        if (boundary_cells):
            hC_x = 0.
            hC_y = 0.

        # divergence of velocity field
        uvC_div = sum_over(Cell > Edge, (uE*nx + vE*ny)
                           * edge_orientation * L) / A

        # build ODE's
        uC_t = -Grav * hC_x
        vC_t = -Grav * hC_y
        hC_t = hC * uvC_div        
