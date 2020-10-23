// atlas functions
#include <atlas/array.h>
#include <atlas/grid.h>
#include <atlas/mesh.h>
#include <atlas/mesh/actions/BuildEdges.h>
#include <atlas/output/Gmsh.h>
#include <atlas/util/CoordinateEnums.h>

// atlas interface for dawn generated code
#include "interface/atlas_interface.hpp"

// driver includes
#include "driver-includes/defs.hpp"

// atlas utilities
#include "atlas_utils/utils/AtlasCartesianWrapper.h"
#include "atlas_utils/utils/GenerateRectAtlasMesh.h"

#include <cmath>
#include <cstdio>
#include <fenv.h>
#include <optional>
#include <vector>

#include "atlas_utils/utils/getAtlasFields.h"

#include <driver-includes/math.hpp>

#include "diffusion_cxx-naive.cpp"

void dumpEdgeFieldOnNodes(const std::string &fname, const atlas::Mesh &mesh,
                          AtlasToCartesian wrapper,
                          atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  const auto &conn = mesh.nodes().edge_connectivity();
  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double h = 0.;
    for (int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
      int eIdx = conn(nodeIdx, nbhIdx);
      h += field(eIdx, 0);
    }
    h /= conn.cols(nodeIdx);
    fprintf(fp, "%f\n", h);
  }
  fclose(fp);
}

void dumpNodeField(const std::string &fname, const atlas::Mesh &mesh,
                   atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    fprintf(fp, "%f\n", field(nodeIdx, level));
  }
  fclose(fp);
}

void dumpMesh4Triplot(const atlas::Mesh &mesh, const std::string prefix,
                      const atlasInterface::Field<double> &field,
                      std::optional<AtlasToCartesian> wrapper) {
  auto xy = atlas::array::make_view<double, 2>(mesh.nodes().xy());
  const atlas::mesh::HybridElements::Connectivity &node_connectivity =
      mesh.cells().node_connectivity();

  {
    char buf[256];
    sprintf(buf, "%sT.txt", prefix.c_str());
    FILE *fp = fopen(buf, "w+");
    for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      int nodeIdx0 = node_connectivity(cellIdx, 0) + 1;
      int nodeIdx1 = node_connectivity(cellIdx, 1) + 1;
      int nodeIdx2 = node_connectivity(cellIdx, 2) + 1;
      fprintf(fp, "%d %d %d\n", nodeIdx0, nodeIdx1, nodeIdx2);
    }
    fclose(fp);
  }

  {
    char buf[256];
    sprintf(buf, "%sP.txt", prefix.c_str());
    FILE *fp = fopen(buf, "w+");
    for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      if (wrapper == std::nullopt) {
        double x = xy(nodeIdx, atlas::LON);
        double y = xy(nodeIdx, atlas::LAT);
        fprintf(fp, "%f %f \n", x, y);
      } else {
        auto [x, y] = wrapper.value().nodeLocation(nodeIdx);
        fprintf(fp, "%f %f \n", x, y);
      }
    }
    fclose(fp);
  }

  {
    char buf[256];
    sprintf(buf, "%sC.txt", prefix.c_str());
    FILE *fp = fopen(buf, "w+");
    for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      fprintf(fp, "%f\n", field(cellIdx, 0));
    }
    fclose(fp);
  }
}

template <typename T> static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

int main() {
  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  // thermal diffusitivity (change for different materials)
  const double kappa = 1;

  // 2 dimensional field
  const int k_size = 1;
  const int level = 0;

  // length of domain
  double lDomain = 10;

  // number of triangles per side
  //    (dt needs to be adapted if changed!)
  const int nPerDim = 20;

  // time stepping
  double t = 0.1; // inital time (needs to be larger than zero because
                  // otherwise solution deteriorates to dirac pulse)
  const double dx = lDomain / nPerDim;
  const double CFL_const = 0.02;
  const double dt = CFL_const * dx * dx / kappa;
  const double t_final = 1.;
  int step = 0;

  // how often do we want to dump an output file?
  int outputFreq = 5;

  // generating a mesh
  auto mesh = AtlasMeshRect(nPerDim * sqrt(3), nPerDim);

  // todo: hide or move this stuff
  atlas::mesh::actions::build_edges(mesh,
                                    atlas::util::Config("pole_edges", false));
  atlas::mesh::actions::build_node_to_edge_connectivity(mesh);
  atlas::mesh::actions::build_element_to_edge_connectivity(mesh);

  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const auto &nodeToEdge = mesh.nodes().edge_connectivity();
    const auto &edgeToCell = mesh.edges().cell_connectivity();
    auto &nodeToCell = mesh.nodes().cell_connectivity();

    std::set<int> nbh;
    for (int nbhEdgeIdx = 0; nbhEdgeIdx < nodeToEdge.cols(nodeIdx);
         nbhEdgeIdx++) {
      int edgeIdx = nodeToEdge(nodeIdx, nbhEdgeIdx);
      if (edgeIdx == nodeToEdge.missing_value()) {
        continue;
      }
      for (int nbhCellIdx = 0; nbhCellIdx < edgeToCell.cols(edgeIdx);
           nbhCellIdx++) {
        int cellIdx = edgeToCell(edgeIdx, nbhCellIdx);
        if (cellIdx == edgeToCell.missing_value()) {
          continue;
        }
        nbh.insert(cellIdx);
      }
    }

    assert(nbh.size() <= 6);
    std::vector<int> initData(nbh.size(), nodeToCell.missing_value());
    nodeToCell.add(1, nbh.size(), initData.data());
    int copyIter = 0;
    for (const int n : nbh) {
      nodeToCell.set(nodeIdx, copyIter++, n);
    }
  }

  printf("mesh stats: #cells %d #nodes %d #edges %d\n", mesh->cells().size(),
         mesh->nodes().size(), mesh->edges().size());

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh, lDomain, false, true);

  // mesh properties
  const int edgesPerVertex = 6;
  const int edgesPerCell = 3;

  // system variables
  //  T: temperature
  //  Tnabla2: laplacian of temperature field
  //  suffix E/C/V -> variabel is located on Edges/Cells/Vertices
  //  {X}_t: temporal derivative of quantity {X}
  auto [TN_F, TN] = MakeAtlasField("TN", mesh.nodes().size(), k_size);
  auto [TNSol_F, TNSol] = MakeAtlasField("TNSol", mesh.nodes().size(),
                                         k_size); // analytical solution
  auto [TC_F, TC] = MakeAtlasField("TC", mesh.cells().size(), k_size);

  auto [TE_F, TE] = MakeAtlasField("TE", mesh.edges().size(), k_size);
  auto [TEinit_F, TEinit] =
      MakeAtlasField("TEinit", mesh.edges().size(), k_size);
  auto [TE_t_F, TE_t] = MakeAtlasField("TEinit", mesh.edges().size(), k_size);
  auto [TEnabla2_F, TEnabla2] =
      MakeAtlasField("TEnabla2", mesh.edges().size(), k_size);

  // Geometrical factors on edges
  auto [inv_L_F, inv_L] = MakeAtlasField("invL", mesh.edges().size(),
                                         k_size); // inverted edge length
  auto [inv_vvL_F, inv_vvL] = MakeAtlasField(
      "invvvL", mesh.edges().size(), k_size); // inverted vert vert length

  // Geometrical factors on vertices
  auto [nnbhV_F, nnbhV] = MakeAtlasField("nnbh", mesh.nodes().size(),
                                         k_size); // number of edge neighbors
                                                  // for simple averaging

  // Masks for boundary conditions
  auto [boundary_edges_F, boundary_edges] =
      MakeAtlasField("boundary_edges", mesh.edges().size(), k_size);
  auto [boundary_cells_F, boundary_cells] =
      MakeAtlasField("boundary_cells", mesh.edges().size(), k_size);

  // These should be globals (currently emulated with fields)
  auto [kappa_Field_F, kappa_Field] =
      MakeAtlasField("kappa_Field", mesh.edges().size(), k_size);
  auto [dt_Field_F, dt_Field] =
      MakeAtlasField("dt_Field", mesh.edges().size(), k_size);

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    inv_L(edgeIdx, level) = 1. / wrapper.edgeLength(mesh, edgeIdx);
    double vert_vert_length = sqrt(3.) * wrapper.edgeLength(mesh, edgeIdx);
    inv_vvL(edgeIdx, level) =
        (vert_vert_length == 0) ? 0 : 1. / vert_vert_length;
  }

  {
    const auto &conn = mesh.nodes().edge_connectivity();
    for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      nnbhV(nodeIdx, level) = conn.cols(nodeIdx);
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize "global" fields
  //===------------------------------------------------------------------------------------------===//
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    dt_Field(edgeIdx, level) = dt;
    kappa_Field(edgeIdx, level) = kappa;
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize temperature fields
  //===------------------------------------------------------------------------------------------===//
  auto solution = [kappa](double r, double t) {
    return exp(-r / (4 * kappa * t)) / (4 * M_PI * kappa * t);
  };

  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [xm, ym] = wrapper.edgeMidpoint(mesh, edgeIdx);
    double r = sqrt(xm * xm + ym * ym);
    TE(edgeIdx, level) = solution(r, t);
    TE_t(edgeIdx, level) = 0.;
  }
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    const auto &conn = mesh.cells().edge_connectivity();
    double T = 0.;
    for (int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
      int edgeIdx = conn(cellIdx, nbhIdx);
      if (edgeIdx == conn.missing_value()) {
        continue;
      }
      T += TE(edgeIdx, level);
    }
    T /= conn.cols(cellIdx);
    TC(cellIdx, level) = T;
  }

  //===------------------------------------------------------------------------------------------===//
  // boundary edges
  //===------------------------------------------------------------------------------------------===//
  {
    const auto &conn = mesh.edges().cell_connectivity();
    for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      if (conn.cols(edgeIdx) < 2) {
        boundary_edges(edgeIdx, level) = 1.;
        continue;
      }

      int cellIdx0 = conn(edgeIdx, 0);
      int cellIdx1 = conn(edgeIdx, 1);
      if (cellIdx0 == conn.missing_value() ||
          cellIdx1 == conn.missing_value()) {
        boundary_edges(edgeIdx, level) = 1.;
        continue;
      }

      boundary_edges(edgeIdx, level) = 0.;
    }
  }

  dumpMesh4Triplot(mesh, "init", TC, wrapper);

  //===------------------------------------------------------------------------------------------===//
  // enter time loop
  //===------------------------------------------------------------------------------------------===//

  while (t < t_final) {
    // increment counters
    t += dt;
    step++;

    // run shallow water solver for a single timestep
    dawn_generated::cxxnaiveico::diffusion<atlasInterface::atlasTag>(
        mesh, k_size, TN, TE, TEinit, TE_t, TEnabla2, inv_L, inv_vvL, nnbhV,
        boundary_edges, kappa_Field, dt_Field)
        .run();

    // dump output
    if (step % outputFreq == 0) {
      char buf[256];
      sprintf(buf, "out/stepT_%04d.txt", step);
      dumpEdgeFieldOnNodes(buf, mesh, wrapper, TE, level);
    }

    // dump analytical solution
    if (step % outputFreq == 0) {
      char buf[256];
      sprintf(buf, "out/solT_%04d.txt", step);
      for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
        auto [x, y] = wrapper.nodeLocation(nodeIdx);
        double r = sqrt(x * x + y * y);
        TNSol(nodeIdx, level) = solution(r, t);
      }
      dumpNodeField(buf, mesh, TNSol, level);
    }

    // diagnostics
    std::cout << "time " << t << " timestep " << step << " dt " << dt << "\n";
  }
}
