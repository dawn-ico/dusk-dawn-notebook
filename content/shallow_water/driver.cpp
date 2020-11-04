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

#include "shallow_water_cxx-naive.cpp"

template <typename T> static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

void dumpCellFieldOnNodes(const std::string &fname, const atlas::Mesh &mesh,
                          AtlasToCartesian wrapper,
                          atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  const auto &conn = mesh.nodes().cell_connectivity();
  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    double h = 0.;
    for (int nbhIdx = 0; nbhIdx < conn.cols(nodeIdx); nbhIdx++) {
      int cIdx = conn(nodeIdx, nbhIdx);
      h += field(cIdx, 0);
    }
    h /= conn.cols(nodeIdx);
    fprintf(fp, "%f\n", h);
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

void dumpNodeField(const std::string &fname, const atlas::Mesh &mesh,
                   atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    fprintf(fp, "%.17g\n", field(nodeIdx, level));
  }
  fclose(fp);
}

void dumpCellField(const std::string &fname, const atlas::Mesh &mesh,
                   atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    fprintf(fp, "%.17g\n", field(cellIdx, level));
  }
  fclose(fp);
}

void dumpEdgeField(const std::string &fname, const atlas::Mesh &mesh,
                   atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    fprintf(fp, "%.17g\n", field(edgeIdx, level));
  }
  fclose(fp);
}

bool compare(const std::vector<double> ref,
             const atlasInterface::Field<double> &view, int level) {
  double Linf = 0.;
  double L1 = 0.;
  double L2 = 0.;
  const int N = ref.size();
  for (int idx = 0; idx < N; idx++) {
    double dif = view(idx, level) - ref[idx];
    Linf = fmax(fabs(dif), Linf);
    L1 += fabs(dif);
    L2 += dif * dif;
  }
  L1 /= N;
  L2 = sqrt(L2) / sqrt(N);
  return L1 < 1e-8 && L2 < 1e-8 && Linf < 1e-8;
}

std::vector<double> readField(const std::string &fname) {
  std::ifstream ifile(fname, std::ios::in);
  double num = 0.0;
  std::vector<double> ret;
  while (ifile >> num) {
    ret.push_back(num);
  }
  return ret;
}

int main(int argc, char *argv[]) {
  // enable floating point exception
  feenableexcept(
      FE_ALL_EXCEPT &
      ~FE_INEXACT); // Enable all floating point exceptions but FE_INEXACT
    
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << "test|run" << std::endl;
    return 1;
  }
  const std::string mode = argv[1];
  if (!(mode == "test" || mode == "run")) {
    std::cerr << "Usage: " << mode << "test|run" << std::endl;
    return 1;
  }

  const bool do_run = (mode == "run");    
    
  // reference level of fluid, make sure to chose this large enough, otherwise
  // initial splash may induce negative fluid height and crash the sim
  const double refHeight = 2.;

  // gravitational constant (change for other planets)
  const double Grav = -9.81;

  // 2 dimensional field
  const int k_size = 1;
  const int level = 0;

  // length of domain
  double lDomain = 10;

  // number of triangles per side
  //    (dt needs to be adapted if changed!)
  const int nPerDim = 20;

  // time stepping
  double t = 0.;
  double dt = 0.002;
  double t_final = 8.;
  int step = 0;

  // how often do we want to dump an output file?
  int outputFreq = 20;

  // how often do you want to make some splashes?
  int slashFreq = 1e3;

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
  //  u,v: velocity components
  //  h: height field
  //  h_x, h_y: gradients in x, y direction of height field h
  //  uv_div: divergence of velocity field
  //  suffix E/C -> variabel is located on Edges/Cells
  //  {X}_t: temporal derivative of quantity {X}
  auto [uC_F, uC] = MakeAtlasField("u", mesh.cells().size(), k_size);
  auto [vC_F, vC] = MakeAtlasField("v", mesh.cells().size(), k_size);
  auto [hC_F, hC] = MakeAtlasField("h", mesh.cells().size(), k_size);

  auto [uE_F, uE] = MakeAtlasField("u", mesh.edges().size(), k_size);
  auto [vE_F, vE] = MakeAtlasField("v", mesh.edges().size(), k_size);
  auto [hE_F, hE] = MakeAtlasField("h", mesh.edges().size(), k_size);

  auto [h_xC_F, hC_x] = MakeAtlasField("h_x", mesh.cells().size(), k_size);
  auto [h_yC_F, hC_y] = MakeAtlasField("h_y", mesh.cells().size(), k_size);
  auto [uvC_div_F, uvC_div] =
      MakeAtlasField("uvC_div", mesh.cells().size(), k_size);

  auto [u_tC_F, uC_t] = MakeAtlasField("u_t", mesh.cells().size(), k_size);
  auto [v_tC_F, vC_t] = MakeAtlasField("v_t", mesh.cells().size(), k_size);
  auto [h_tC_F, hC_t] = MakeAtlasField("h_t", mesh.cells().size(), k_size);

  auto [uinitC_F, uinitC] =
      MakeAtlasField("uinit", mesh.cells().size(), k_size);
  auto [vinitC_F, vinitC] =
      MakeAtlasField("vinit", mesh.cells().size(), k_size);
  auto [hinitC_F, hinitC] =
      MakeAtlasField("hinit", mesh.cells().size(), k_size);

  // Geometrical factors on edges
  auto [L_F, L] =
      MakeAtlasField("L", mesh.edges().size(), k_size); // edge length
  auto [nx_F, nx] =
      MakeAtlasField("nx", mesh.edges().size(), k_size); // normals
  auto [ny_F, ny] = MakeAtlasField("ny", mesh.edges().size(), k_size);
  auto [alpha_F, alpha] =
      MakeAtlasField("alpha", mesh.edges().size(),
                     k_size); // linear interpolation coefficients to cells

  // Geometrical factors on cells
  auto [A_F, A] =
      MakeAtlasField("A", mesh.cells().size(), k_size); // cell areas
  auto [edge_orientation_cell_F,
        edge_orientation_cell] = // edge orientations to flip normals locally
                                 // for FV ops
      MakeAtlasSparseField("edge_orientation_cell", mesh.cells().size(),
                           edgesPerCell, k_size);

  // Masks for boundary conditions
  auto [boundary_edges_F, boundary_edges] =
      MakeAtlasField("boundary_edges", mesh.edges().size(), k_size);
  auto [boundary_cells_F, boundary_cells] =
      MakeAtlasField("boundary_cells", mesh.cells().size(), k_size);

  // These should be globals (currently emulated with fields)
  auto [Grav_Field_F, Grav_Field] =
      MakeAtlasField("Grav_Field", mesh.cells().size(), k_size);

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    L(edgeIdx, level) = wrapper.edgeLength(mesh, edgeIdx);
    auto [nxe, nye] = wrapper.primalNormal(mesh, edgeIdx);
    nx(edgeIdx, level) = nxe;
    ny(edgeIdx, level) = nye;
  }

  {
    const auto &conn = mesh.edges().cell_connectivity();
    for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      int cellIdx1 = conn(edgeIdx, 0);
      int cellIdx2 = conn(edgeIdx, 1);
      double d1 = (cellIdx1 >= 0)
                      ? wrapper.distanceToCircumcenter(mesh, cellIdx1, edgeIdx)
                      : 0.;
      double d2 = (cellIdx2 >= 0)
                      ? wrapper.distanceToCircumcenter(mesh, cellIdx2, edgeIdx)
                      : 0.;
      alpha(edgeIdx, level) = d2 / (d1 + d2);
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on cells
  //===------------------------------------------------------------------------------------------===//
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    A(cellIdx, level) = wrapper.cellArea(mesh, cellIdx);
  }

  auto dot = [](const Vector &v1, const Vector &v2) {
    return std::get<0>(v1) * std::get<0>(v2) +
           std::get<1>(v1) * std::get<1>(v2);
  };
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    const atlas::mesh::HybridElements::Connectivity &cellEdgeConnectivity =
        mesh.cells().edge_connectivity();
    auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);

    const int missingVal = cellEdgeConnectivity.missing_value();
    int numNbh = cellEdgeConnectivity.cols(cellIdx);
    assert(numNbh == edgesPerCell);

    for (int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = cellEdgeConnectivity(cellIdx, nbhIdx);
      auto [emX, emY] = wrapper.edgeMidpoint(mesh, edgeIdx);
      Vector toOutsdie{emX - xm, emY - ym};
      Vector primal = {nx(edgeIdx, level), ny(edgeIdx, level)};
      edge_orientation_cell(cellIdx, nbhIdx, level) =
          sgn(dot(toOutsdie, primal));
    }
    // explanation: the vector cellMidpoint -> edgeMidpoint is guaranteed to
    // point outside. The dot product checks if the edge normal has the same
    // orientation. edgeMidpoint is arbitrary, any point on e would work just as
    // well
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize "global" fields
  //===------------------------------------------------------------------------------------------===//
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    Grav_Field(cellIdx, level) = Grav;
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize height and velocity fields
  //===------------------------------------------------------------------------------------------===//
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
    xm -= 0;
    ym -= 0;
    double v = sqrt(xm * xm + ym * ym);
    hC(cellIdx, level) = exp(-5 * v * v) + refHeight;
  }

  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    uC(cellIdx, level) = 0.;
    vC(cellIdx, level) = 0.;

    hC_t(cellIdx, level) = 0.;
    uC_t(cellIdx, level) = 0.;
    vC_t(cellIdx, level) = 0.;
  }

  dumpMesh4Triplot(mesh, "init", hC, wrapper);

  //===------------------------------------------------------------------------------------------===//
  // set up boundary masks
  //===------------------------------------------------------------------------------------------===//
  std::set<int> boundaryEdgesSet;
  {
    const auto &conn = mesh.edges().cell_connectivity();
    for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      if (conn.cols(edgeIdx) < 2) {
        boundary_edges(edgeIdx, level) = 1.;
        boundaryEdgesSet.insert(edgeIdx);
        continue;
      }

      int cellIdx0 = conn(edgeIdx, 0);
      int cellIdx1 = conn(edgeIdx, 1);
      if (cellIdx0 == conn.missing_value() ||
          cellIdx1 == conn.missing_value()) {
        boundary_edges(edgeIdx, level) = 1.;
        boundaryEdgesSet.insert(edgeIdx);
        continue;
      }

      boundary_edges(edgeIdx, level) = 0.;
    }
  }
  {
    const auto &conn = mesh.cells().edge_connectivity();
    for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      int edgeIdx0 = conn(cellIdx, 0);
      int edgeIdx1 = conn(cellIdx, 1);
      int edgeIdx2 = conn(cellIdx, 2);

      boundary_cells(cellIdx, level) =
          (boundaryEdgesSet.count(edgeIdx0) ||
           boundaryEdgesSet.count(edgeIdx1) || boundaryEdgesSet.count(edgeIdx2))
              ? 1.
              : 0.;
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // enter time loop
  //===------------------------------------------------------------------------------------------===//

  // auto globalMesh = atlasToGlobalGpuTriMesh(mesh);
  if (!do_run) {
    // save old result
    for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      hinitC(cellIdx, level) = hC(cellIdx, level);
      vinitC(cellIdx, level) = vC(cellIdx, level);
      uinitC(cellIdx, level) = uC(cellIdx, level);
    }

    // predict
    for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      hC(cellIdx, level) =
          hinitC(cellIdx, level) + 0.5 * dt * hC_t(cellIdx, level);
      uC(cellIdx, level) =
          uinitC(cellIdx, level) + 0.5 * dt * uC_t(cellIdx, level);
      vC(cellIdx, level) =
          vinitC(cellIdx, level) + 0.5 * dt * vC_t(cellIdx, level);
    }

    // run shallow water solver for a single timestep
    dawn_generated::cxxnaiveico::shallow_water<atlasInterface::atlasTag>(
        mesh, k_size, hC, hC_t, vC, vC_t, uC, uC_t, hC_x, hC_y, uvC_div, hE, vE,
        uE, nx, ny, L, alpha, boundary_edges, boundary_cells, A,
        edge_orientation_cell, Grav_Field)
        .run();

    auto checkField =
        [&level](const std::string &fname,
                 const atlasInterface::Field<double> &view) -> bool {
      auto field = readField("ref/" + fname + "_ref.txt");
      if (!compare(field, view, level)) {
        std::cout << "looks like there is a mistake in " << fname << "\n";
        return false;
      } else {
        return true;
      }
    };

    bool allCorrect = true;
    allCorrect &= checkField("hE", hE);
    allCorrect &= checkField("uE", uE);
    allCorrect &= checkField("vE", vE);
    allCorrect &= checkField("hC_x", hC_x);
    allCorrect &= checkField("hC_y", hC_y);
    allCorrect &= checkField("uvC_div", uvC_div);
    allCorrect &= checkField("hC_t", hC_t);
    allCorrect &= checkField("uC_t", uC_t);
    allCorrect &= checkField("vC_t", vC_t);

    if (allCorrect) {
      std::cout << "congratulations!, your stencil is correct, you can "
                   "visualize now!\n";
    }

  } else {
    while (t < t_final) {
      // make some splashes
      if (step > 0 && step % slashFreq == 0) {
        double splash_x = (rand() / ((double)RAND_MAX) - 0.5) * 2;
        double splash_y = (rand() / ((double)RAND_MAX) - 0.5) * 2;
        printf("splashing at %f %f\n", splash_x, splash_y);
        for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
          auto [xm, ym] = wrapper.cellCircumcenter(mesh, cellIdx);
          xm += splash_x;
          ym += splash_y;
          double v = sqrt(xm * xm + ym * ym);
          hC(cellIdx, level) += exp(-5 * v * v);
        }
      }

      // increment counters
      t += dt;
      step++;

      // save old result
      for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        hinitC(cellIdx, level) = hC(cellIdx, level);
        vinitC(cellIdx, level) = vC(cellIdx, level);
        uinitC(cellIdx, level) = uC(cellIdx, level);
      }

      // predict
      for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        hC(cellIdx, level) =
            hinitC(cellIdx, level) + 0.5 * dt * hC_t(cellIdx, level);
        uC(cellIdx, level) =
            uinitC(cellIdx, level) + 0.5 * dt * uC_t(cellIdx, level);
        vC(cellIdx, level) =
            vinitC(cellIdx, level) + 0.5 * dt * vC_t(cellIdx, level);
      }

      // run shallow water solver for a single timestep
      dawn_generated::cxxnaiveico::shallow_water<atlasInterface::atlasTag>(
          mesh, k_size, hC, hC_t, vC, vC_t, uC, uC_t, hC_x, hC_y, uvC_div, hE,
          vE, uE, nx, ny, L, alpha, boundary_edges, boundary_cells, A,
          edge_orientation_cell, Grav_Field)
          .run();

      // correct
      for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
        hC(cellIdx, level) = hinitC(cellIdx, level) + dt * hC_t(cellIdx, level);
        uC(cellIdx, level) = uinitC(cellIdx, level) + dt * uC_t(cellIdx, level);
        vC(cellIdx, level) = vinitC(cellIdx, level) + dt * vC_t(cellIdx, level);
      }

      // dump output
      if (step % outputFreq == 0) {
        char buf[256];
        sprintf(buf, "out/stepH_%04d.txt", step);
        dumpCellFieldOnNodes(buf, mesh, wrapper, hC, level);
      }

      // diagnostics
      std::cout << "time " << t << " timestep " << step << " dt " << dt << "\n";
    }
  }
}

// dumpEdgeField("hE_ref.txt", mesh, hE, level);
// dumpEdgeField("uE_ref.txt", mesh, uE, level);
// dumpEdgeField("vE_ref.txt", mesh, vE, level);
// dumpCellField("hC_x_ref.txt", mesh, hC_x, level);
// dumpCellField("hC_y_ref.txt", mesh, hC_y, level);
// dumpCellField("uvC_div_ref.txt", mesh, uvC_div, level);
// dumpCellField("hC_t_ref.txt", mesh, hC_t, level);
// dumpCellField("uC_t_ref.txt", mesh, uC_t, level);
// dumpCellField("vC_t_ref.txt", mesh, vC_t, level);