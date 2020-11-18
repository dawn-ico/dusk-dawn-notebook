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

#include "diff_ops_cxx-naive.cpp"

void dumpCellField(const std::string &fname, const atlas::Mesh &mesh,
                   atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    fprintf(fp, "%.17g\n", field(cellIdx, level));
  }
  fclose(fp);
}

void dumpNodeFieldOnCells(const std::string &fname, const atlas::Mesh &mesh,
                          atlasInterface::Field<double> &field, int level) {
  FILE *fp = fopen(fname.c_str(), "w+");
  const auto &conn = mesh.cells().node_connectivity();
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    double h = 0.;
    for (int nbhIdx = 0; nbhIdx < conn.cols(cellIdx); nbhIdx++) {
      int cIdx = conn(cellIdx, nbhIdx);
      h += std::isfinite(field(cIdx, 0)) ? field(cIdx, 0) : 0.;
    }
    h /= conn.cols(cellIdx);
    fprintf(fp, "%f\n", h);
  }
  fclose(fp);
}

void dumpMesh4Triplot(const atlas::Mesh &mesh, const std::string prefix,
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
}

template <typename T> static int sgn(T val) {
  return (T(0) < val) - (val < T(0));
}

std::tuple<double, double, double>
compare(const atlasInterface::Field<double> &ref,
        const atlasInterface::Field<double> &view,
        const std::vector<bool> &filter, int level) {
  double Linf = 0.;
  double L1 = 0.;
  double L2 = 0.;
  const int N = ref.numElements();
  for (int idx = 0; idx < N; idx++) {
    if (filter[idx]) {
      continue;
    }
    double dif = view(idx, level) - ref(idx, level);
    Linf = fmax(fabs(dif), Linf);
    L1 += fabs(dif);
    L2 += dif * dif;
  }
  L1 /= N;
  L2 = sqrt(L2) / sqrt(N);
  return {L1, L2, Linf};
}

enum class diff_type { gradient, divergence, curl };

int main(int argc, char *argv[]) {
  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);
  if (argc != 2) {
    std::cerr << "Usage: " << argv[0] << "gradient|divergence|curl"
              << std::endl;
    return 1;
  }
  const std::string mode = argv[1];
  if (!(mode == "gradient" || mode == "divergence" || mode == "curl")) {
    std::cerr << "Usage: " << mode << "gradient|divergence|curl" << std::endl;
    return 1;
  }

  diff_type op;
  if (mode == "gradient") {
    op = diff_type::gradient;
  }
  if (mode == "divergence") {
    op = diff_type::divergence;
  }
  if (mode == "curl") {
    op = diff_type::curl;
  }

  // 2 dimensional field
  const int k_size = 1;
  const int level = 0;

  // length of domain
  double lDomain = M_PI;

  // number of triangles per side
  //    (dt needs to be adapted if changed!)
  const int nPerDim = 32;

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
  const int edgesPerNode = 6;
  const int edgesPerCell = 3;

  // system variables
  //  u,v: velocity components
  //  h: height field
  //  h_x, h_y: gradients in x, y direction of height field h
  //  uv_div: divergence of velocity field
  //  suffix E/C -> variabel is located on Edges/Cells
  //  {X}_t: temporal derivative of quantity {X}
  auto [fE_F, fE] = MakeAtlasField("f", mesh.edges().size(), k_size);
  auto [uE_F, uE] = MakeAtlasField("u", mesh.edges().size(), k_size);
  auto [vE_F, vE] = MakeAtlasField("v", mesh.edges().size(), k_size);

  auto [f_xC_F, fC_x] = MakeAtlasField("f_x", mesh.cells().size(), k_size);
  auto [f_yC_F, fC_y] = MakeAtlasField("f_y", mesh.cells().size(), k_size);
  auto [uvC_div_F, uvC_div] =
      MakeAtlasField("f_div", mesh.cells().size(), k_size);
  auto [uvN_curl_F, uvN_curl] =
      MakeAtlasField("f_curl", mesh.nodes().size(), k_size);

  // fields to hold analytical solutions
  auto [f_xC_F_Sol, fC_x_Sol] =
      MakeAtlasField("f_x_Sol", mesh.cells().size(), k_size);
  auto [f_yC_F_Sol, fC_y_Sol] =
      MakeAtlasField("f_y_Sol", mesh.cells().size(), k_size);
  auto [uvC_div_F_Sol, uvC_div_Sol] =
      MakeAtlasField("f_div_Sol", mesh.cells().size(), k_size);
  auto [uvN_curl_F_Sol, uvN_curl_Sol] =
      MakeAtlasField("f_curl_Sol", mesh.nodes().size(), k_size);

  // Geometrical factors on edges
  auto [L_F, L] =
      MakeAtlasField("L", mesh.edges().size(), k_size); // edge length
  auto [dualL_F, dualL] =
      MakeAtlasField("dualL", mesh.edges().size(), k_size); // edge length
  auto [nx_F, nx] =
      MakeAtlasField("nx", mesh.edges().size(), k_size); // normals
  auto [ny_F, ny] = MakeAtlasField("ny", mesh.edges().size(), k_size);

  // Geometrical factors on cells
  auto [A_F, A] =
      MakeAtlasField("A", mesh.cells().size(), k_size); // cell areas
  auto [edge_orientation_cell_F,
        edge_orientation_cell] = // edge orientations to flip normals locally
                                 // for FV ops
      MakeAtlasSparseField("edge_orientation_cell", mesh.cells().size(),
                           edgesPerCell, k_size);

  // Geometrical factors on nodes
  auto [dualA_F, dualA] =
      MakeAtlasField("dualA", mesh.nodes().size(), k_size); // cell areas
  auto [edge_orientation_node_F,
        edge_orientation_node] = // edge orientations to flip normals locally
                                 // for FV ops
      MakeAtlasSparseField("edge_orientation_node", mesh.nodes().size(),
                           edgesPerNode, k_size);

  // Masks for boundary conditions
  std::vector<bool> boundary_edges(mesh.edges().size());
  std::vector<bool> boundary_nodes(mesh.nodes().size());
  std::vector<bool> boundary_cells(mesh.cells().size());

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on edges
  //===------------------------------------------------------------------------------------------===//
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    L(edgeIdx, level) = wrapper.edgeLength(mesh, edgeIdx);
    dualL(edgeIdx, level) = wrapper.dualEdgeLength(mesh, edgeIdx);
    auto [nxe, nye] = wrapper.primalNormal(mesh, edgeIdx);
    nx(edgeIdx, level) = nxe;
    ny(edgeIdx, level) = nye;
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
      Vector primal = {nx(edgeIdx), ny(edgeIdx)};
      edge_orientation_cell(cellIdx, nbhIdx) = sgn(dot(toOutsdie, primal));
    }
    // explanation: the vector cellMidpoint -> edgeMidpoint is guaranteed to
    // point outside. The dot product checks if the edge normal has the same
    // orientation. edgeMidpoint is arbitrary, any point on e would work just
    // as well
  }
  //===------------------------------------------------------------------------------------------===//
  // initialize boundar sets to exclude boundary entitiets from error
  // measurement
  //===------------------------------------------------------------------------------------------===//
  {
    const auto &conn = mesh.edges().cell_connectivity();
    for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      if (conn.cols(edgeIdx) != 2) {
        boundary_edges[edgeIdx] = true;
        continue;
      }
      if (conn(edgeIdx, 0) == conn.missing_value() ||
          conn(edgeIdx, 1) == conn.missing_value()) {
        boundary_edges[edgeIdx] = true;
        continue;
      }
      boundary_edges[edgeIdx] = false;
    }
  }

  {
    const auto &conn = mesh.nodes().edge_connectivity();
    for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
      if (conn.cols(nodeIdx) != 6) {
        boundary_nodes[nodeIdx] = true;
        continue;
      }
      for (int nbh = 0; nbh < 6; nbh++) {
        if (conn(nodeIdx, nbh) == conn.missing_value()) {
          boundary_nodes[nodeIdx] = true;
        }
      }
      if (!boundary_edges[nodeIdx]) {
        boundary_nodes[nodeIdx] = false;
      }
    }
  }

  {
    const auto &conn = mesh.cells().node_connectivity();
    for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
      boundary_cells[cellIdx] = boundary_nodes[conn(cellIdx, 0)] ||
                                boundary_nodes[conn(cellIdx, 1)] ||
                                boundary_nodes[conn(cellIdx, 2)];
    }
  }

  //===------------------------------------------------------------------------------------------===//
  // initialize geometrical info on vertices
  //===------------------------------------------------------------------------------------------===//
  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    dualA(nodeIdx) = wrapper.dualCellArea(mesh, nodeIdx);
  }

  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    const auto &nodeEdgeConnectivity = mesh.nodes().edge_connectivity();
    const auto &edgeNodeConnectivity = mesh.edges().node_connectivity();

    const int missingVal = nodeEdgeConnectivity.missing_value();
    int numNbh = nodeEdgeConnectivity.cols(nodeIdx);

    // arbitrary val at boundary
    bool anyMissing = false;
    for (int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      anyMissing |= nodeEdgeConnectivity(nodeIdx, nbhIdx) == missingVal;
    }
    if (numNbh != 6 || anyMissing) {
      for (int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
        edge_orientation_node(nodeIdx, nbhIdx) = -1;
      }
      continue;
    }

    for (int nbhIdx = 0; nbhIdx < numNbh; nbhIdx++) {
      int edgeIdx = nodeEdgeConnectivity(nodeIdx, nbhIdx);

      int n0 = edgeNodeConnectivity(edgeIdx, 0);
      int n1 = edgeNodeConnectivity(edgeIdx, 1);

      int centerIdx = (n0 == nodeIdx) ? n0 : n1;
      int farIdx = (n0 == nodeIdx) ? n1 : n0;

      auto [xLo, yLo] = wrapper.nodeLocation(centerIdx);
      auto [xHi, yHi] = wrapper.nodeLocation(farIdx);

      Vector edge = {xHi - xLo, yHi - yLo};

      auto [nx, ny] = wrapper.primalNormal(mesh, edgeIdx);
      Vector dualNormal = {ny, -nx};

      double dbg = dot(edge, dualNormal);
      int systemSign = sgn(
          dot(edge, dualNormal)); // geometrical factor "corrects" normal such
                                  // that the resulting system is left handed
      edge_orientation_node(nodeIdx, nbhIdx) = systemSign;
    }
  }

  auto wave = [](double x, double y) -> double { return sin(x) * sin(y); };
  auto waveGrad = [](double x, double y) -> std::tuple<double, double> {
    return {cos(x) * sin(y), sin(x) * cos(y)};
  };

  auto sphericalHarmonic = [](double x,
                              double y) -> std::tuple<double, double> {
    return {0.25 * sqrt(105. / (2 * M_PI)) * cos(2 * x) * cos(y) * cos(y) *
                sin(y),
            0.5 * sqrt(15. / (2 * M_PI)) * cos(x) * cos(y) * sin(y)};
  };
  auto analyticalDivergence = [](double x, double y) {
    return -0.5 * (sqrt(105. / (2 * M_PI))) * sin(2 * x) * cos(y) * cos(y) *
               sin(y) +
           0.5 * sqrt(15. / (2 * M_PI)) * cos(x) *
               (cos(y) * cos(y) - sin(y) * sin(y));
  };
  auto analyticalCurl = [](double x, double y) {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    double dudy =
        c1 * cos(2 * x) * cos(y) * (cos(y) * cos(y) - 2 * sin(y) * sin(y));
    double dvdx = -c2 * cos(y) * sin(x) * sin(y);
    return dvdx - dudy;
  };

  // initialize test functions
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
    fE(edgeIdx, level) = sin(x) * sin(y);
  }

  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [uu, vv] = sphericalHarmonic(x, y);
    uE(edgeIdx, level) = uu;
    vE(edgeIdx, level) = vv;
  }

  // collect analytical solutions
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [x, y] = wrapper.cellMidpoint(mesh, cellIdx);

    auto [fx, fy] = waveGrad(x, y);
    fC_x_Sol(cellIdx, level) = fx;
    fC_y_Sol(cellIdx, level) = fy;

    double uv_div = analyticalDivergence(x, y);
    uvC_div_Sol(cellIdx, level) = uv_div;
  }

  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    auto [x, y] = wrapper.nodeLocation(nodeIdx);
    double uv_curl = analyticalCurl(x, y);
    uvN_curl_Sol(nodeIdx, level) = uv_curl;
  }

  dumpMesh4Triplot(mesh, "out/mesh", wrapper);
  if (op == diff_type::gradient) {
    dawn_generated::cxxnaiveico::gradient<atlasInterface::atlasTag>(
        mesh, k_size, fE, nx, ny, L, A, edge_orientation_cell, fC_x, fC_y)
        .run();
    {
      auto [L1, L2, Linf] = compare(fC_x_Sol, fC_x, boundary_cells, level);
      printf("grad_x L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);
    }
    {
      auto [L1, L2, Linf] = compare(fC_y_Sol, fC_y, boundary_cells, level);
      printf("grad_y L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);
    }
    dumpCellField("out/f_x.txt", mesh, fC_x, level);
    dumpCellField("out/f_y.txt", mesh, fC_y, level);
    dumpCellField("out/f_x_Sol.txt", mesh, fC_x_Sol, level);
    dumpCellField("out/f_y_Sol.txt", mesh, fC_y_Sol, level);
    exit(0);
  }

  if (op == diff_type::divergence) {
    dawn_generated::cxxnaiveico::divergence<atlasInterface::atlasTag>(
        mesh, k_size, uE, vE, nx, ny, L, A, edge_orientation_cell, uvC_div)
        .run();
    auto [L1, L2, Linf] = compare(uvC_div_Sol, uvC_div, boundary_cells, level);
    printf("div_uv L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);
    dumpCellField("out/uv_div.txt", mesh, uvC_div, level);
    dumpCellField("out/uv_div_Sol.txt", mesh, uvC_div_Sol, level);
    exit(0);
  }

  if (op == diff_type::curl) {
    dawn_generated::cxxnaiveico::curl<atlasInterface::atlasTag>(
        mesh, k_size, uE, vE, nx, ny, dualL, dualA, edge_orientation_node,
        uvN_curl)
        .run();
    auto [L1, L2, Linf] =
        compare(uvN_curl_Sol, uvN_curl, boundary_nodes, level);
    for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
        if (boundary_nodes[nodeIdx]) {
            uvN_curl(nodeIdx, level) = 0.;
            uvN_curl_Sol(nodeIdx, level) = 0.;
        }
    }
    printf("curl_uv L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);      
    dumpNodeFieldOnCells("out/uv_curl.txt", mesh, uvN_curl, level);
    dumpNodeFieldOnCells("out/uv_curl_Sol.txt", mesh, uvN_curl_Sol, level);
    exit(0);
  }
}