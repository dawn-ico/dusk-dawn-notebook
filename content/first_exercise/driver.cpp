// WARNING: this driver is very ugly and full of hacks!
// You'll probably have better luck with other drivers...

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
#include <iostream>
#include <fstream>
#include <string>

#include "atlas_utils/utils/getAtlasFields.h"

#include <driver-includes/math.hpp>

#if STENCIL_TYPE == 0
#include "simple_stencil_cxx-naive.cpp"
#elif STENCIL_TYPE == 1
#include "image_stencil_cxx-naive.cpp"
#else
static_assert(false);
#endif

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

int main(int argc, char *argv[]) {
  // enable floating point exception
  feenableexcept(FE_INVALID | FE_OVERFLOW);

  // 2 dimensional field
  const int k_size = 1;
  const int level = 0;

#if STENCIL_TYPE == 0
  int w = 54, h = 32;
#elif STENCIL_TYPE == 1
  std::string filename = "./picture.txt";
  std::ifstream picture(filename);
  if (!picture.is_open()) {
    std::cout << "Failed to open '" << filename << "'!" << std::endl;
    exit(-1);
  }
  int w, h;
  picture >> h >> w;
#endif

  // generating a mesh
  auto mesh = AtlasMeshRect(w, h);

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

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh, w, false, false);

  // mesh properties
  const int edgesPerVertex = 6;
  const int edgesPerNode = 6;
  const int edgesPerCell = 3;

  // system variables
#if STENCIL_TYPE == 0
  auto [a, a_view] = MakeAtlasField("a", mesh.cells().size(), k_size);
  auto [b, b_view] = MakeAtlasField("b", mesh.cells().size(), k_size);
  auto [c, c_view] = MakeAtlasField("c", mesh.cells().size(), k_size);
  const auto N = a_view.numElements();

  for (int idx = 0; idx < N; idx++) {

    c_view(idx, level) = -1;

    const auto [x, y] = wrapper.cellMidpoint(mesh, idx);
    a_view(idx, level) = x;
    b_view(idx, level) = y;
  }

  dumpCellField("out/a.txt", mesh, a_view, level);
  dumpCellField("out/b.txt", mesh, b_view, level);

  dawn_generated::cxxnaiveico::simple_stencil<atlasInterface::atlasTag>(
      mesh,
      k_size,
      a_view,
      b_view,
      c_view
    )
      .run();

  dumpMesh4Triplot(mesh, "out/mesh", wrapper);
  dumpCellField("out/c.txt", mesh, c_view, level);

  printf("Finished simple stencil successfully.\n");
#elif STENCIL_TYPE == 1
  auto [r, r_view] = MakeAtlasField("r", mesh.cells().size(), k_size);
  auto [g, g_view] = MakeAtlasField("g", mesh.cells().size(), k_size);
  auto [b, b_view] = MakeAtlasField("b", mesh.cells().size(), k_size);
  auto [x, x_view] = MakeAtlasField("x", mesh.cells().size(), k_size);
  auto [y, y_view] = MakeAtlasField("y", mesh.cells().size(), k_size);
  auto [val1, val1_view] = MakeAtlasField("val1", mesh.cells().size(), k_size);
  auto [val2, val2_view] = MakeAtlasField("val2", mesh.cells().size(), k_size);
  auto [val3, val3_view] = MakeAtlasField("val3", mesh.cells().size(), k_size);
  auto [val4, val4_view] = MakeAtlasField("val4", mesh.cells().size(), k_size);
  auto [val5, val5_view] = MakeAtlasField("val5", mesh.cells().size(), k_size);
  auto [val6, val6_view] = MakeAtlasField("val6", mesh.cells().size(), k_size);
  const auto N = r_view.numElements();

  // initialize fields
  for (int idx = 0; idx < N; idx++) {

    picture
      >> r_view(idx, level)
      >> g_view(idx, level)
      >> b_view(idx, level)
      ;

    val1_view(idx, level) = -1;
    val2_view(idx, level) = -1;
    val3_view(idx, level) = -1;
    val4_view(idx, level) = -1;
    val5_view(idx, level) = -1;
    val6_view(idx, level) = -1;

    //const auto [x, y] = wrapper.nodeLocation(idx);
    const auto [x, y] = wrapper.cellMidpoint(mesh, idx);
    x_view(idx, level) = x;
    y_view(idx, level) = y;
  }

  dawn_generated::cxxnaiveico::image_stencil<atlasInterface::atlasTag>(
      mesh,
      k_size,
      r_view,
      g_view,
      b_view,
      val1_view,
      val2_view,
      val3_view,
      val4_view,
      val5_view,
      val6_view,
      x_view,
      y_view
    )
      .run();

  dumpMesh4Triplot(mesh, "out/mesh", wrapper);
  dumpCellField("out/r.txt", mesh, r_view, level);
  dumpCellField("out/g.txt", mesh, g_view, level);
  dumpCellField("out/b.txt", mesh, b_view, level);

  printf("Finished image stencil successfully.\n");

#endif

  exit(0);
}
