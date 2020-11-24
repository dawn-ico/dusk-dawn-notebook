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

#include "intro_cxx-naive.cpp"

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

  // wrapper with various atlas helper functions
  AtlasToCartesian wrapper(mesh, lDomain, false, true);

  // mesh properties
  const int edgesPerVertex = 6;
  const int edgesPerNode = 6;
  const int edgesPerCell = 3;

  // system variables
  auto [a, a_view] = MakeAtlasField("a", mesh.edges().size(), k_size);
  const auto N = a_view.numElements();

  // initialize a with not `1`
  for (int idx = 0; idx < N; idx++) {
    a_view(idx, level) = -1;
  }

  dawn_generated::cxxnaiveico::simple_example_stencil<atlasInterface::atlasTag>(
      mesh, k_size, a_view)
      .run();

  dumpCellField("out/a.txt", mesh, a_view, level);

  for (int idx = 0; idx < N; idx++) {
    const auto val = a_view(idx, level);
    if (val != 1.0f) {
      printf("For some cell `a` == %f (should be 1.0)!\n", val);
      exit(-1);
    }
  }
  printf("As expected, `a` == 1 everywhere.\n");

  exit(0);
}
