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

#include "laplacian_fvm_cxx-naive.cpp"

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

void dumpEdgeFieldOnCells(const std::string &fname, const atlas::Mesh &mesh,
                          AtlasToCartesian wrapper,
                          atlasInterface::Field<double> &field, int level, std::vector<bool> boundary) {
  FILE *fp = fopen(fname.c_str(), "w+");
  const auto &conn = mesh.cells().edge_connectivity();
  for (int cIdx = 0; cIdx < mesh.cells().size(); cIdx++) {      
    double h = 0.;
    for (int nbhIdx = 0; nbhIdx < conn.cols(cIdx); nbhIdx++) {
      int eIdx = conn(cIdx, nbhIdx);
      assert(eIdx != conn.missing_value());
      h += std::isfinite(field(eIdx, 0)) ? field(eIdx, 0) : 0.;
    }
    h /= conn.cols(cIdx);
    if (!boundary[cIdx]) {
        fprintf(fp, "%f\n", h);
    } else {
        fprintf(fp, "%f\n", 0.);
    }
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

std::vector<double> readField(const std::string &fname) {
  std::ifstream ifile(fname, std::ios::in);
  if (!ifile.good())
      std::cout << "could not open " << fname << "\n";
  double num = 0.0;
  std::vector<double> ret;
  while (ifile >> num) {
    ret.push_back(num);
  }
  return ret;
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


int main(int argc, char *argv[]) {
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
  auto [uE_F, uE] = MakeAtlasField("u", mesh.edges().size(), k_size);
  auto [vE_F, vE] = MakeAtlasField("v", mesh.edges().size(), k_size);
  
  auto [uvC_div_F, uvC_div] =
      MakeAtlasField("uv_div", mesh.cells().size(), k_size);
  auto [uvN_curl_F, uvN_curl] =
      MakeAtlasField("uv_curl", mesh.nodes().size(), k_size);
  auto [uvE_nabla2_F, uvE_nabla2] =
      MakeAtlasField("uv_nabla2", mesh.edges().size(), k_size);
    
  auto [gard_of_curl_F, grad_of_curl] =
      MakeAtlasField("grad_of_curl", mesh.edges().size(), k_size);    
  auto [grad_of_div_F, grad_of_div] =
      MakeAtlasField("grad_of_div", mesh.edges().size(), k_size);        

  // fields to hold analytical solutions  
  auto [uvC_div_F_Sol, uvC_div_Sol] =
      MakeAtlasField("uv_div_Sol", mesh.cells().size(), k_size);
  auto [uvN_curl_F_Sol, uvN_curl_Sol] =
      MakeAtlasField("uv_curl_Sol", mesh.nodes().size(), k_size);
  auto [uvN_nabla2_F_Sol, uvE_nabla2_Sol] =
      MakeAtlasField("uv_nabla2_Sol", mesh.edges().size(), k_size);

  // Geometrical factors on edges
  auto [L_F, L] =
      MakeAtlasField("L", mesh.edges().size(), k_size); // edge length
  auto [dualL_F, dualL] =
      MakeAtlasField("dualL", mesh.edges().size(), k_size); // edge length
  auto [nx_F, nx] =
      MakeAtlasField("nx", mesh.edges().size(), k_size); // normals
  auto [ny_F, ny] = MakeAtlasField("ny", mesh.edges().size(), k_size);  
  auto [tangent_orientation_F, tangent_orientation] =
      MakeAtlasField("tangent_orientation", mesh.edges().size(), 1);    

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
  auto dot = [](const Vector &v1, const Vector &v2) {
    return std::get<0>(v1) * std::get<0>(v2) +
           std::get<1>(v1) * std::get<1>(v2);
  };
  
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    L(edgeIdx, level) = wrapper.edgeLength(mesh, edgeIdx);
    tangent_orientation(edgeIdx) = wrapper.tangentOrientation(mesh, edgeIdx);        
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
  auto analyticalLaplacian = [](double x, double y) -> std::tuple<double, double> {
    double c1 = 0.25 * sqrt(105. / (2 * M_PI));
    double c2 = 0.5 * sqrt(15. / (2 * M_PI));
    return {-4 * c1 * cos(2 * x) * cos(y) * cos(y) * sin(y), -4 * c2 * cos(x) * sin(y) * cos(y)};
  };
  

  // initialize test functions  
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
    auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [uu, vv] = sphericalHarmonic(x, y);
    uE(edgeIdx, level) = uu;
    vE(edgeIdx, level) = vv;
  }

  // collect analytical solutions
  for (int cellIdx = 0; cellIdx < mesh.cells().size(); cellIdx++) {
    auto [x, y] = wrapper.cellMidpoint(mesh, cellIdx);    
    double uv_div = analyticalDivergence(x, y);
    uvC_div_Sol(cellIdx, level) = uv_div;
  }
  for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
    auto [x, y] = wrapper.nodeLocation(nodeIdx);
    double uv_curl = analyticalCurl(x, y);
    uvN_curl_Sol(nodeIdx, level) = uv_curl;
  }
  for (int edgeIdx = 0; edgeIdx < mesh.edges().size(); edgeIdx++) {
      uvE_nabla2_Sol(edgeIdx, level) = 0.;
  }
  for (auto edgeIdx : wrapper.innerEdges(mesh)) {
    auto [x, y] = wrapper.edgeMidpoint(mesh, edgeIdx);
    auto [uu, vv] = analyticalLaplacian(x, y);
    uvE_nabla2_Sol(edgeIdx, level) = uu*nx(edgeIdx, level) + vv*ny(edgeIdx, level);
  }
    
  dumpMesh4Triplot(mesh, "out/mesh", wrapper);  
  dawn_generated::cxxnaiveico::laplacian_fvm<atlasInterface::atlasTag>(
        mesh, k_size, uE, vE, nx, ny, 
        uvC_div, uvN_curl, grad_of_curl, grad_of_div, uvE_nabla2, 
        L, dualL, A, dualA, tangent_orientation, edge_orientation_node, edge_orientation_cell)
        .run();
   {    
    auto [L1, L2, Linf] =
        compare(uvN_curl_Sol, uvN_curl, boundary_nodes, level);
    for (int nodeIdx = 0; nodeIdx < mesh.nodes().size(); nodeIdx++) {
        if (boundary_nodes[nodeIdx]) {
            uvN_curl(nodeIdx, level) = 0.;
            uvN_curl_Sol(nodeIdx, level) = 0.;
        }
    }
              
    printf("curl_uv L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);
    if (L1 < 0.002 && L2 < 0.002 && Linf < 0.008) {
        printf("which looks right! :)\n");
    } else {
        printf("which looks off, please recheck your curl computation. You can copy the one from the previous exercise\n");
    }
    dumpNodeFieldOnCells("out/uv_curl.txt", mesh, uvN_curl, level);
    dumpNodeFieldOnCells("out/uv_curl_Sol.txt", mesh, uvN_curl_Sol, level);
   }

    
   {   
    auto [L1, L2, Linf] = compare(uvC_div_Sol, uvC_div, boundary_cells, level);
    printf("div_uv L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);
    if (L1 < 0.02 && L2 < 0.02 && Linf < 0.06) {
        printf("which looks right! :)\n");
    } else {
        printf("which looks off, please recheck your divergence computation. You can copy the one from the previous exercise\n");
    }
    dumpCellField("out/uv_div.txt", mesh, uvC_div, level);
    dumpCellField("out/uv_div_Sol.txt", mesh, uvC_div_Sol, level);    
  }
   
  {
    for (int eIdx = 0; eIdx < mesh.edges().size(); eIdx++) {
        if (!std::isfinite(uvE_nabla2(eIdx, level))) {
            uvE_nabla2(eIdx, level) = uvE_nabla2_Sol(eIdx,level);
        }
    }
    auto innerEdges = wrapper.innerEdges(mesh);
      for(int eIdx = 0; eIdx < mesh.edges().size(); eIdx++) {
        if(std::find(innerEdges.begin(), innerEdges.end(), eIdx) == innerEdges.end()) {
          uvE_nabla2(eIdx, 0) = uvE_nabla2_Sol(eIdx, 0);
        }
      }
    //dumpEdgeField("ref/nabla2_ref.txt", mesh, uvE_nabla2, level);        
    
    auto [L1, L2, Linf] = compare(uvE_nabla2_Sol, uvE_nabla2, boundary_edges, level);        
    printf("nabal2_uv L1: %f, L2: %f, Linf: %f\n", L1, L2, Linf);
    printf("which may or may not be correct, checking aginst reference...");
    auto nabla2ref = readField("ref/nabla2_ref.txt");
    bool correct = compare(nabla2ref, uvE_nabla2, level);
    if (correct) {
        printf("...check against reference shows that this is indeed correct! :)");
    } else {
        printf("...check against reference unfortunately shows that something is off");
    }
    dumpEdgeFieldOnCells("out/uv_nabla2.txt", mesh, wrapper, uvE_nabla2, level, boundary_cells);    
    dumpEdgeFieldOnCells("out/uv_nabla2_Sol.txt", mesh, wrapper, uvE_nabla2_Sol, level, boundary_cells);        
  }
  
  return 0;

  
}