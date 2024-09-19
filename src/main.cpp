#include "geometrycentral/surface/integer_coordinates_intrinsic_triangulation.h"
#include "geometrycentral/surface/manifold_surface_mesh.h"
#include "geometrycentral/surface/meshio.h"
#include "geometrycentral/surface/vertex_position_geometry.h"

#include "polyscope/curve_network.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#include "args/args.hxx"
#include "imgui.h"

#include <fstream>

#include "utils.h"

using namespace geometrycentral;
using namespace geometrycentral::surface;

// == Geometry-central data
std::unique_ptr<ManifoldSurfaceMesh> mesh;
std::unique_ptr<VertexPositionGeometry> geom;
std::unique_ptr<IntegerCoordinatesIntrinsicTriangulation> intTri;

// Poly scope visualization handle, to quickly add data to the surface
polyscope::SurfaceMesh* psMesh;

void updateTriagulationViz() {
    // Get the edge traces
    EdgeData<std::vector<SurfacePoint>> traces =
        intTri->traceAllIntrinsicEdgesAlongInput();

    // Convert to 3D positions
    std::vector<Vector3> traces3Dp;
    std::vector<std::array<size_t, 2>> traces3De;
    std::vector<size_t> edgeId, edgeHighlight;
    size_t i = 0;
    for (Edge e : intTri->mesh.edges()) {
        size_t N = traces3Dp.size();
        for (size_t iP = 0; iP + 1 < traces[e].size(); iP++) {
            traces3De.push_back({N + iP, N + iP + 1});
            edgeId.push_back(e.getIndex());
            edgeHighlight.push_back((e.getIndex() == 216) ? 1 : 0);
        }
        for (SurfacePoint& p : traces[e]) {
            traces3Dp.push_back(p.interpolate(geom->inputVertexPositions));
        }
        i++;
    }

    // Register with polyscope
    auto psEdges = polyscope::registerCurveNetwork("intrinsic edges", traces3Dp,
                                                   traces3De);
    psEdges->addEdgeScalarQuantity("id", edgeId);
    psEdges->addEdgeScalarQuantity("highlight", edgeHighlight);
}

// A user-defined callback, for creating control panels (etc)
// Use ImGUI commands to build whatever you want here, see
// https://github.com/ocornut/imgui/blob/master/imgui.h
void myCallback() {}

void doRandomFlips(IntegerCoordinatesIntrinsicTriangulation& intTri,
                   size_t nFlips = 500) {
    size_t iFlip = 0;
    while (iFlip < nFlips) {
        size_t iE = rand() % intTri.intrinsicMesh->nEdges();
        // std::cout << iFlip << ": \t" << iE << vendl;
        Edge eFlip = intTri.intrinsicMesh->edge(iE);

        bool flipped = intTri.flipEdgeIfPossible(eFlip);
        if (flipped) iFlip++;
    }
}

bool flipBackTopological(ManifoldSurfaceMesh& mesh, NormalCoordinates& n,
                         size_t* nFlipsOut = nullptr) {
    size_t nFlips = 0;

    for (Vertex v : mesh.vertices()) {
        bool done = false;
        while (!done) {
            done = true;
            for (Corner c : v.adjacentCorners()) {
                if (n.strictDegree(c) > 0) {
                    done            = false;
                    Halfedge heCurr = c.halfedge().next();
                    int p           = n.strictDegree(heCurr.corner());

                    bool curveEnded = false;
                    while (!curveEnded) {
                        Halfedge heToFlip = heCurr;
                        curveEnded        = n.stepTopologicalCurve(heCurr, p);
                        // Ideally, we would just move to heCurr after flipping.
                        // However, half edges transform in weird ways when
                        // performing edge flips, so instead we record whether
                        // it's the left or the right halfedge, and recompute it
                        // afterwards
                        bool wentRight = (heCurr == heToFlip.twin().next());

                        // flip rotates heToFlip *clockwise*
                        auto flipData = n.computeFlippedData(heToFlip.edge());
                        bool flipSuccess = mesh.flip(
                            heToFlip.edge(), false /* allow self edges */);
                        n.applyFlippedData(heToFlip.edge(), flipData);
                        nFlips++;

                        if (!flipSuccess) {
                            WATCH(heToFlip.edge());
                            if (nFlipsOut) *nFlipsOut = nFlips;
                            return false;
                        }

                        heCurr = (wentRight) ? heToFlip.twin().next().next()
                                             : heToFlip.next();
                    }
                }
                if (!done) break;
            }
        }
    }

    for (Edge e : mesh.edges()) {
        verbose_assert(n[e] < 0, "non-shared edge");
    }

    if (nFlipsOut) *nFlipsOut = nFlips;

    return true;
}

bool flipBack(IntegerCoordinatesIntrinsicTriangulation& intTri) {
    size_t nFlips = 0;

    NormalCoordinates& n = intTri.normalCoordinates;

    for (Vertex v : intTri.intrinsicMesh->vertices()) {
        bool done = false;
        while (!done) {
            done = true;
            for (Corner c : v.adjacentCorners()) {
                if (n.strictDegree(c) > 0) {
                    done            = false;
                    Halfedge heCurr = c.halfedge().next();
                    int p           = n.strictDegree(heCurr.corner());

                    bool curveEnded = false;
                    while (!curveEnded) {
                        Halfedge heToFlip = heCurr;
                        curveEnded        = n.stepTopologicalCurve(heCurr, p);
                        // Ideally, we would just move to heCurr after flipping.
                        // However, half edges transform in weird ways when
                        // performing edge flips, so instead we record whether
                        // it's the left or the right halfedge, and recompute it
                        // afterwards
                        bool wentRight = (heCurr == heToFlip.twin().next());

                        // flip rotates heToFlip *clockwise*
                        bool flipSuccess =
                            intTri.flipEdgeIfPossible(heToFlip.edge());

                        if (!flipSuccess) {
                            WATCH(heToFlip.edge());
                            return false;
                        }

                        heCurr = (wentRight) ? heToFlip.twin().next().next()
                                             : heToFlip.next();
                    }
                }
                if (!done) break;
            }
        }
    }

    return true;
}

int main(int argc, char** argv) {

    // Configure the argument parser
    args::ArgumentParser parser("Triangle Flipper");
    args::Positional<std::string> inputFilename(parser, "mesh",
                                                "Mesh to be processed.");

    std::cout << std::boolalpha; // print booleans as "true" or "false"

    // Parse args
    try {
        parser.ParseCLI(argc, argv);
    } catch (const args::Help&) {
        std::cout << parser;
        return 0;
    } catch (const args::ParseError& e) {
        std::cerr << e.what() << std::endl;
        std::cerr << parser;
        return 1;
    }

    std::string filename = "../../meshes/bunny_small.obj";
    // Make sure a mesh name was given
    if (inputFilename) {
        filename = args::get(inputFilename);
    }

    // Initialize polyscope
    polyscope::init();

    // Set the callback function
    polyscope::state::userCallback = myCallback;

    // Load mesh
    std::tie(mesh, geom) = readManifoldSurfaceMesh(filename);

    srand(0); // set seed to zero
    std::ofstream out("log.tsv");
    out << "nOrigFlips\tsum of normal "
           "coordinatess\tflips required\ttopologicalSuccess\tge"
           "ometricSuccess"
        << std::endl;
    std::vector<size_t> flipCounts{100,   250,   500,    1000,   2000,
                                   3000,  4000,  5000,   7500,   10000,
                                   20000, 50000, 100000, 1000000};
    for (size_t nF : flipCounts) {
        for (size_t rep = 0; rep < 10; rep++) {
            intTri.reset(
                new IntegerCoordinatesIntrinsicTriangulation(*mesh, *geom));

            doRandomFlips(*intTri, nF);
            if (rep + 1 == 10) updateTriagulationViz();

            // topological
            std::unique_ptr<ManifoldSurfaceMesh> intMesh =
                intTri->intrinsicMesh->copy();

            NormalCoordinates n(*intMesh);
            n.edgeCoords = intTri->normalCoordinates.edgeCoords;

            size_t normalCoordinateSum = 0;
            for (Edge e : intMesh->edges())
                normalCoordinateSum += fmax(n[e], 0);
            size_t nMosherFlips = 0;
            bool topologicalSuccess =
                flipBackTopological(*intMesh, n, &nMosherFlips);

            // geometric
            bool geometricSuccess = flipBack(*intTri);


            out << nF << "\t" << normalCoordinateSum << "\t" << nMosherFlips
                << "\t" << topologicalSuccess << "\t" << geometricSuccess
                << std::endl;
            std::cout << nF << "\t" << normalCoordinateSum << "\t"
                      << nMosherFlips << "\t" << topologicalSuccess << "\t"
                      << geometricSuccess << std::endl;
        }
    }

    if (false) {
        intTri.reset(
            new IntegerCoordinatesIntrinsicTriangulation(*mesh, *geom));

        srand(0); // set seed to zero
        doRandomFlips(*intTri, 5000);
        // intTri->flipEdgeIfPossible(intTri->intrinsicMesh->edge(1970));
        // intTri->flipEdgeIfPossible(intTri->intrinsicMesh->edge(1969));
        // intTri->flipEdgeIfPossible(intTri->intrinsicMesh->edge(123));
        // intTri->flipEdgeIfPossible(intTri->intrinsicMesh->edge(124));

        updateTriagulationViz();

        // WATCH(flipBack(*intTri));
        WATCH(flipBackTopological(*intTri->intrinsicMesh,
                                  intTri->normalCoordinates));
    }

    // Register the mesh with polyscope
    psMesh = polyscope::registerSurfaceMesh(
        polyscope::guessNiceNameFromPath(filename), geom->vertexPositions,
        mesh->getFaceVertexList(), polyscopePermutations(*mesh));

    // Give control to the polyscope gui
    polyscope::show();

    return EXIT_SUCCESS;
}
