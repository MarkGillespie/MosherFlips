bool flipBack(IntegerCoordinatesIntrinsicTriangulation& intTri) {
    size_t nFlips = 0;

    NormalCoordinates& n = intTri.normalCoordinates;

    std::deque<Vertex> verticesToConsider;
    for (Vertex v : intTri.intrinsicMesh->vertices()) {
        size_t nMissingEdges = 0;
        for (Corner c : v.adjacentCorners()) {
            nMissingEdges += n.strictDegree(c);
        }
        if (nMissingEdges > 0) verticesToConsider.push_back(v);
    }

    while (!verticesToConsider.empty()) {
        bool doneWithVertex = false;
        bool vertexFailed   = false;
        bool madeProgress   = true;
        while (!vertexFailed && !doneWithVertex) {
            Vertex v = verticesToConsider.front();
            verticesToConsider.pop_front();
            doneWithVertex = true;
            for (Corner c : v.adjacentCorners()) {
                if (n.strictDegree(c) > 0) {
                    doneWithVertex  = false;
                    Halfedge heCurr = c.halfedge().next();
                    int p           = n.strictDegree(heCurr.corner());

                    bool curveEnded = false;
                    while (!curveEnded) {
                        Halfedge heToFlip = heCurr;
                        curveEnded        = n.stepTopologicalCurve(heCurr, p);
                        // Ideally, we would just move to heCurr after
                        // flipping. However, half edges transform in weird
                        // ways when performing edge flips, so instead we
                        // record whether it's the left or the right
                        // halfedge, and recompute it afterwards
                        bool wentRight = (heCurr == heToFlip.twin().next());

                        // flip rotates heToFlip *clockwise*
                        bool flipSuccess =
                            intTri.flipEdgeIfPossible(heToFlip.edge());

                        if (!flipSuccess) {
                            verticesToConsider.emplace_back(v);
                            vertexFailed = true;
                            break;
                        }

                        heCurr = (wentRight) ? heToFlip.twin().next().next()
                                             : heToFlip.next();
                    }
                }
                if (vertexFailed || !doneWithVertex) break;
            }
        }
    }


    return true;
}
