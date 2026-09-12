#pragma once

#include <brepslicer/Shape.h>

#include <utility>
#include <vector>

namespace brepslicer {

// Manifold-ish edge↔face adjacency for one shape (all solids). Built by OCC once per
// model; used by assemble to link slice chains by shared-edge order instead of pure
// geometric nearest tip.
struct SolidAdjacency {
    // edge_id -> incident face_ids (usually 2; 1 = boundary; 3+ = non-manifold).
    std::vector<std::vector<int>> edge_faces;
    // face_id -> (edge_id, neighbor_face_id) in no particular order.
    // neighbor_face_id = -1 for boundary edges.
    std::vector<std::vector<std::pair<int, int>>> face_links;

    bool empty() const { return edge_faces.empty(); }

    bool facesShareEdge(int fa, int fb, int* edge_id = nullptr) const {
        if (fa < 0 || fb < 0 || fa == fb) return false;
        if (fa >= static_cast<int>(face_links.size())) return false;
        for (const auto& lk : face_links[static_cast<size_t>(fa)]) {
            if (lk.second == fb) {
                if (edge_id) *edge_id = lk.first;
                return true;
            }
        }
        return false;
    }

    // Faces adjacent to `face` (unique).
    std::vector<int> neighbors(int face) const {
        std::vector<int> out;
        if (face < 0 || face >= static_cast<int>(face_links.size())) return out;
        for (const auto& lk : face_links[static_cast<size_t>(face)]) {
            if (lk.second < 0) continue;
            bool dup = false;
            for (int x : out) {
                if (x == lk.second) {
                    dup = true;
                    break;
                }
            }
            if (!dup) out.push_back(lk.second);
        }
        return out;
    }
};

}  // namespace brepslicer
