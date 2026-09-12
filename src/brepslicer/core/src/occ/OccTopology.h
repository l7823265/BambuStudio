#pragma once

#include <topo/SolidAdjacency.h>

#include <memory>
#include <vector>

namespace brepslicer {

// Build edge↔face adjacency from OCC TopoDS (IsSame edge map). Returns empty on failure.
SolidAdjacency buildSolidAdjacency(const std::vector<FaceRecord>& faces);

}  // namespace brepslicer
