#include <brepslicer/Algo.h>
#include <topo/SolidAdjacency.h>

namespace brepslicer {

SolidAdjacency IBrepKernel::buildSolidAdjacency(const std::vector<FaceRecord>&) const {
    return {};
}

}  // namespace brepslicer
