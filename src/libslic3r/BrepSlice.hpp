#ifndef slic3r_BrepSlice_hpp_
#define slic3r_BrepSlice_hpp_

#include "ExPolygon.hpp"
#include "ObjectID.hpp"
#include "Point.hpp"

#include <functional>
#include <string>
#include <vector>

namespace Slic3r {

class ModelObject;

struct BrepVolumeSlices {
    ObjectID                  volume_id;
    std::vector<ExPolygons>   layers;
};

// True when this object was imported from a STEP file that is still on disk,
// so BrepSlicer can cut the original B-rep instead of the tessellated mesh.
bool model_object_has_brep_step(const ModelObject &object);

// Slice each model-part volume of a STEP object at the given object-space Z
// values (same Zs the mesh slicer would use). On success, out_by_volume holds
// one entry per sliced model part with matching volume_id.
// Returns false to fall back to mesh slicing (missing file, modifiers, etc.).
bool slice_model_object_brep(const ModelObject             &object,
                             const Transform3d             &trafo_centered,
                             const std::vector<float>      &zs,
                             double                         chord_error,
                             const std::function<void()>   &throw_on_cancel,
                             std::vector<BrepVolumeSlices> &out_by_volume);

} // namespace Slic3r

#endif
