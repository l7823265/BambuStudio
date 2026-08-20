#include <occ/OccStrategy.h>
#include <brepslicer/Engine.h>

#include <memory>

namespace brepslicer {

// Defined in the other OCC translation units.
IBrepKernel* makeOccKernel();
IShapeBuilder* makeOccBuilder();
ISectionRef* makeOccSectionRef();

IBrepKernel* OccStrategy::GetKernel() { return makeOccKernel(); }
IShapeBuilder* OccStrategy::GetShapeBuilder() { return makeOccBuilder(); }
ISectionRef* OccStrategy::GetSectionRef() { return makeOccSectionRef(); }

#define USING_OCCENGINE

#ifdef USING_OCCENGINE
std::unique_ptr<IEngineStrategy> ShapeEngine::strategy_ = std::make_unique<OccStrategy>();
#endif

std::unique_ptr<IBrepKernel> ShapeEngine::kernel_{strategy_->GetKernel()};
std::unique_ptr<IShapeBuilder> ShapeEngine::builder_{strategy_->GetShapeBuilder()};
std::unique_ptr<ISectionRef> ShapeEngine::section_{strategy_->GetSectionRef()};

IBrepKernel& ShapeEngine::Kernel() { return *kernel_; }
IShapeBuilder& ShapeEngine::Builder() { return *builder_; }
ISectionRef& ShapeEngine::SectionRef() { return *section_; }

}  // namespace brepslicer
