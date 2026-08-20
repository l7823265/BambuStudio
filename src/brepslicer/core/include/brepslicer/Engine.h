#pragma once

#include <brepslicer/Algo.h>

#include <memory>

namespace brepslicer {

class ShapeEngine {
public:
    static IBrepKernel& Kernel();
    static IShapeBuilder& Builder();
    static ISectionRef& SectionRef();

private:
    static std::unique_ptr<IEngineStrategy> strategy_;
    static std::unique_ptr<IBrepKernel> kernel_;
    static std::unique_ptr<IShapeBuilder> builder_;
    static std::unique_ptr<ISectionRef> section_;
};

}  // namespace brepslicer
