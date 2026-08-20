#pragma once

#include <brepslicer/Algo.h>

namespace brepslicer {

class OccKernel;
class OccShapeBuilder;
class OccSectionRef;

class OccStrategy : public IEngineStrategy {
public:
    IBrepKernel* GetKernel() override;
    IShapeBuilder* GetShapeBuilder() override;
    ISectionRef* GetSectionRef() override;
};

}  // namespace brepslicer
