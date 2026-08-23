#pragma once

#include <brepslicer/Types.h>
#include <brepslicer/Shape.h>
#include <brepslicer/Engine.h>

#include <memory>
#include <string>

namespace brepslicer {

std::shared_ptr<IShape> readStep(const std::string& path);

SliceResult sliceShape(const std::shared_ptr<IShape>& shape, const SliceOptions& opt);
SliceResult sliceFile(const std::string& stepPath, const SliceOptions& opt);

std::string writeJson(const SliceResult& result);
void writeJsonFile(const SliceResult& result, const std::string& path);

void writeLayerSvg(const Layer& layer, const Vec3& normal, const std::string& path);
void writeAllLayerSvg(const SliceResult& result, const std::string& dir);

void writeDxfFile(const SliceResult& result, const std::string& path);
void writeAllLayerDxf(const SliceResult& result, const std::string& path);
void writeOpenLayerDxf(const SliceResult& result,
                       const std::vector<std::pair<double, int>>& unclosed_layers,
                       const std::string& dir);
void writeSeedDxfFile(const SliceResult& result, const std::string& path);
void writeSeedSummaryFile(const SliceResult& result, const std::string& path);

VerifyReport verifySlice(const SliceResult& result);
VerifyReport verifyAgainstRef(const SliceResult& result, const std::shared_ptr<IShape>& shape);

}  // namespace brepslicer
