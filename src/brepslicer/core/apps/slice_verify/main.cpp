#include <brepslicer/Slicer.h>

#include <iostream>
#include <stdexcept>
#include <string>

static void usage() {
    std::cerr
        << "Usage:\n"
        << "  slice_verify model.step --normal 0,0,1 --layer-height 0.05 [--out slices.json] [--svg dir]\n";
}

static brepslicer::Vec3 parseVec(const std::string& s) {
    const auto c1 = s.find(',');
    const auto c2 = s.find(',', c1 + 1);
    if (c1 == std::string::npos || c2 == std::string::npos) {
        throw std::runtime_error("expected x,y,z");
    }
    return {std::stod(s.substr(0, c1)), std::stod(s.substr(c1 + 1, c2 - c1 - 1)),
            std::stod(s.substr(c2 + 1))};
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            usage();
            return 1;
        }
        std::string step;
        std::string out;
        brepslicer::SliceOptions opt;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto need = [&](const char*) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error("missing value");
                return argv[++i];
            };
            if (a == "-h" || a == "--help") {
                usage();
                return 0;
            }
            if (!a.empty() && a[0] != '-') {
                step = a;
                continue;
            }
            if (a == "--normal") opt.normal = parseVec(need("--normal"));
            else if (a == "--layer-height") opt.layer_height = std::stod(need("--layer-height"));
            else if (a == "--start-height") {
                opt.start_height = std::stod(need("--start-height"));
                opt.start_height_set = true;
            } else if (a == "--layers") opt.layer_count = std::stoi(need("--layers"));
            else if (a == "--tolerance") opt.tolerance = std::stod(need("--tolerance"));
            else if (a == "--out") out = need("--out");
            else if (a == "--svg") opt.svg_dir = need("--svg");
            else throw std::runtime_error("unknown option: " + a);
        }
        if (step.empty()) throw std::runtime_error("missing STEP input");

        const auto shape = brepslicer::readStep(step);
        auto result = brepslicer::sliceShape(shape, opt);
        if (!out.empty()) brepslicer::writeJsonFile(result, out);

        const auto report = brepslicer::verifyAgainstRef(result, shape);
        for (const auto& e : report.errors) std::cerr << "ERROR: " << e << "\n";
        for (const auto& w : report.warnings) std::cerr << "WARN:  " << w << "\n";
        for (const auto& n : report.notes) std::cout << n << "\n";
        std::cout << (report.ok ? "PASS" : "FAIL") << "  layers=" << result.layers.size() << "\n";
        return report.ok ? 0 : 2;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        usage();
        return 1;
    }
}
