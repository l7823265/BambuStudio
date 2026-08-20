#include <brepslicer/Slicer.h>

#include <cstdlib>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using brepslicer::SliceOptions;
using brepslicer::Vec3;

static std::vector<std::string> split(const std::string& s, char sep) {
    std::vector<std::string> out;
    std::string cur;
    for (char c : s) {
        if (c == sep) {
            out.push_back(cur);
            cur.clear();
        } else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

static Vec3 parseVec(const std::string& s) {
    const auto p = split(s, ',');
    if (p.size() != 3) throw std::runtime_error("expected x,y,z got: " + s);
    return {std::stod(p[0]), std::stod(p[1]), std::stod(p[2])};
}

static void usage() {
    std::cerr
        << "Usage:\n"
        << "  brep_slice model.step --normal 0,0,1 --layer-height 0.05 --out slices.json\n"
        << "Options:\n"
        << "  --normal x,y,z         slice plane normal (default 0,0,1)\n"
        << "  --layer-height h       layer thickness\n"
        << "  --start-height z       first plane offset n·p (default: bbox mid-first layer)\n"
        << "  --layers n             number of layers\n"
        << "  --heights z0,z1,...    explicit offsets (overrides height/count)\n"
        << "  --tolerance t          stitch / JSON tolerance (default 1e-4)\n"
        << "  --out file.json        output path (default slices.json)\n"
        << "  --svg dir              one SVG per layer (arcs use SVG A)\n"
        << "  --dxf path             one .dxf (ARC/CIRCLE), or a directory of per-layer DXF\n";
}

int main(int argc, char** argv) {
    try {
        if (argc < 2) {
            usage();
            return 1;
        }
        std::string step;
        std::string out = "slices.json";
        SliceOptions opt;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            auto need = [&](const char* name) -> std::string {
                if (i + 1 >= argc) throw std::runtime_error(std::string("missing value for ") + name);
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
            else if (a == "--heights") {
                for (const auto& h : split(need("--heights"), ',')) {
                    if (!h.empty()) opt.explicit_heights.push_back(std::stod(h));
                }
            } else if (a == "--tolerance") opt.tolerance = std::stod(need("--tolerance"));
            else if (a == "--out") out = need("--out");
            else if (a == "--svg") opt.svg_dir = need("--svg");
            else if (a == "--dxf") opt.dxf_path = need("--dxf");
            else throw std::runtime_error("unknown option: " + a);
        }
        if (step.empty()) throw std::runtime_error("missing STEP input");

        auto result = brepslicer::sliceFile(step, opt);
        brepslicer::writeJsonFile(result, out);
        std::cout << "wrote " << out << "  layers=" << result.layers.size() << "\n";
        if (!opt.svg_dir.empty()) std::cout << "wrote SVG in " << opt.svg_dir << "\n";
        if (!opt.dxf_path.empty()) std::cout << "wrote DXF " << opt.dxf_path << "\n";
        for (const auto& log : result.logs) std::cerr << log << "\n";
        return 0;
    } catch (const std::exception& ex) {
        std::cerr << "error: " << ex.what() << "\n";
        usage();
        return 1;
    }
}
