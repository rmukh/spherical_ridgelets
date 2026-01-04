#include "SphericalRidgelets/include/CLIParser.h"
#include "cxxopts.hpp"
#include <fstream>
#include <iostream>

using namespace std;

static bool file_exists(const string &path) {
    if (path.empty()) return false;
    ifstream f(path.c_str());
    return f.good();
}

namespace ridgelets {

ParseResult parse_cli(int argc, char** argv) {
    ParseResult res;
    res.ok = false;

    DATA_SOURCE::input_parse p;
    // Defaults (match previous defaults)
    p.input_dmri = "";
    p.input_mask = "";
    p.output_ridgelets = "";
    p.output_fiber_max_odf = "";
    p.output_odf = "";
    p.signal_recon = "";
    p.output_A = "";
    p.external_gradients = "";
    p.ext_signal_recon = "";
    p.max_odf_thresh = 0.7;
    p.fista_lambda = 0.01;
    p.sph_rho = 3.125;
    p.lvl = 4;
    p.sph_J = 2;
    p.n_splits = -1;
    p.nth = -1;
    p.is_compress = false;
    p.fista_iterations = 2000;
    p.fista_tolerance = 0.001;

    try {
        cxxopts::Options options("sphridg", "Spherical Ridgelets processing");
        options
            .allow_unrecognised_options()
            .add_options()
            ("i,input", "Input dMRI file", cxxopts::value<string>())
            ("m,mask", "Mask file", cxxopts::value<string>())
            ("ridg", "Output ridgelet file", cxxopts::value<string>())
            ("sr", "Signal reconstruction output", cxxopts::value<string>())
            ("ext_sr", "External gradients signal reconstruction", cxxopts::value<string>())
            ("odf", "ODF values output", cxxopts::value<string>())
            ("omd", "ODF maxima dirs & values output", cxxopts::value<string>())
            ("A", "A basis output (matrix file)", cxxopts::value<string>())
            ("c,compress", "Enable compression", cxxopts::value<bool>()->default_value("false"))
            ("sj", "Spherical ridgelets J (integer)", cxxopts::value<int>()->default_value(to_string(p.sph_J)))
            ("srho", "Spherical ridgelets rho (float)", cxxopts::value<double>()->default_value(to_string(p.sph_rho)))
            ("lvl", "Icosahderon tessellation order (int)", cxxopts::value<int>()->default_value(to_string(p.lvl)))
            ("nspl", "Splits coefficient (int)", cxxopts::value<int>()->default_value(to_string(p.n_splits)))
            ("nth", "Number of threads", cxxopts::value<int>()->default_value(to_string(p.nth)))
            ("ext_grads", "External gradients file", cxxopts::value<string>())
            ("fi", "FISTA iterations", cxxopts::value<int>()->default_value(to_string(p.fista_iterations)))
            ("ft", "FISTA tolerance", cxxopts::value<double>()->default_value(to_string(p.fista_tolerance)))
            ("lmd", "FISTA lambda", cxxopts::value<double>()->default_value(to_string(p.fista_lambda)))
            ("mth", "Maxima ODF threshold", cxxopts::value<double>()->default_value(to_string(p.max_odf_thresh)))
            ("help", "Print help");

        auto result = options.parse(argc, argv);

        if (result.count("help") || !result.count("input")) {
            cout << options.help() << endl;
            res.ok = false;
            res.error_message = "";
            return res;
        }

        // required
        p.input_dmri = result["input"].as<string>();

        if (result.count("mask")) p.input_mask = result["mask"].as<string>();
        if (result.count("ridg")) p.output_ridgelets = result["ridg"].as<string>();
        if (result.count("sr")) p.signal_recon = result["sr"].as<string>();
        if (result.count("ext_sr")) p.ext_signal_recon = result["ext_sr"].as<string>();
        if (result.count("odf")) p.output_odf = result["odf"].as<string>();
        if (result.count("omd")) p.output_fiber_max_odf = result["omd"].as<string>();
        if (result.count("A")) p.output_A = result["A"].as<string>();
        if (result.count("compress")) p.is_compress = result["compress"].as<bool>();

        if (result.count("sj")) p.sph_J = static_cast<unsigned>(result["sj"].as<int>());
        if (result.count("srho")) p.sph_rho = static_cast<float>(result["srho"].as<double>());
        if (result.count("lvl")) p.lvl = static_cast<unsigned>(result["lvl"].as<int>());
        if (result.count("nspl")) p.n_splits = result["nspl"].as<int>();
        if (result.count("nth")) p.nth = result["nth"].as<int>();
        if (result.count("ext_grads")) p.external_gradients = result["ext_grads"].as<string>();
        if (result.count("fi")) p.fista_iterations = result["fi"].as<int>();
        if (result.count("ft")) p.fista_tolerance = static_cast<precisionType>(result["ft"].as<double>());
        if (result.count("lmd")) p.fista_lambda = static_cast<precisionType>(result["lmd"].as<double>());
        if (result.count("mth")) p.max_odf_thresh = static_cast<precisionType>(result["mth"].as<double>());

        // Simple validation
        if (!file_exists(p.input_dmri)) {
            res.ok = false;
            res.error_message = "Input dMRI file not found: " + p.input_dmri;
            return res;
        }
        if (!p.input_mask.empty() && !file_exists(p.input_mask)) {
            res.ok = false;
            res.error_message = "Mask file not found: " + p.input_mask;
            return res;
        }

        res.ok = true;
        res.params = p;
        return res;
    } catch (const cxxopts::OptionException &e) {
        res.ok = false;
        res.error_message = string("Error parsing options: ") + e.what();
        return res;
    }
}

} // namespace ridgelets
