#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>
#include <vector>

#include "GpuHydroOpenCL.h"
#include "HydroProfile.h"

namespace {

struct Options {
    int samples = 1000000;
    int repeats = 5;
    int ebe_hydro = 0;
    std::string cent = "0-5";
};

void printUsage(const char *prog) {
    std::cerr << "Usage: " << prog
              << " [--samples N] [--repeats N] [--ebe-hydro 0|1] [--cent C]\n"
              << "Run from a directory containing the hydro files expected by HydroProfile.\n";
}

bool parseInt(const char *s, int &out) {
    char *end = nullptr;
    long v = std::strtol(s, &end, 10);
    if (!end || *end != '\0') return false;
    out = static_cast<int>(v);
    return true;
}

bool parseArgs(int argc, char **argv, Options &opt) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto requireValue = [&](const char *name) -> const char* {
            if (i + 1 >= argc) {
                std::cerr << "Missing value for " << name << "\n";
                return nullptr;
            }
            return argv[++i];
        };
        if (arg == "--samples") {
            const char *value = requireValue("--samples");
            if (!value || !parseInt(value, opt.samples)) return false;
        } else if (arg == "--repeats") {
            const char *value = requireValue("--repeats");
            if (!value || !parseInt(value, opt.repeats)) return false;
        } else if (arg == "--ebe-hydro") {
            const char *value = requireValue("--ebe-hydro");
            if (!value || !parseInt(value, opt.ebe_hydro)) return false;
        } else if (arg == "--cent") {
            const char *value = requireValue("--cent");
            if (!value) return false;
            opt.cent = value;
        } else if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            std::exit(0);
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return false;
        }
    }
    return opt.samples > 0 && opt.repeats > 0;
}

double checksum(const std::vector<double> &a,
                const std::vector<double> &b,
                const std::vector<double> &c) {
    double sum = 0.0;
    const size_t stride = std::max<size_t>(1, a.size() / 10000);
    for (size_t i = 0; i < a.size(); i += stride) {
        sum += a[i] * 1.0 + b[i] * 0.25 + c[i] * 0.125;
    }
    return sum;
}

double maxAbsDiff(const std::vector<double> &a, const std::vector<double> &b) {
    double diff = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        diff = std::max(diff, std::abs(a[i] - b[i]));
    }
    return diff;
}

}  // namespace

int main(int argc, char **argv) {
    Options opt;
    if (!parseArgs(argc, argv, opt)) {
        printUsage(argv[0]);
        return 2;
    }

    HydroProfile hydro;
    hydro.loadHydro(opt.ebe_hydro, opt.cent);
    if (hydro.hydroT().empty()) {
        std::cerr << "Hydro profile did not load any cells\n";
        return 1;
    }

    std::vector<double> tau(opt.samples);
    std::vector<double> x(opt.samples);
    std::vector<double> y(opt.samples);
    std::mt19937_64 rng(12345);
    std::uniform_real_distribution<double> tau_dist(hydro.hydroTau0(),
                                                    std::max(hydro.hydroTau0(), hydro.hydroTauMax() - hydro.hydroDtau()));
    std::uniform_real_distribution<double> xy_dist(-0.95 * hydro.hydroXmax(), 0.95 * hydro.hydroXmax());
    for (int i = 0; i < opt.samples; ++i) {
        tau[i] = tau_dist(rng);
        x[i] = xy_dist(rng);
        y[i] = xy_dist(rng);
    }

    std::vector<double> cpu_t(opt.samples), cpu_vx(opt.samples), cpu_vy(opt.samples);
    double cpu_ms_total = 0.0;
    double cpu_sum = 0.0;
    for (int r = 0; r < opt.repeats; ++r) {
        const auto start = std::chrono::steady_clock::now();
        for (int i = 0; i < opt.samples; ++i) {
            hydro.getValues(tau[i], x[i], y[i], cpu_t[i], cpu_vx[i], cpu_vy[i]);
        }
        const auto stop = std::chrono::steady_clock::now();
        cpu_ms_total += std::chrono::duration<double, std::milli>(stop - start).count();
        cpu_sum += checksum(cpu_t, cpu_vx, cpu_vy);
    }

    GpuHydroOpenCL gpu;
    std::string error;
    if (!gpu.initialize(hydro, error)) {
        std::cerr << "OpenCL GPU initialization failed: " << error << "\n";
        return 3;
    }

    std::vector<double> gpu_t, gpu_vx, gpu_vy;
    double gpu_ms_total = 0.0;
    double gpu_sum = 0.0;
    std::string device;
    for (int r = 0; r < opt.repeats; ++r) {
        GpuHydroTiming timing;
        if (!gpu.interpolate(tau, x, y, gpu_t, gpu_vx, gpu_vy, &timing, error)) {
            std::cerr << "OpenCL GPU interpolation failed: " << error << "\n";
            return 4;
        }
        gpu_ms_total += timing.milliseconds;
        gpu_sum += checksum(gpu_t, gpu_vx, gpu_vy);
        device = timing.device_name;
    }

    const double cpu_ms = cpu_ms_total / opt.repeats;
    const double gpu_ms = gpu_ms_total / opt.repeats;
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "samples " << opt.samples << "\n";
    std::cout << "repeats " << opt.repeats << "\n";
    std::cout << "device " << device << "\n";
    std::cout << "cpu_ms_per_repeat " << cpu_ms << "\n";
    std::cout << "gpu_ms_per_repeat " << gpu_ms << "\n";
    std::cout << "speedup_cpu_over_gpu " << (gpu_ms > 0.0 ? cpu_ms / gpu_ms : 0.0) << "\n";
    std::cout << "max_abs_diff_temp " << maxAbsDiff(cpu_t, gpu_t) << "\n";
    std::cout << "max_abs_diff_vx " << maxAbsDiff(cpu_vx, gpu_vx) << "\n";
    std::cout << "max_abs_diff_vy " << maxAbsDiff(cpu_vy, gpu_vy) << "\n";
    std::cout << "cpu_checksum " << cpu_sum << "\n";
    std::cout << "gpu_checksum " << gpu_sum << "\n";
    return 0;
}
