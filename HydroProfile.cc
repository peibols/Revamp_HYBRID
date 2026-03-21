#include "HydroProfile.h"

#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

HydroProfile::HydroProfile() = default;
HydroProfile::~HydroProfile() = default;

size_t HydroProfile::index(int it, int iy, int ix, int ietamax, int ixmax) {
    return (static_cast<size_t>(it) * static_cast<size_t>(ietamax) + static_cast<size_t>(iy)) * static_cast<size_t>(ixmax) + static_cast<size_t>(ix);
}

int HydroProfile::clampIndex(int idx, int maxCount) {
    if (idx < 0) return 0;
    const int maxIdx = std::max(1, maxCount - 1);
    if (idx > maxIdx - 1) return maxIdx - 1;
    return idx;
}

double HydroProfile::fracFromIndex(double coord, double origin, double spacing, int idx) {
    return (coord - (origin + static_cast<double>(idx) * spacing)) / spacing;
}

void HydroProfile::loadHydro(int mode, const std::string &cent) {
    mode_ = mode;
    if (mode == 0) {
        tempScalingFactor_ = 0.2;
        loadPlaintextHydro(cent);
    } else {
        tempScalingFactor_ = 1.0;
        loadIpsatBinary(cent, "evolution_all_xyeta.dat");
    }
}

void HydroProfile::loadIpsatBinary(const std::string &cent, const std::string &filename) {
    std::ostringstream oss;
    oss << filename;

    FILE *hydro = std::fopen(oss.str().c_str(), "rb");
    if (!hydro) {
        std::cerr << "Hydro open fail = " << oss.str() << "\n";
        std::exit(1);
    }

    std::cout << "Reading Hydro..." << std::endl;

    float header[16];
    int status = std::fread(&header, sizeof(float), 16, hydro);
    if (status == 0) {
        std::cerr << "[HydroProfile::loadIpsat]: ERROR: Can not read the evolution file header" << std::endl;
        std::exit(1);
    }

    hydroTau0_ = double(header[0]);
    hydroDtau_ = double(header[1]);
    ixmax_ = static_cast<int>(header[2]);
    hydroDx_ = double(header[3]);
    hydroXmax_ = double(std::abs(header[4]));
    ietamax_ = static_cast<int>(header[8]);
    // header[9-10] are eta spacing / max, currently unused in this implementation

    float cell_info[16];
    int itau_max = 0;
    double maxtemp = 0.0;

    // Pre-allocate with at least one element in each dimension.
    itaumax_ = 1;
    hydrot_.clear();
    hydrox_.clear();
    hydroy_.clear();

    while (true) {
        status = std::fread(&cell_info, sizeof(float), 16, hydro);
        if (status == 0) break;
        if (status != 16) {
            std::cerr << "[HydroProfile::loadIpsat]: ERROR: the evolution file format is not correct" << std::endl;
            std::exit(1);
        }

        itau_max = std::max(itau_max, static_cast<int>(cell_info[0]));
        int itau = static_cast<int>(cell_info[0]);
        int ix = static_cast<int>(cell_info[1]);
        int iy = static_cast<int>(cell_info[2]);

        float temperature = cell_info[6];
        float ux = cell_info[8];
        float uy = cell_info[9];
        float uz = cell_info[10];

        float u2 = ux * ux + uy * uy + uz * uz;
        float gamma = std::sqrt(1 + u2);
        float vx = ux / gamma;
        float vy = uy / gamma;

        // Ensure we have enough storage for this time slice
        int requiredTime = itau + 1;
        if (requiredTime > itaumax_) {
            itaumax_ = requiredTime;
            size_t totalSize = static_cast<size_t>(itaumax_) * static_cast<size_t>(ietamax_) * static_cast<size_t>(ixmax_);
            hydrot_.resize(totalSize, 0.0);
            hydrox_.resize(totalSize, 0.0);
            hydroy_.resize(totalSize, 0.0);
        }

        size_t idx = index(itau, iy, ix, ietamax_, ixmax_);
        hydrot_[idx] = double(temperature);
        hydrox_[idx] = double(vx);
        hydroy_[idx] = double(vy);

        maxtemp = std::max(maxtemp, static_cast<double>(temperature));
    }

    hydroTauMax_ = hydroTau0_ + hydroDtau_ * (itaumax_ - 1);

    std::cout << "Read hydro: itaumax=" << itaumax_ << " max temp=" << maxtemp << "\n";
    std::fclose(hydro);
}

void HydroProfile::loadPlaintextHydro(const std::string &cent) {
    // Event-averaged hydro: read plaintext file ./hydroinfoPlaintxtHuichaoFormat.dat
    // Format: hor ver tou enedat tdat vxdat vydat (horizontal, vertical, time, energy_density, temperature, vel_x, vel_y)
    // Fixed parameters (boost-invariant: same temperature/velocity for all eta/rapidity)
    const int maxx = 100;
    const int maxy = 100;
    const double tau0 = 0.6;
    const double deltat = 0.1;
    const double deltax = 0.3;
    const double deltay = 0.3;
    
    std::string hydFile = "./hydroinfoPlaintxtHuichaoFormat.dat";
    std::ifstream hydro(hydFile);
    if (!hydro.is_open()) {
        std::cerr << "Hydro plaintext file open fail: " << hydFile << "\n";
        std::exit(1);
    }

    std::cout << "Reading plaintext Hydro from " << hydFile << "..." << std::endl;
    std::cout.flush();
    
    // Initialize grid parameters
    hydroTau0_ = tau0;
    hydroDtau_ = deltat;
    hydroDx_ = deltax;
    hydroXmax_ = deltax * maxx / 2.0;
    ixmax_ = maxx + 1;
    ietamax_ = maxy + 1;  // The "eta" dimension is really the y-dimension for plaintext
    itaumax_ = 0;  // Start at 0, will allocate on first data point
    
    hydrot_.clear();
    hydrox_.clear();
    hydroy_.clear();

    double hor, ver, tou, enedat, tdat, vxdat, vydat;
    double maxtemp = 0.0;
    
    int line_count = 0;

    while (hydro >> hor >> ver >> tou >> enedat >> tdat >> vxdat >> vydat) {
        line_count++;

        // Discretize continuous coordinates into grid indices
        int it = static_cast<int>((tou + deltat / 2.0 - tau0) / deltat);
        int ix = static_cast<int>((hor + deltax * maxx / 2.0 + deltax / 2.0) / deltax);
        int iy = static_cast<int>((ver + deltay * maxy / 2.0 + deltay / 2.0) / deltay);

        // Bounds checking
        if (it < 0 || ix < 0 || ix > maxx || iy < 0 || iy > maxy) {
            continue;  // Skip out-of-bounds points
        }

        maxtemp = std::max(maxtemp, tdat);

        // Ensure we have enough storage
        int requiredTime = it + 1;
        if (requiredTime > itaumax_) {
            itaumax_ = requiredTime;
            size_t totalSize = static_cast<size_t>(itaumax_) * static_cast<size_t>(ietamax_) * static_cast<size_t>(ixmax_);
            
            // Check for reasonable sizes
            if (totalSize > 1000000000) {  // 1 billion
                std::cerr << "[HydroProfile::loadPlaintextHydro] ERROR: Allocate size too large: " << totalSize << std::endl;
                continue;
            }
            
            hydrot_.resize(totalSize, 0.0);
            hydrox_.resize(totalSize, 0.0);
            hydroy_.resize(totalSize, 0.0);
        }

        // Store at (it, iy, ix) - no eta loop needed. Boost invariance is implicit.
        size_t idx = index(it, iy, ix, ietamax_, ixmax_);
        hydrot_[idx] = tdat;
        hydrox_[idx] = vxdat;
        hydroy_[idx] = vydat;
    }
    
    hydroTauMax_ = hydroTau0_ + hydroDtau_ * (itaumax_ - 1);
    
    std::cout << "Read plaintext hydro: itaumax=" << itaumax_ << " ixmax=" << ixmax_ 
              << " ietamax=" << ietamax_ << " lines_processed=" << line_count << " maxtemp=" << maxtemp << "\n";
    std::cout.flush();
    hydro.close();
}

double HydroProfile::getValue(const std::vector<double> &data, double tau, double x, double y) const {
    if (tau >= hydroTauMax_ || tau < hydroTau0_) {
        return 0.0;
    }

    int it = static_cast<int>(std::floor((tau - hydroTau0_) / hydroDtau_));
    it = clampIndex(it, itaumax_);
    double dt = fracFromIndex(tau, hydroTau0_, hydroDtau_, it);

    double xgrid = (hydroXmax_ + x) / hydroDx_;
    int ix = static_cast<int>(std::floor(xgrid));
    ix = clampIndex(ix, ixmax_);
    double dx = xgrid - static_cast<double>(ix);

    double ygrid = (hydroXmax_ + y) / hydroDx_;
    int iy = static_cast<int>(std::floor(ygrid));
    iy = clampIndex(iy, ietamax_);
    double dy = ygrid - static_cast<double>(iy);
    
    size_t base = index(it, iy, ix, ietamax_, ixmax_);
    size_t strideY = static_cast<size_t>(ixmax_);
    size_t strideT = static_cast<size_t>(ietamax_) * static_cast<size_t>(ixmax_);

    auto fetch = [&](int dt_i, int dy_i, int dx_i) {
        size_t idx = base + dt_i * strideT + dy_i * strideY + dx_i;
        if (idx >= data.size()) {
            std::cerr << "ERROR: fetch idx out of bounds! " << idx << " >= " << data.size() << std::endl;
            return 0.0;
        }
        return data[idx];
    };

    double value = 0.0;
    value += fetch(0, 0, 0) * (1. - dt) * (1. - dx) * (1. - dy);
    value += fetch(1, 0, 0) * dt * (1. - dx) * (1. - dy);
    value += fetch(0, 0, 1) * (1. - dt) * dx * (1. - dy);
    value += fetch(0, 1, 0) * (1. - dt) * (1. - dx) * dy;
    value += fetch(1, 0, 1) * dt * dx * (1. - dy);
    value += fetch(0, 1, 1) * (1. - dt) * dx * dy;
    value += fetch(1, 1, 0) * dt * (1. - dx) * dy;
    value += fetch(1, 1, 1) * dt * dx * dy;

    return value;
}

double HydroProfile::temperature(double tau, double x, double y) const {
    return getValue(hydrot_, tau, x, y) * tempScalingFactor_;
}

double HydroProfile::velocityX(double tau, double x, double y) const {
    return getValue(hydrox_, tau, x, y);
}

double HydroProfile::velocityY(double tau, double x, double y) const {
    return getValue(hydroy_, tau, x, y);
}
