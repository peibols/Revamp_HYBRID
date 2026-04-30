#include "HydroProfile.h"

#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>

namespace {

int legacyAxisIndex(double coord, double delta, int maxCount) {
    if (coord >= 0.) {
        return static_cast<int>(coord / delta) + maxCount / 2;
    }
    return static_cast<int>(coord / delta) + maxCount / 2 - 1;
}

template <typename FetchFn>
double interpolateLegacyAveraged(double tau, double x, double y, double eta, FetchFn fetch) {
    constexpr int kMaxX = 100;
    constexpr int kMaxY = 100;
    constexpr int kMaxEta = 64;
    constexpr double kDeltaX = 0.3;
    constexpr double kDeltaY = 0.3;
    constexpr double kDeltaEta = 0.203125;
    constexpr double kDeltaT = 0.1;
    constexpr double kTau0 = 0.6;
    constexpr double kTau1 = 18.5;
    constexpr double kEta1 = 6.5;

    if (tau >= kTau1 || std::abs(eta) >= kEta1 || tau < kTau0) {
        return 0.0;
    }

    const int it = static_cast<int>((tau - kTau0) / kDeltaT);
    const double dt = (tau - kTau0 - static_cast<double>(it) * kDeltaT) / kDeltaT;

    const int ix = legacyAxisIndex(x, kDeltaX, kMaxX);
    const double dx = (x - static_cast<double>(ix - kMaxX / 2) * kDeltaX) / kDeltaX;

    const int iy = legacyAxisIndex(y, kDeltaY, kMaxY);
    const double dy = (y - static_cast<double>(iy - kMaxY / 2) * kDeltaY) / kDeltaY;

    int ieta;
    double deta;
    if (eta >= 0.) {
        ieta = static_cast<int>(eta / kDeltaEta) + kMaxEta / 2;
        deta = (eta - static_cast<double>(ieta - kMaxEta / 2) * kDeltaEta) / kDeltaEta;
    } else {
        ieta = static_cast<int>(eta / kDeltaEta) + kMaxEta / 2 - 1;
        deta = (eta - static_cast<double>(ieta - kMaxEta / 2) * kDeltaEta) / kDeltaEta;
    }

    if (ix < 0 || ix >= kMaxX || iy < 0 || iy >= kMaxY || ieta < 0 || ieta >= kMaxEta - 1) {
        return 0.0;
    }

    const double base =
        fetch(it, iy, ix) * (1. - dt) * (1. - dx) * (1. - dy) +
        fetch(it + 1, iy, ix) * dt * (1. - dx) * (1. - dy) +
        fetch(it, iy, ix + 1) * (1. - dt) * dx * (1. - dy) +
        fetch(it, iy + 1, ix) * (1. - dt) * (1. - dx) * dy +
        fetch(it + 1, iy, ix + 1) * dt * dx * (1. - dy) +
        fetch(it, iy + 1, ix + 1) * (1. - dt) * dx * dy +
        fetch(it + 1, iy + 1, ix) * dt * (1. - dx) * dy +
        fetch(it + 1, iy + 1, ix + 1) * dt * dx * dy;

    (void)deta;
    return base;
}

}

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
        tempScalingFactor_ = 0.197327;  // hbar*c [GeV·fm]: converts fm^-1 -> GeV
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

    // Pre-allocate storage for the first time slice so itau=0 writes are valid.
    itaumax_ = 1;
    size_t totalSize = static_cast<size_t>(itaumax_) * static_cast<size_t>(ietamax_) * static_cast<size_t>(ixmax_);
    hydrot_.assign(totalSize, 0.0);
    hydrox_.assign(totalSize, 0.0);
    hydroy_.assign(totalSize, 0.0);

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
        assert(idx < data.size());
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

double HydroProfile::temperatureElasticLegacy(double tau, double x, double y, double eta) const {
    if (mode_ != 0) return temperature(tau, x, y);
    if (itaumax_ < 2 || ixmax_ < 101 || ietamax_ < 101) return 0.0;
    return interpolateLegacyAveraged(
        tau, x, y, eta,
        [&](int it, int iy, int ix) {
            if (it < 0 || it + 1 >= itaumax_ || iy < 0 || iy + 1 >= ietamax_ || ix < 0 || ix + 1 >= ixmax_) {
                return 0.0;
            }
            return hydrot_[index(it, iy, ix, ietamax_, ixmax_)];
        }
    ) * tempScalingFactor_;
}

double HydroProfile::velocityXElasticLegacy(double tau, double x, double y, double eta) const {
    if (mode_ != 0) return velocityX(tau, x, y);
    if (itaumax_ < 2 || ixmax_ < 101 || ietamax_ < 101) return 0.0;
    return interpolateLegacyAveraged(
        tau, x, y, eta,
        [&](int it, int iy, int ix) {
            if (it < 0 || it + 1 >= itaumax_ || iy < 0 || iy + 1 >= ietamax_ || ix < 0 || ix + 1 >= ixmax_) {
                return 0.0;
            }
            return hydrox_[index(it, iy, ix, ietamax_, ixmax_)];
        }
    );
}

double HydroProfile::velocityYElasticLegacy(double tau, double x, double y, double eta) const {
    if (mode_ != 0) return velocityY(tau, x, y);
    if (itaumax_ < 2 || ixmax_ < 101 || ietamax_ < 101) return 0.0;
    return interpolateLegacyAveraged(
        tau, x, y, eta,
        [&](int it, int iy, int ix) {
            if (it < 0 || it + 1 >= itaumax_ || iy < 0 || iy + 1 >= ietamax_ || ix < 0 || ix + 1 >= ixmax_) {
                return 0.0;
            }
            return hydroy_[index(it, iy, ix, ietamax_, ixmax_)];
        }
    );
}

void HydroProfile::getValues(double tau, double x, double y, double &temp, double &vx, double &vy) const {
    if (tau >= hydroTauMax_ || tau < hydroTau0_) {
        temp = vx = vy = 0.0;
        return;
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

    // Compute the 8 interpolation weights once
    double w000 = (1.-dt)*(1.-dx)*(1.-dy);
    double w100 = dt    *(1.-dx)*(1.-dy);
    double w010 = (1.-dt)*dx    *(1.-dy);
    double w001 = (1.-dt)*(1.-dx)*dy;
    double w110 = dt    *dx    *(1.-dy);
    double w011 = (1.-dt)*dx    *dy;
    double w101 = dt    *(1.-dx)*dy;
    double w111 = dt    *dx    *dy;

    size_t i000 = base;
    size_t i100 = base + strideT;
    size_t i010 = base + 1;
    size_t i001 = base + strideY;
    size_t i110 = base + strideT + 1;
    size_t i011 = base + strideY + 1;
    size_t i101 = base + strideT + strideY;
    size_t i111 = base + strideT + strideY + 1;

    assert(i111 < hydrot_.size());

    temp = (hydrot_[i000]*w000 + hydrot_[i100]*w100 + hydrot_[i010]*w010 + hydrot_[i001]*w001 +
            hydrot_[i110]*w110 + hydrot_[i011]*w011 + hydrot_[i101]*w101 + hydrot_[i111]*w111) * tempScalingFactor_;
    vx   =  hydrox_[i000]*w000 + hydrox_[i100]*w100 + hydrox_[i010]*w010 + hydrox_[i001]*w001 +
            hydrox_[i110]*w110 + hydrox_[i011]*w011 + hydrox_[i101]*w101 + hydrox_[i111]*w111;
    vy   =  hydroy_[i000]*w000 + hydroy_[i100]*w100 + hydroy_[i010]*w010 + hydroy_[i001]*w001 +
            hydroy_[i110]*w110 + hydroy_[i011]*w011 + hydroy_[i101]*w101 + hydroy_[i111]*w111;
}
