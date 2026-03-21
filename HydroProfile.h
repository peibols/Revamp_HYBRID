#pragma once

#include <string>
#include <vector>

// Encapsulates hydro-field data and interpolation. Designed to replace global arrays
// and free functions with a single owning object.
class HydroProfile {
public:
    HydroProfile();
    ~HydroProfile();

    // Load hydro evolution data. Supports both event-by-event (ebe_hydro=1) and 
    // event-averaged (ebe_hydro=0) modes.
    // Mode 0: reads "./hydroinfoPlaintxtHuichaoFormat.dat" (plaintext, boost-invariant)
    // Mode 1: reads "evolution_all_xyeta.dat" (binary, event-by-event)
    void loadHydro(int mode, const std::string &cent);

    // Legacy interface: directs to loadHydro(1, cent) for backward compatibility
    void loadIpsat(int nhyd, const std::string &cent, const std::string &filename = "evolution_all_xyeta.dat") {
        loadHydro(1, cent);
    }

    // Evaluate fields at (tau,x,y).
    double temperature(double tau, double x, double y) const;
    double velocityX(double tau, double x, double y) const;
    double velocityY(double tau, double x, double y) const;

private:
    int mode_ = 1;  // 0 = event-averaged (plaintext), 1 = event-by-event (IPSAT binary)
    double tempScalingFactor_ = 1.0;  // Temperature scaling: 0.2 for mode 0, 1.0 for mode 1

    int ixmax_ = 0;
    int ietamax_ = 0;
    int itaumax_ = 0;

    double hydroDx_ = 0;
    double hydroXmax_ = 0;
    double hydroTau0_ = 0;
    double hydroDtau_ = 0;
    double hydroTauMax_ = 0;

    std::vector<double> hydrot_;
    std::vector<double> hydrox_;
    std::vector<double> hydroy_;

    // Mode-specific loaders
    void loadPlaintextHydro(const std::string &cent);
    void loadIpsatBinary(const std::string &cent, const std::string &filename);

    static size_t index(int it, int iy, int ix, int ietamax, int ixmax);
    static int clampIndex(int idx, int maxCount);
    static double fracFromIndex(double coord, double origin, double spacing, int idx);

    double getValue(const std::vector<double> &data, double tau, double x, double y) const;
};
