#pragma once

#include <string>
#include <vector>
#include "Random.h"

// Encapsulates Glauber model for nuclear geometry and impact parameter sampling.
// Replaces global arrays and parameters in Glauber_PbPb.cc with class-based approach.
class GlauberModel {
public:
    GlauberModel();
    ~GlauberModel();

    // Load nuclear overlap function (TA) from file
    void readNuclear(int nhyd, const std::string &cent);

    // Sample impact parameter and position in transverse plane
    void sampleXY(double &x, double &y, numrand &nr);

    // Get centrality bounds
    void getCentralityBounds(const std::string &cent, double &bmin, double &bmax) const;

    // IPSAT methods (from Glauber_Chun.cc)
    int readNuclearIPSAT(int nhyd, const std::string &cent);
    void sampleXYIPSAT(double &x, double &y, int n, numrand &nr);

private:
    // Global arrays/parameters (previously static/global in Glauber_PbPb.cc)
    static const int GLAUB_SIZE = 200;
    static const int TA_SIZE = 4000;

    double glaub_[GLAUB_SIZE][GLAUB_SIZE];  // Initial energy density grid

    std::vector<double> TA_;  // Nuclear overlap function
    double step_;
    double norm_;

    double bmin_;  // Impact parameter bounds for current centrality
    double bmax_;

    // IPSAT data (from Glauber_Chun.cc)
    std::vector<std::pair<double, double>> binpos_;

    // Helper functions (previously global free functions)
    double gTAA(double x, double y, double b) const;
};
