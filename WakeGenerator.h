#pragma once

#include <array>
#include <vector>
#include "Parton.h"
#include "Quench.h"
#include "Wake.h"
#include "Random.h"

// Encapsulates hadron wake/back-reaction generation from quenched partons.
// Replaces global parameters in HadWake.cc with class-based configuration.
class WakeGenerator {
public:
    WakeGenerator();
    ~WakeGenerator();

    // Generate wake particles from energy loss.
    void generate(const std::vector<Quench> &quenched, 
                  const std::vector<Parton> &partons, 
                  std::vector<Wake> &wake, 
                  numrand &nr);

    // Configuration setters (all parameters that were global)
    void setTransCut(double val) { transcut_ = val; }
    void setBaseSigma(double val) { basesig_ = val; }
    void setNrun(int val) { Nrun_ = val; }
    void setMaxPtSq(double val) { maxptsq_ = val; }
    void setMaxRapidity(double val) { maxrap_ = val; }
    void setTolerance(double val) { tole_ = val; }
    void setMassPion(double val) { masspi_ = val; }
    void setMassProton(double val) { masspro_ = val; }

    // Getters
    double getTransCut() const { return transcut_; }
    double getBaseSigma() const { return basesig_; }
    int getNrun() const { return Nrun_; }
    double getMaxPtSq() const { return maxptsq_; }
    double getMaxRapidity() const { return maxrap_; }
    double getTolerance() const { return tole_; }
    int getTooMuch() const { return toomuch_; }

private:
    // Parameters (previously global)
    double transcut_;      // Lower threshold for transverse mass
    double basesig_;       // Starting gaussian width for Metropolis
    int Nrun_;             // Maximum iterations of Metropolis
    double maxptsq_;       // Maximum squared pt for MC
    double maxrap_;        // Maximum absolute rapidity for MC
    double tole_;          // Non-conservation tolerance in GeV
    double masspi_;        // Pion mass
    double masspro_;       // Proton mass

    // Arrays (previously global)
    double masstra_[2];
    double normcoop_[2];
    double maxcooper_[2];
    int toomuch_;

    // Private methods
    double rapid(double pt, double pz) const;
    int set_charge(int spe, numrand &nr) const;
    double thermal(int spe, double ptrand) const;
    void one_body(std::vector<Wake> &wake, const std::array<double,4>& delta, const std::array<double,4>& momback, 
                  double ptlost, double mtlost, double raplost, numrand &nr, int spe, int mode);
    bool tryAddResidualWake(std::vector<Wake> &wake,
                            const std::array<double,4>& delta,
                            std::array<double,4>& momback,
                            numrand &nr) const;
    std::array<double,4> vec_abs(const std::array<double,4>& p) const;
};
