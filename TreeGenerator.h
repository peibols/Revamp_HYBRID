#pragma once

#include <string>
#include <vector>

#include "Parton.h"

namespace Pythia8 {
class Pythia;
}

// Encapsulates Pythia initialization + event generation into a reusable object.
// This avoids global state and keeps the Pythia instance inside a class.
class TreeGenerator {
public:
    TreeGenerator();
    ~TreeGenerator();

    // Initialize Pythia with an explicit random seed and command file.
    void init(int seed, const std::string &cmndFile = "setup_pythia.cmnd");

    // Set trigger options (optional)
    void setTrigger(double pt, double eta, int id);

    // Generate one event and fill the output parton list.
    // Returns false if the event generation failed.
    // If trigger is set, checks for trigger particle.
    bool nextEvent(std::vector<Parton> &partons, double &weight, double &cross, double &cross_err);

private:
    Pythia8::Pythia *pythia_ = nullptr;
    bool use_trigger_ = false;
    double trigger_pt_ = 0.0;
    double trigger_eta_ = 0.0;
    int trigger_id_ = 0;
};
