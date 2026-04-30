#pragma once

#include <vector>
#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"

// Encapsulates Pythia-based hadronization (Lund string model).
// Replaces global Pythia instance with a class-based approach.
class LundGenerator {
public:
    LundGenerator();
    ~LundGenerator();

    // Initialize Pythia for hadronization with an explicit random seed.
    void init(int seed);

    // Hadronize vacuum partons
    void hadronizeVacuum(const std::vector<Parton> &partons, 
                         std::vector<Hadron> &vhadrons);

    // Hadronize medium-modified partons
    bool hadronizeMedium(const std::vector<Quench> &quenched, 
                         std::vector<Hadron> &qhadrons,
                         int hadro_type);

private:
    class Pythia8Impl;
    Pythia8Impl *pimpl_;  // Pointer to implementation (Pythia instance hidden)

    // Helper functions
    void processVacuumPartons(const std::vector<Parton> &partons, 
                              std::vector<Hadron> &vhadrons);
    bool processQuenchedPartons(const std::vector<Quench> &quenched, 
                                std::vector<Hadron> &qhadrons,
                                int hadro_type);
};
