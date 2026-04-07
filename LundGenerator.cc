#include "LundGenerator.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

#include <iostream>
#include <cmath>
#include <string>

using namespace Pythia8;
using std::vector;

// Pimpl: private implementation to hide Pythia dependency
class LundGenerator::Pythia8Impl {
public:
    Pythia8Impl() : pythia(new Pythia8::Pythia()) {}
    ~Pythia8Impl() { delete pythia; }
    Pythia8::Pythia *pythia;
};

LundGenerator::LundGenerator() : pimpl_(new Pythia8Impl()) {
    // Don't call init() here; it will be called explicitly in HYBRID::run()
}

LundGenerator::~LundGenerator() {
    delete pimpl_;
}

void LundGenerator::init(int seed) {
    if (!pimpl_ || !pimpl_->pythia) return;
    
    Pythia8::Pythia *pythia = pimpl_->pythia;
    
    // Hadronization-only mode: skip hard process and parton shower,
    // run string fragmentation on manually inserted partons.
    pythia->readString("Random:setSeed = on");
    pythia->readString("Random:seed = " + std::to_string(seed));
    pythia->readString("ProcessLevel:all = off");
    pythia->init();
}

void LundGenerator::hadronizeVacuum(const std::vector<Parton> &partons, 
                                    std::vector<Hadron> &vhadrons) {
    if (!pimpl_ || !pimpl_->pythia) return;
    
    Pythia8::Pythia *pythia = pimpl_->pythia;
    Pythia8::ParticleData &pdt = pythia->particleData;

    pythia->event.reset();  // clear event record before each call

    int colsum = 0;
    for (unsigned int i = 0; i < partons.size(); i++) {
        if (partons[i].GetD1() != -1) continue;
        
        int ide = partons[i].GetId();
        auto p = partons[i].vGetP();
        double px = p[0];
        double py = p[1];
        double pz = p[2];
        double mm = pdt.m0(ide);
        int col = partons[i].GetCol();
        int acol = partons[i].GetAcol();
        double ee = sqrt(px * px + py * py + pz * pz + mm * mm);
        
        // Insert into Pythia event record
        pythia->event.append(ide, 23, col, acol, px, py, pz, ee, mm);
        
        if (col != 0) colsum += col;
        if (acol != 0) colsum -= acol;
    }
    
    // Generate hadronization
    if (colsum == 0 && pythia->event.size() > 1) {
        pythia->next();
        
        for (int i = 0; i < pythia->event.size(); ++i) {
            if (!pythia->event[i].isFinal()) continue;
            
            vector<double> hp;
            hp.push_back(pythia->event[i].px());
            hp.push_back(pythia->event[i].py());
            hp.push_back(pythia->event[i].pz());
            hp.push_back(pythia->event[i].e());
            
            vhadrons.emplace_back(
                Parton(hp, 0., pythia->event[i].m(), 0, -1, -1, 
                       pythia->event[i].id(), "lund", 0, 0, true),
                pythia->event[i].charge(),
                -1
            );
        }
    }
}

void LundGenerator::hadronizeMedium(const std::vector<Quench> &quenched, 
                                    std::vector<Hadron> &qhadrons) {
    if (!pimpl_ || !pimpl_->pythia) return;
    
    Pythia8::Pythia *pythia = pimpl_->pythia;
    Pythia8::ParticleData &pdt = pythia->particleData;

    pythia->event.reset();  // clear event record before each call

    int colsum = 0;
    for (unsigned int i = 0; i < quenched.size(); i++) {
        if (quenched[i].GetD1() != -1) continue;
        
        int ide = quenched[i].GetId();
        auto p = quenched[i].vGetP();
        
        double px = p[0];
        double py = p[1];
        double pz = p[2];
        double mm = pdt.m0(ide);
        int col = quenched[i].GetCol();
        int acol = quenched[i].GetAcol();
        double ee = sqrt(px * px + py * py + pz * pz + mm * mm);
        
        pythia->event.append(ide, 23, col, acol, px, py, pz, ee, mm);
        
        if (col != 0) colsum += col;
        if (acol != 0) colsum -= acol;
    }
    
    if (colsum == 0 && pythia->event.size() > 1) {
        pythia->next();
        
        for (int i = 0; i < pythia->event.size(); ++i) {
            if (!pythia->event[i].isFinal()) continue;
            
            vector<double> hp;
            hp.push_back(pythia->event[i].px());
            hp.push_back(pythia->event[i].py());
            hp.push_back(pythia->event[i].pz());
            hp.push_back(pythia->event[i].e());
            
            qhadrons.emplace_back(
                Parton(hp, 0., pythia->event[i].m(), 0, -1, -1, 
                       pythia->event[i].id(), "lund", 0, 0, true),
                pythia->event[i].charge(),
                -1
            );
        }
    }
}

void LundGenerator::processVacuumPartons(const std::vector<Parton> &partons, 
                                         std::vector<Hadron> &vhadrons) {
    hadronizeVacuum(partons, vhadrons);
}

void LundGenerator::processQuenchedPartons(const std::vector<Quench> &quenched, 
                                           std::vector<Hadron> &qhadrons) {
    hadronizeMedium(quenched, qhadrons);
}
