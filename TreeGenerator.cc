#include "TreeGenerator.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

#include <iostream>
#include <sstream>

using namespace Pythia8;

TreeGenerator::TreeGenerator() = default;
TreeGenerator::~TreeGenerator() {
    delete pythia_;
}

void TreeGenerator::init(int njob, const std::string &cmndFile) {
    if (pythia_) {
        delete pythia_;
        pythia_ = nullptr;
    }

    pythia_ = new Pythia();

    std::ostringstream pythiaset;
    pythiaset << cmndFile;
    pythia_->readFile(pythiaset.str());

    int seed = 33 + njob;
    std::ostringstream seedstring;
    seedstring << "Random:seed = " << seed;
    pythia_->readString(seedstring.str().c_str());

    pythia_->init();
}

void TreeGenerator::setTrigger(double pt, double eta, int id) {
    use_trigger_ = true;
    trigger_pt_ = pt;
    trigger_eta_ = eta;
    trigger_id_ = id;
}

bool TreeGenerator::nextEvent(std::vector<Parton> &partons, double &weight, double &cross, double &cross_err) {
    if (!pythia_) {
        std::cerr << "TreeGenerator not initialized (call init() before generating events).\n";
        return false;
    }

    const Info &info = pythia_->info;

    if (!pythia_->next()) {
        return false;
    }

    weight = info.weight();
    cross = info.sigmaGen();
    cross_err = info.sigmaErr();

    // Check for trigger
    bool have_trigger = false;
    if (use_trigger_) {
        for (int i = 0; i < pythia_->event.size(); ++i) {
            const Particle &p = pythia_->event[i];
            if (!p.isFinal()) continue;
            if (p.id() != trigger_id_) continue;
            if (p.pT() < trigger_pt_) continue;
            if (std::abs(p.eta()) > trigger_eta_) continue;
            have_trigger = true;
            break;
        }
        if (!have_trigger) {
            // If trigger is set but not found, skip this event
            return false;
        }
    }

    partons.clear();

    for (int i = 0; i < pythia_->event.size(); ++i) {
        if (!pythia_->event[i].isFinal()) continue;

        std::vector<double> p;
        for (unsigned int j = 1; j < 4; ++j) p.push_back(pythia_->event[i].p()[j]);
        p.push_back(pythia_->event[i].p()[0]);

        // Simply store remnants (virt 0, mother 0, daughters -1)
        if (pythia_->event[i].status() == 63) {
            partons.emplace_back(p, 0., pythia_->event[i].m(), 0, -1, -1, pythia_->event[i].id(), "rem", pythia_->event[i].col(), pythia_->event[i].acol(), true);
            continue;
        }

        int use = i;
        int m1 = 0;
        int m2 = 0;
        do {
            m1 = pythia_->event[use].mother1();
            m2 = pythia_->event[use].mother2();
            if (m1 == m2) use = m1;
        } while (m1 == m2);

        double virt = sqrt(std::abs(pow(pythia_->event[i].e(), 2.) - pow(pythia_->event[i].px(), 2.) - pow(pythia_->event[i].py(), 2.) - pow(pythia_->event[i].pz(), 2.) - pythia_->event[i].m2()));

        partons.emplace_back(p, virt, pythia_->event[i].m(), m1, -1, -1, pythia_->event[i].id(), "ps", pythia_->event[i].col(), pythia_->event[i].acol(), false);

        if (pythia_->event[m1].status() == -41) {
            partons.back().SetMom(-1);
            partons.back().SetOrig("isr");
            partons.back().SetIsDone(true);
        }

        if (pythia_->event[m1].status() == -21) {
            partons.back().SetMom(-1);
            partons.back().SetOrig("hs");
            partons.back().SetIsDone(true);
        }
    }

    // Reconstruct tree, using momentum of daughters to get the one for mothers
    int changes = 0;
    do {
        changes = 1;
        const unsigned int ps = partons.size();
        for (unsigned int i = 0; i < ps; ++i) {
            if (!partons[i].GetIsDone()) {
                for (unsigned int j = 0; j < ps; ++j) {
                    if (partons[i].GetMom() == partons[j].GetMom() && i != j && !partons[j].GetIsDone()) {
                        int mom = partons[i].GetMom();

                        std::vector<double> p = partons[i].vGetP();
                        const std::vector<double> &p_j = partons[j].vGetP();
                        for (size_t k = 0; k < p.size() && k < p_j.size(); ++k) {
                            p[k] += p_j[k];
                        }
                        double virt = sqrt(std::abs(pow(p[3], 2.) - pow(p[0], 2.) - pow(p[1], 2.) - pow(p[2], 2.) - pow(pythia_->event[mom].m(), 2.)));

                        int use = mom;
                        int m1 = 0;
                        int m2 = 0;
                        do {
                            m1 = pythia_->event[use].mother1();
                            m2 = pythia_->event[use].mother2();
                            if (m1 == m2) use = m1;
                        } while (m1 == m2);

                        partons.emplace_back(p, virt, pythia_->event[mom].m(), m1, i, j, pythia_->event[mom].id(), "ps", pythia_->event[mom].col(), pythia_->event[mom].acol(), false);

                        partons[i].SetMom(partons.size() - 1);
                        partons[j].SetMom(partons.size() - 1);
                        partons[i].SetIsDone(true);
                        partons[j].SetIsDone(true);

                        if (pythia_->event[m1].status() == -41) {
                            partons.back().SetMom(-1);
                            partons.back().SetOrig("isr");
                            partons.back().SetIsDone(true);
                        }

                        if (pythia_->event[m1].status() == -21) {
                            partons.back().SetMom(-1);
                            partons.back().SetOrig("hs");
                            partons.back().SetIsDone(true);
                        }

                        changes = 0;
                        break;
                    }
                }
            }
        }
    } while (changes == 0);

    return true;
}