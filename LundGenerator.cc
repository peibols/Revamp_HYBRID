#include "LundGenerator.h"

#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"

#include <iostream>
#include <cmath>
#include <cstring>
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

bool LundGenerator::hadronizeMedium(const std::vector<Quench> &quenched,
                                    std::vector<Hadron> &qhadrons,
                                    int hadro_type) {
    if (!pimpl_ || !pimpl_->pythia) return false;

    Pythia8::Pythia *pythia = pimpl_->pythia;
    Pythia8::ParticleData &pdt = pythia->particleData;

    pythia->event.reset();

    auto append_final_hadrons = [&]() {
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
    };

    if (hadro_type == 0) {
        int colsum = 0;
        for (const auto &q : quenched) {
            if (q.GetD1() != -1) continue;

            int ide = q.GetId();
            auto p = q.vGetP();
            double px = p[0];
            double py = p[1];
            double pz = p[2];
            double mm = pdt.m0(ide);
            int col = q.GetCol();
            int acol = q.GetAcol();
            double ee = sqrt(px * px + py * py + pz * pz + mm * mm);

            pythia->event.append(ide, 23, col, acol, px, py, pz, ee, mm);
            colsum += col - acol;
        }

        if (colsum != 0) {
            std::cout << " Sumadecolor= " << colsum << std::endl;
        }
        if (!pythia->next()) {
            std::cout << " Event generation aborted prematurely, owing to error" << std::endl;
            return false;
        }

        append_final_hadrons();
        return true;
    }

    if (hadro_type != 1) {
        std::cout << " Unsupported hadro_type= " << hadro_type << std::endl;
        return false;
    }

    double rempx = 0.2;
    double rempy = 0.2;
    double p_fake = 2500.;
    double rempz = p_fake;
    double reme = std::sqrt(rempx * rempx + rempy * rempy + rempz * rempz);

    std::vector<Quench> pIn;
    std::vector<Quench> uncolored;
    for (const auto &q : quenched) {
        if (q.GetD1() == -1 && q.vGetP()[3] != 0. && q.GetOrig() != "rem") {
            if (q.GetId() == 21 || std::abs(q.GetId()) <= 6) pIn.push_back(q);
            else uncolored.push_back(q);
        }
    }

    std::vector<int> col(pIn.size() + 2, 0);
    std::vector<int> acol(pIn.size() + 2, 0);
    std::vector<int> isdone(pIn.size() + 2, 0);

    int nquarks = 0;
    std::vector<int> isquark(pIn.size() + 2, 0);
    for (unsigned int ipart = 0; ipart < pIn.size(); ++ipart) {
        if (std::abs(pIn[ipart].GetId()) <= 6) {
            isquark[nquarks] = int(ipart);
            nquarks += 1;
        }
    }

    int nstrings = std::max(int(double(nquarks) / 2. + 0.6), 1);
    int istring = 0;
    std::vector<int> one_end(nstrings, -1);
    std::vector<int> two_end(nstrings, -1);
    if (nquarks == 0) {
        vector<double> pfirstq = {rempx, rempy, rempz, reme};
        pIn.emplace_back(Parton(pfirstq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true));
        isquark[nquarks] = int(pIn.size()) - 1;
        nquarks += 1;
        isdone[pIn.size() - 1] = 1;
        one_end[0] = int(pIn.size()) - 1;

        vector<double> psecq = {rempx, rempy, -rempz, reme};
        pIn.emplace_back(Parton(psecq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true));
        isquark[nquarks] = int(pIn.size()) - 1;
        nquarks += 1;
        isdone[pIn.size() - 1] = 1;
        two_end[istring] = int(pIn.size()) - 1;
    }

    for (int iquark = 0; iquark < nquarks; ++iquark) {
        if (isdone[isquark[iquark]] != 0) continue;

        isdone[isquark[iquark]] = 1;
        one_end[istring] = isquark[iquark];
        double min_delR = 1000000.;
        int partner = -2;
        for (int jquark = 0; jquark < nquarks; ++jquark) {
            if (iquark == jquark) continue;
            int d_jquark = isquark[jquark];
            if (isdone[d_jquark] == 0) {
                double delR = pIn[isquark[iquark]].delta_R(pIn[d_jquark]);
                if (delR < min_delR) {
                    min_delR = delR;
                    partner = jquark;
                }
            }
        }
        if (partner != -2) {
            isdone[isquark[partner]] = 1;
            two_end[istring] = isquark[partner];
            istring += 1;
        } else {
            vector<double> prem = {rempx, rempy, rempz, reme};
            pIn.emplace_back(Parton(prem, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true));
            isquark[nquarks] = int(pIn.size()) - 1;
            nquarks += 1;
            isdone[pIn.size() - 1] = 1;
            two_end[istring] = int(pIn.size()) - 1;
            std::cout << "Attached quark remnant flying down +Pz beam" << std::endl;
        }
    }

    std::vector<int> my_string(pIn.size(), 0);
    for (unsigned int ipart = 0; ipart < pIn.size(); ++ipart) {
        if (pIn[ipart].GetId() != 21) continue;
        double min_delR = 100000.;
        for (int ns = 0; ns < nstrings; ++ns) {
            int fq = one_end[ns];
            int sq = two_end[ns];
            double f_delR = pIn[ipart].delta_R(pIn[fq]);
            double s_delR = pIn[ipart].delta_R(pIn[sq]);
            double delR = (f_delR + s_delR) / 2.;
            if (delR < min_delR) {
                min_delR = delR;
                my_string[ipart] = ns;
            }
        }
    }

    auto append_parton = [&](const Quench &q, int c, int ac) {
        int ide = q.GetId();
        auto p = q.vGetP();
        double px = p[0];
        double py = p[1];
        double pz = p[2];
        double mm = pdt.m0(int(ide));
        double ee = std::sqrt(px * px + py * py + pz * pz + mm * mm);
        if (c == 0 && ac == 0 && (ide == 21 || std::abs(ide) <= 6)) {
            std::cout << "Stopping because of colorless parton trying to be introduced in PYTHIA string" << std::endl;
            return false;
        }
        pythia->event.append(int(ide), 23, c, ac, px, py, pz, ee, mm);
        return true;
    };

    int lab_col = 102;
    for (int ns = 0; ns < nstrings; ++ns) {
        int tquark = one_end[ns];
        if (pIn[tquark].GetId() > 0) col[tquark] = lab_col;
        else acol[tquark] = lab_col;
        lab_col += 1;
        int link = tquark;

        if (!append_parton(pIn[tquark], col[tquark], acol[tquark])) return false;

        int changes = 1;
        do {
            changes = 0;
            double min_delR = 100000.;
            int next_link = 0;
            for (unsigned int ipart = 0; ipart < pIn.size(); ++ipart) {
                if (pIn[ipart].GetId() == 21 && isdone[ipart] == 0 && my_string[ipart] == ns) {
                    changes = 1;
                    double delR = pIn[link].delta_R(pIn[ipart]);
                    if (delR < min_delR) {
                        min_delR = delR;
                        next_link = int(ipart);
                    }
                }
            }
            if (changes == 1) {
                isdone[next_link] = 1;
                if (col[link] == lab_col - 1) {
                    col[next_link] = lab_col;
                    acol[next_link] = lab_col - 1;
                } else {
                    col[next_link] = lab_col - 1;
                    acol[next_link] = lab_col;
                }
                lab_col += 1;
                link = next_link;
                if (!append_parton(pIn[next_link], col[next_link], acol[next_link])) return false;
            }
        } while (changes == 1);

        if (col[link] == lab_col - 1) {
            col[two_end[ns]] = 0;
            acol[two_end[ns]] = lab_col - 1;
        } else {
            col[two_end[ns]] = lab_col - 1;
            acol[two_end[ns]] = 0;
        }

        if (col[two_end[ns]] != 0) {
            if (pIn[two_end[ns]].GetId() < 0) pIn[two_end[ns]].SetId(-pIn[two_end[ns]].GetId());
        } else {
            if (pIn[two_end[ns]].GetId() > 0) pIn[two_end[ns]].SetId(-pIn[two_end[ns]].GetId());
        }

        if (!append_parton(pIn[two_end[ns]], col[two_end[ns]], acol[two_end[ns]])) return false;
    }

    for (const auto &q : uncolored) {
        int ide = q.GetId();
        auto p = q.vGetP();
        double px = p[0];
        double py = p[1];
        double pz = p[2];
        double mm = pdt.m0(int(ide));
        double ee = std::sqrt(px * px + py * py + pz * pz + mm * mm);
        pythia->event.append(int(ide), 23, 0, 0, px, py, pz, ee, mm);
    }

    if (!pythia->next()) {
        std::cout << " GOT ISSUES " << std::endl;
        return false;
    }

    append_final_hadrons();
    return true;
}

void LundGenerator::processVacuumPartons(const std::vector<Parton> &partons, 
                                         std::vector<Hadron> &vhadrons) {
    hadronizeVacuum(partons, vhadrons);
}

bool LundGenerator::processQuenchedPartons(const std::vector<Quench> &quenched, 
                                           std::vector<Hadron> &qhadrons,
                                           int hadro_type) {
    return hadronizeMedium(quenched, qhadrons, hadro_type);
}
