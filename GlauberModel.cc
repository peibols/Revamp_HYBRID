#include "GlauberModel.h"

#include <fstream>
#include <iostream>
#include <cmath>
#include <cassert>
#include <sstream>

using std::vector;

GlauberModel::GlauberModel()
    :
      step_(0.),
      norm_(1.),
      bmin_(0.),
      bmax_(0.) {
    // Initialize arrays
    for (int i = 0; i < GLAUB_SIZE; i++) {
        for (int j = 0; j < GLAUB_SIZE; j++) {
            glaub_[i][j] = 0.;
        }
    }
    TA_.resize(TA_SIZE, 0.);
}

GlauberModel::~GlauberModel() = default;

void GlauberModel::readNuclear(int nhyd, const std::string &cent) {
    const char *glauFile = "./TAb2LL.dat";
    std::ifstream initial(glauFile);
    
    if (initial.fail()) {
        std::cerr << "Cannot open Glauber file: " << glauFile << std::endl;
        return;
    }
    
    std::cout << "Reading Initial Energy Density..." << std::endl;
    
    getCentralityBounds(cent, bmin_, bmax_);
    std::cout << "Bmin= " << bmin_ << " Bmax= " << bmax_ << std::endl;
    
    // Read nuclear overlap function
    double b2;
    for (unsigned int a = 0; a < TA_SIZE; a++) {
        initial >> b2 >> TA_[a];
        if (a == 1) step_ = b2;
    }
    
    initial.close();
}

void GlauberModel::getCentralityBounds(const std::string &cent, double &bmin, double &bmax) const {
    // LHC centrality boundaries
    if (cent == "0-5") {
        bmin = 0.;
        bmax = 3.5;
    } else if (cent == "5-10") {
        bmin = 3.5;
        bmax = 4.94;
    } else if (cent == "10-20") {
        bmin = 4.94;
        bmax = 6.98;
    } else if (cent == "20-30") {
        bmin = 6.98;
        bmax = 8.55;
    } else if (cent == "30-40") {
        bmin = 8.55;
        bmax = 9.88;
    } else if (cent == "40-50") {
        bmin = 9.88;
        bmax = 11.04;
    } else if (cent == "50-60") {
        bmin = 11.04;
        bmax = 12.09;
    } else if (cent == "60-70") {
        bmin = 12.09;
        bmax = 13.05;
    } else {
        std::cerr << "Unrecognized centrality: " << cent << std::endl;
        bmin = 0.;
        bmax = 20.;
    }
}

double GlauberModel::gTAA(double x, double y, double b) const {
    int il, irr;
    double rho2, use;
    
    rho2 = pow(x + b / 2., 2.) + y * y;
    il = static_cast<int>(rho2 / step_);
    rho2 = pow(x - b / 2., 2.) + y * y;
    irr = static_cast<int>(rho2 / step_);
    
    if (il >= static_cast<int>(TA_.size())) il = TA_.size() - 1;
    if (irr >= static_cast<int>(TA_.size())) irr = TA_.size() - 1;
    
    use = TA_[il] * TA_[irr] / TA_[0];
    return use;
}

void GlauberModel::sampleXY(double &x, double &y, numrand &nr) {
    double rho, phi;
    double P;
    double b;
    
    // Metropolis sampling
    bool accepted = false;
    while (!accepted) {
        b = sqrt((bmax_ * bmax_ - bmin_ * bmin_) * nr.rando() + bmin_ * bmin_);
        norm_ = 1.;
        norm_ = gTAA(0., 0., bmin_);
        
        rho = sqrt(150. * nr.rando());
        phi = 2. * 3.141592654 * nr.rando();
        x = rho * cos(phi);
        y = rho * sin(phi);
        
        P = nr.rando();
        if (P <= gTAA(x, y, b)) {
            accepted = true;
        }
    }
}

int GlauberModel::readNuclearIPSAT(int nhyd, const std::string &cent) {
    std::ostringstream filename;
    filename << "./NcollList.dat";
    std::ifstream initial(filename.str().c_str());
    if (initial.fail()) {
        std::cout << "Initial open fail " << filename.str() << std::endl;
        exit(1);
    }

    std::string s;
    // First line: some crap
    std::getline(initial, s);
    // Rest of lines: x, y;
    double x, y;
    binpos_.clear();
    while (std::getline(initial, s)) {
        std::istringstream iss(s);
        iss >> x >> y;
        binpos_.emplace_back(x, y);
    }
    return static_cast<int>(binpos_.size());
}

void GlauberModel::sampleXYIPSAT(double &x, double &y, int n, numrand &nr) {
    int index = static_cast<int>(n * nr.rando());
    if (index >= n) index = n - 1;
    x = binpos_[index].first;
    y = binpos_[index].second;
}
