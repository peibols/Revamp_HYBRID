#include "ResponseLedger.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <vector>

namespace {
bool close(double a, double b) {
    return std::abs(a - b) < 1.e-12;
}
}

int main() {
    Parton shower(std::array<double,4>{10., 0., 0., 10.}, 0., 0., -1, -1, -1,
                  1, "ps", 0, 0, true);
    Quench quenched(shower);
    quenched.vSetP(std::array<double,4>{8., 0., 0., 8.});

    Parton recoil_part(std::array<double,4>{2., 0., 0., 2.}, 0., 0., -1, -1, -1,
                       2, "recoiler", 0, 0, true);
    Quench recoil(recoil_part);
    recoil.setOrigEn({3., 0., 0., 3.});

    Parton hole_part(std::array<double,4>{-1., 0., 0., 1.}, 0., 0., -1, -1, -1,
                     1, "hole", 0, 0, true);
    Quench hole(hole_part);
    hole.setOrigEn({4., 0., 0., 4.});

    const auto ledger = build_response_ledger(
        std::vector<Quench>{quenched}, std::vector<Parton>{shower},
        std::vector<Quench>{recoil, hole});
    assert(ledger.size() == 2);
    assert(ledger[0].kind == ResponseDepositKind::ShowerLoss);
    assert(ledger[1].kind == ResponseDepositKind::RecoilerLoss);
    assert(close(ledger[0].momentum[3], 2.));
    assert(close(ledger[1].momentum[3], 1.));

    std::cout << "response ledger unit tests passed\n";
    return 0;
}
