#pragma once

#include <algorithm>
#include <vector>

#include "Parton.h"
#include "Quench.h"
#include "ResponseDeposit.h"
#include "vector_operators.h"

// Build the positive medium-deposition ledger exactly once. Shower losses and
// propagated-recoiler losses are wake sources; holes remain separate negative
// thermal removals and are deliberately excluded.
inline std::vector<ResponseDeposit> build_response_ledger(
    const std::vector<Quench> &quenched,
    const std::vector<Parton> &partons,
    const std::vector<Quench> &recoiled) {
    std::vector<ResponseDeposit> deposits;
    const size_t shower_count = std::min(quenched.size(), partons.size());
    for (size_t i = 0; i < shower_count; ++i) {
        if (partons[i].GetD1() != -1) continue;
        deposits.push_back({partons[i].vGetP() - quenched[i].vGetP(),
                            static_cast<int>(i), ResponseDepositKind::ShowerLoss});
    }
    for (size_t i = 0; i < recoiled.size(); ++i) {
        const auto &response = recoiled[i];
        if (response.GetOrig() != "recoiler") continue;
        deposits.push_back({response.origEn() - response.vGetP(),
                            static_cast<int>(i), ResponseDepositKind::RecoilerLoss});
    }
    return deposits;
}
