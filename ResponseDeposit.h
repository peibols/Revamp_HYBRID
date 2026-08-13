#pragma once

#include <array>

enum class ResponseDepositKind {
    ShowerLoss,
    RecoilerLoss
};

// Positive four-momentum transferred from a propagated hard object to the
// medium. Holes are separate negative thermal removals and never enter here.
struct ResponseDeposit {
    std::array<double,4> momentum = {0., 0., 0., 0.};
    int source_index = -1;
    ResponseDepositKind kind = ResponseDepositKind::ShowerLoss;
};
