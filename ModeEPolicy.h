#pragma once

namespace modee {

// A coherent Moliere source must itself carry color. In particular, a failed
// daughter probe below gamma -> q qbar cannot be reinterpreted as a photon
// scattering. The colored daughters may still generate resolving probes.
constexpr bool is_colored_coherent_source(int pdg_id) noexcept {
    return pdg_id == 21 ||
           (pdg_id >= 1 && pdg_id <= 6) ||
           (pdg_id <= -1 && pdg_id >= -6);
}

}  // namespace modee
