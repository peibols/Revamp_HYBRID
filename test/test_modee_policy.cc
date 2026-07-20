#include "ModeEPolicy.h"

#include <cassert>
#include <initializer_list>

int main() {
    for (int flavor = 1; flavor <= 6; ++flavor) {
        assert(modee::is_colored_coherent_source(flavor));
        assert(modee::is_colored_coherent_source(-flavor));
    }
    assert(modee::is_colored_coherent_source(21));

    for (int neutral_id : {0, 11, -11, 22, 23, 111, 211}) {
        assert(!modee::is_colored_coherent_source(neutral_id));
    }
    return 0;
}
