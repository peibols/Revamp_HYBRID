#include "Config.h"
#include "HYBRID.h"

#include <iostream>

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <hybrid_input_file>\n";
        return 1;
    }

    Config cfg;
    if (!cfg.load(argv[1])) {
        std::cerr << "Failed to read configuration file: " << argv[1] << "\n";
        return 2;
    }

    HYBRID sim(cfg);
    sim.run();

    return 0;
}
