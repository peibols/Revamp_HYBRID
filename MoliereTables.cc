#include "MoliereTables.h"

#include <iostream>
#include <stdexcept>

#include "read_tables.hpp"

bool MoliereTables::loaded_ = false;
std::string MoliereTables::loaded_path_;

void MoliereTables::ensureLoaded(const std::string &tables_path) {
    if (tables_path.empty()) {
        throw std::runtime_error("Elastic scattering requested but tables_path is empty");
    }
    if (loaded_) {
        if (loaded_path_ != tables_path) {
            throw std::runtime_error("Moliere tables already loaded from a different path");
        }
        return;
    }

    std::cout << "Loading Moliere tables from " << tables_path << std::endl;
    read_tables(tables_path);
    use_tables = true;
    loaded_path_ = tables_path;
    loaded_ = true;
}

bool MoliereTables::loaded() {
    return loaded_;
}
