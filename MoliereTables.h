#pragma once

#include <string>

class MoliereTables {
public:
    static void ensureLoaded(const std::string &tables_path);
    static bool loaded();

private:
    static bool loaded_;
    static std::string loaded_path_;
};
