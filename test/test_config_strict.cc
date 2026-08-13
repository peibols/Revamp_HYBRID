#include "Config.h"

#include <cassert>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char **argv) {
    assert(argc == 2);
    const std::string path = argv[1];
    {
        std::ofstream out(path);
        out << "Nev = 10events\n";
    }
    Config bad_number;
    assert(bad_number.load(path));
    bool rejected = false;
    try {
        (void)bad_number.getIntOr("Nev", 1);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    assert(rejected);

    {
        std::ofstream out(path);
        out << "do_lres = maybe\n";
    }
    Config bad_bool;
    assert(bad_bool.load(path));
    rejected = false;
    try {
        (void)bad_bool.getBoolOr("do_lres", false);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    assert(rejected);

    {
        std::ofstream out(path);
        out << "kappa = nan\n";
    }
    Config bad_finite;
    assert(bad_finite.load(path));
    rejected = false;
    try {
        (void)bad_finite.getDoubleOr("kappa", 1.);
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    assert(rejected);

    {
        std::ofstream out(path);
        out << "Nev = 4\n";
    }
    Config reloaded;
    assert(reloaded.load(path));
    {
        std::ofstream out(path);
        out << "do_lres = true\n";
    }
    assert(reloaded.load(path));
    assert(!reloaded.getInt("Nev"));
    assert(reloaded.getBoolOr("do_lres", false));

    {
        std::ofstream out(path);
        out << "typo_key = 1\n";
    }
    Config unknown;
    assert(unknown.load(path));
    rejected = false;
    try {
        unknown.validateKnownKeys({"Nev"});
    } catch (const std::invalid_argument &) {
        rejected = true;
    }
    assert(rejected);

    std::cout << "strict config unit tests passed\n";
    return 0;
}
