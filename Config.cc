#include "Config.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cctype>

static inline std::string trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c); };
    s.erase(s.begin(), std::find_if_not(s.begin(), s.end(), is_space));
    s.erase(std::find_if_not(s.rbegin(), s.rend(), is_space).base(), s.end());
    return s;
}

bool Config::load(const std::string &path) {
    std::ifstream f(path);
    if (!f) return false;

    std::string line;
    while (std::getline(f, line)) {
        line = trim(line);
        if (line.empty() || line.rfind("#", 0) == 0 || line.rfind("//", 0) == 0) continue;

        auto eq = line.find('=');
        if (eq == std::string::npos) continue;

        auto key = trim(line.substr(0, eq));
        auto value = trim(line.substr(eq + 1));
        if (!key.empty()) {
            entries[key] = value;
        }
    }

    return true;
}

std::optional<std::string> Config::getString(const std::string &key) const {
    auto it = entries.find(key);
    if (it == entries.end()) return std::nullopt;
    return it->second;
}

std::optional<int> Config::getInt(const std::string &key) const {
    if (auto v = getString(key)) {
        try {
            return std::stoi(*v);
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<double> Config::getDouble(const std::string &key) const {
    if (auto v = getString(key)) {
        try {
            return std::stod(*v);
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<bool> Config::getBool(const std::string &key) const {
    if (auto v = getString(key)) {
        const auto &s = *v;
        if (s == "1" || s == "true" || s == "True" || s == "TRUE") return true;
        if (s == "0" || s == "false" || s == "False" || s == "FALSE") return false;
    }
    return std::nullopt;
}

std::string Config::getStringOr(const std::string &key, const std::string &fallback) const {
    if (auto v = getString(key)) return *v;
    return fallback;
}

int Config::getIntOr(const std::string &key, int fallback) const {
    if (auto v = getInt(key)) return *v;
    return fallback;
}

double Config::getDoubleOr(const std::string &key, double fallback) const {
    if (auto v = getDouble(key)) return *v;
    return fallback;
}

bool Config::getBoolOr(const std::string &key, bool fallback) const {
    if (auto v = getBool(key)) return *v;
    return fallback;
}
