#include "Config.h"

#include <fstream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <cctype>
#include <stdexcept>

static inline std::string trim(std::string s) {
    auto is_space = [](unsigned char c) { return std::isspace(c); };
    s.erase(s.begin(), std::find_if_not(s.begin(), s.end(), is_space));
    s.erase(std::find_if_not(s.rbegin(), s.rend(), is_space).base(), s.end());
    return s;
}

bool Config::load(const std::string &path) {
    std::ifstream f(path);
    if (!f) return false;

    entries.clear();

    std::string line;
    int line_number = 0;
    while (std::getline(f, line)) {
        ++line_number;
        line = trim(line);
        if (line.empty() || line.rfind("#", 0) == 0 || line.rfind("//", 0) == 0) continue;

        auto eq = line.find('=');
        if (eq == std::string::npos) {
            throw std::invalid_argument("Malformed config line " + std::to_string(line_number) +
                                        ": expected key = value");
        }

        auto key = trim(line.substr(0, eq));
        auto value = trim(line.substr(eq + 1));
        if (key.empty() || value.empty()) {
            throw std::invalid_argument("Empty config key or value on line " +
                                        std::to_string(line_number));
        }
        if (!entries.emplace(key, value).second) {
            throw std::invalid_argument("Duplicate config key: " + key);
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
            size_t pos = 0;
            const int parsed = std::stoi(*v, &pos);
            if (pos != v->size()) return std::nullopt;
            return parsed;
        } catch (...) {
            return std::nullopt;
        }
    }
    return std::nullopt;
}

std::optional<double> Config::getDouble(const std::string &key) const {
    if (auto v = getString(key)) {
        try {
            size_t pos = 0;
            const double parsed = std::stod(*v, &pos);
            if (pos != v->size() || !std::isfinite(parsed)) return std::nullopt;
            return parsed;
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
    if (entries.find(key) == entries.end()) return fallback;
    if (auto v = getInt(key)) return *v;
    throw std::invalid_argument("Invalid integer for config key " + key + ": " + entries.at(key));
}

double Config::getDoubleOr(const std::string &key, double fallback) const {
    if (entries.find(key) == entries.end()) return fallback;
    if (auto v = getDouble(key)) return *v;
    throw std::invalid_argument("Invalid number for config key " + key + ": " + entries.at(key));
}

bool Config::getBoolOr(const std::string &key, bool fallback) const {
    if (entries.find(key) == entries.end()) return fallback;
    if (auto v = getBool(key)) return *v;
    throw std::invalid_argument("Invalid boolean for config key " + key + ": " + entries.at(key));
}

void Config::validateKnownKeys(const std::unordered_set<std::string> &known) const {
    for (const auto &entry : entries) {
        if (known.find(entry.first) == known.end()) {
            throw std::invalid_argument("Unknown config key: " + entry.first);
        }
    }
}
