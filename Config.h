#pragma once

#include <string>
#include <unordered_map>
#include <optional>

// A simple key/value parser for "key = value" text files.
// Lines beginning with '#' or '//' are ignored.
// Keys are case-sensitive.

struct Config {
    std::unordered_map<std::string, std::string> entries;

    bool load(const std::string &path);

    std::optional<std::string> getString(const std::string &key) const;
    std::optional<int> getInt(const std::string &key) const;
    std::optional<double> getDouble(const std::string &key) const;
    std::optional<bool> getBool(const std::string &key) const;

    // Convenience getters with default values
    std::string getStringOr(const std::string &key, const std::string &fallback) const;
    int getIntOr(const std::string &key, int fallback) const;
    double getDoubleOr(const std::string &key, double fallback) const;
    bool getBoolOr(const std::string &key, bool fallback) const;
};
