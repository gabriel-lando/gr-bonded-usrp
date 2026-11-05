//
// Copyright 2026 Gabriel Lando
// SPDX-License-Identifier: GPL-3.0-or-later
//
// Comma-separated parameter parsing shared by the block constructors.
//

#pragma once

#include <cctype>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

namespace bonded { namespace detail {

//! Split on commas, trim whitespace around each item, drop empty items.
inline std::vector<std::string> split_csv(const std::string& csv)
{
    std::vector<std::string> out;
    size_t pos = 0;
    while (pos <= csv.size()) {
        size_t end = csv.find(',', pos);
        if (end == std::string::npos) {
            end = csv.size();
        }
        size_t a = pos, b = end;
        while (a < b && std::isspace(static_cast<unsigned char>(csv[a]))) {
            a++;
        }
        while (b > a && std::isspace(static_cast<unsigned char>(csv[b - 1]))) {
            b--;
        }
        if (b > a) {
            out.push_back(csv.substr(a, b - a));
        }
        pos = end + 1;
    }
    return out;
}

//! Parse a list of numbers; throws std::invalid_argument naming `what` on a
//! malformed item (a typo must not silently drop a per-radio value).
inline std::vector<double> parse_csv_doubles(const std::string& csv, const std::string& what)
{
    std::vector<double> out;
    for (const auto& item : split_csv(csv)) {
        size_t used = 0;
        double v    = 0.0;
        try {
            v = std::stod(item, &used);
        } catch (const std::exception&) {
            used = 0;
        }
        if (used != item.size()) {
            throw std::invalid_argument(what + ": not a number: '" + item + "'");
        }
        out.push_back(v);
    }
    return out;
}

}} // namespace bonded::detail
