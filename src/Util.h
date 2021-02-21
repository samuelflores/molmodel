#ifndef UTIL_H
#define UTIL_H

#include <string>

inline
std::string trim_both(const std::string &s) {
    const auto first = s.find_first_not_of(" ");
    const auto last = s.find_last_not_of(" ");

    if (first == std::string::npos)
        return "";
    return s.substr(first, last - first + 1);
}

#endif // UTIL_H
