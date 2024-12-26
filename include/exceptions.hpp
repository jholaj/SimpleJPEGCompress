#pragma once

#include <stdexcept>
#include <string>

class ImageException : public std::runtime_error {
public:
    explicit ImageException(const std::string& message);
};
