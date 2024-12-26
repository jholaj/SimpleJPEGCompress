#pragma once

#include <vector>
#include <array>
#include "constants.hpp"

struct PixelBlock {
    std::vector<std::vector<std::array<int, 3>>> data;

    PixelBlock();
};
