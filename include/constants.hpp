#pragma once

#include <vector>

namespace JPEGConstants {
    constexpr int BLOCK_SIZE = 8;
    constexpr int COLOR_RANGE = 256;
    constexpr double QUALITY_FACTOR = 1.0;

    extern const std::vector<std::vector<int>> LUMINANCE_QUANT_MATRIX;
    extern const std::vector<std::vector<int>> CHROMINANCE_QUANT_MATRIX;
}
