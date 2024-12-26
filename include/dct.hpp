#pragma once

#include <vector>
#include "pixel_block.hpp"

class DCTTransform {
public:
    static void initializeMatrices();
    static void applyDCT(PixelBlock& block);
    static void applyIDCT(PixelBlock& block);

private:
    static std::vector<std::vector<double>> dctMatrix;
    static std::vector<std::vector<double>> dctMatrixTranspose;

    static std::vector<std::vector<double>> createDCTMatrix();
    static std::vector<std::vector<double>> transpose(const std::vector<std::vector<double>>& matrix);
};
