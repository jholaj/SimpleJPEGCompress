#include "dct.hpp"
#include "image.hpp"
#include <cmath>
#include <algorithm>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

std::vector<std::vector<double>> DCTTransform::dctMatrix;
std::vector<std::vector<double>> DCTTransform::dctMatrixTranspose;

void DCTTransform::initializeMatrices() {
    if (!dctMatrix.empty()) return;  // Early return if already initialized

    dctMatrix = createDCTMatrix();
    dctMatrixTranspose = transpose(dctMatrix);
}

std::vector<std::vector<double>> DCTTransform::createDCTMatrix() {
    const int size = JPEGConstants::BLOCK_SIZE;
    std::vector<std::vector<double>> matrix(size, std::vector<double>(size));
    const double normalizationFactor = std::sqrt(2.0 / size);

    // Pre-calculate common values
    const double piDivSize = M_PI / (2.0 * size);

    for (int i = 0; i < size; ++i) {
        const double alpha = (i == 0) ? normalizationFactor / std::sqrt(2.0) : normalizationFactor;
        for (int j = 0; j < size; ++j) {
            matrix[i][j] = alpha * std::cos((2.0 * j + 1.0) * i * piDivSize);
        }
    }

    return matrix;
}

std::vector<std::vector<double>> DCTTransform::transpose(const std::vector<std::vector<double>>& matrix) {
    const size_t rows = matrix.size();
    const size_t cols = matrix[0].size();
    std::vector<std::vector<double>> result(cols, std::vector<double>(rows));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            result[j][i] = matrix[i][j];
        }
    }
    return result;
}

// Helper function for matrix multiplication
static void multiplyMatrices(const std::vector<std::vector<double>>& a,
                           const std::vector<std::vector<double>>& b,
                           std::vector<std::vector<double>>& result) {
    const int size = JPEGConstants::BLOCK_SIZE;
    for (int i = 0; i < size; ++i) {
        for (int j = 0; j < size; ++j) {
            double sum = 0.0;
            for (int k = 0; k < size; ++k) {
                sum += a[i][k] * b[k][j];
            }
            result[i][j] = sum;
        }
    }
}

void DCTTransform::applyDCT(PixelBlock& block) {
    const int size = JPEGConstants::BLOCK_SIZE;
    std::vector<std::vector<double>> channelData(size, std::vector<double>(size));
    std::vector<std::vector<double>> temp(size, std::vector<double>(size));
    std::vector<std::vector<double>> result(size, std::vector<double>(size));

    for (int channel = 0; channel < 3; ++channel) {
        // Extract channel data and apply Y channel offset
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                channelData[i][j] = block.data[i][j][channel] - (channel == 0 ? 128 : 0);
            }
        }

        // Perform DCT: result = dctMatrix * channelData * dctMatrixTranspose
        multiplyMatrices(dctMatrix, channelData, temp);
        multiplyMatrices(temp, dctMatrixTranspose, result);

        // Store results
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                block.data[i][j][channel] = std::round(result[i][j]);
            }
        }
    }
}

void DCTTransform::applyIDCT(PixelBlock& block) {
    const int size = JPEGConstants::BLOCK_SIZE;
    std::vector<std::vector<double>> channelData(size, std::vector<double>(size));
    std::vector<std::vector<double>> temp(size, std::vector<double>(size));
    std::vector<std::vector<double>> result(size, std::vector<double>(size));

    for (int channel = 0; channel < 3; ++channel) {
        // Extract channel data
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                channelData[i][j] = block.data[i][j][channel];
            }
        }

        // Perform IDCT: result = dctMatrixTranspose * channelData * dctMatrix
        multiplyMatrices(dctMatrixTranspose, channelData, temp);
        multiplyMatrices(temp, dctMatrix, result);

        // Store results and apply Y channel offset
        for (int i = 0; i < size; ++i) {
            for (int j = 0; j < size; ++j) {
                double value = result[i][j] + (channel == 0 ? 128 : 0);
                block.data[i][j][channel] = Image::clamp(std::round(value), 0, 255);
            }
        }
    }
}
