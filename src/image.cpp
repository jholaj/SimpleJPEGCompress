#include "image.hpp"
#include "dct.hpp"
#include "exceptions.hpp"
#include <fstream>
#include <cmath>

Image::Image(const std::string& filename) {
    DCTTransform::initializeMatrices();
    loadImage(filename);
}

void Image::loadImage(const std::string& filename) {
    std::ifstream file(filename, std::ios::binary);
    if (!file) {
        throw ImageException("Cannot find input file: " + filename);
    }

    // Read BMP header
    std::vector<char> header(54);
    if (!file.read(header.data(), 54)) {
        throw ImageException("Failed to read BMP header");
    }

    // Verify BMP signature
    if (header[0] != 'B' || header[1] != 'M') {
        throw ImageException("Invalid BMP file format");
    }

    // Extract image dimensions
    width = *reinterpret_cast<int*>(&header[18]);
    height = *reinterpret_cast<int*>(&header[22]);

    if (width <= 0 || height <= 0) {
        throw ImageException("Invalid image dimensions");
    }

    // Calculate padding
    int padding = (4 - (width * 3) % 4) % 4;

    // Allocate memory for pixels
    pixels.resize(height, std::vector<std::array<int, 3>>(width));

    // Read pixel data
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            std::array<unsigned char, 3> colors;
            if (!file.read(reinterpret_cast<char*>(colors.data()), 3)) {
                throw ImageException("Failed to read pixel data");
            }
            pixels[y][x] = {colors[2], colors[1], colors[0]}; // BGR to RGB
        }
        file.seekg(padding, std::ios::cur);
    }
}

void Image::saveImage(const std::string& filename,
                     const std::vector<std::vector<std::array<int, 3>>>& imageData) {
    std::ofstream file(filename, std::ios::binary);
    if (!file) {
        throw ImageException("Cannot create output file: " + filename);
    }

    // Calculate padding and file size
    int padding = (4 - (width * 3) % 4) % 4;
    int fileSize = 54 + (3 * width + padding) * height;

    // Prepare BMP header
    std::vector<char> header(54, 0);
    // Signature
    header[0] = 'B';
    header[1] = 'M';
    // File size
    *reinterpret_cast<int*>(&header[2]) = fileSize;
    // Data offset
    *reinterpret_cast<int*>(&header[10]) = 54;
    // Header size
    *reinterpret_cast<int*>(&header[14]) = 40;
    // Width and height
    *reinterpret_cast<int*>(&header[18]) = width;
    *reinterpret_cast<int*>(&header[22]) = height;
    // Planes
    *reinterpret_cast<short*>(&header[26]) = 1;
    // Bits per pixel
    *reinterpret_cast<short*>(&header[28]) = 24;

    // Write header
    file.write(header.data(), 54);

    // Write pixel data
    std::vector<char> paddingBytes(padding, 0);
    for (int y = height - 1; y >= 0; --y) {
        for (int x = 0; x < width; ++x) {
            // Write BGR values
            unsigned char bgr[] = {
                static_cast<unsigned char>(imageData[y][x][2]), // B
                static_cast<unsigned char>(imageData[y][x][1]), // G
                static_cast<unsigned char>(imageData[y][x][0])  // R
            };
            file.write(reinterpret_cast<char*>(bgr), 3);
        }
        if (padding > 0) {
            file.write(paddingBytes.data(), padding);
        }
    }
}

std::vector<std::vector<std::array<int, 3>>> Image::convertToYCbCr() {
    std::vector<std::vector<std::array<int, 3>>> ycbcrPixels(height,
                                                            std::vector<std::array<int, 3>>(width));

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            double R = pixels[y][x][0];
            double G = pixels[y][x][1];
            double B = pixels[y][x][2];

            int Y  = std::round(0.299 * R + 0.587 * G + 0.114 * B);
            int Cb = std::round(128 - 0.168736 * R - 0.331264 * G + 0.5 * B);
            int Cr = std::round(128 + 0.5 * R - 0.418688 * G - 0.081312 * B);

            ycbcrPixels[y][x] = {
                static_cast<int>(clamp(Y, 0, 255)),
                static_cast<int>(clamp(Cb, 0, 255)),
                static_cast<int>(clamp(Cr, 0, 255))
            };
        }
    }

    return ycbcrPixels;
}

std::vector<PixelBlock> Image::divideIntoBlocks(
    const std::vector<std::vector<std::array<int, 3>>>& imageData) {
    std::vector<PixelBlock> blocks;
    int blockCountX = (width + JPEGConstants::BLOCK_SIZE - 1) / JPEGConstants::BLOCK_SIZE;
    int blockCountY = (height + JPEGConstants::BLOCK_SIZE - 1) / JPEGConstants::BLOCK_SIZE;

    blocks.reserve(blockCountX * blockCountY);

    for (int blockY = 0; blockY < blockCountY; ++blockY) {
        for (int blockX = 0; blockX < blockCountX; ++blockX) {
            PixelBlock block;

            for (int y = 0; y < JPEGConstants::BLOCK_SIZE; ++y) {
                for (int x = 0; x < JPEGConstants::BLOCK_SIZE; ++x) {
                    int imageX = blockX * JPEGConstants::BLOCK_SIZE + x;
                    int imageY = blockY * JPEGConstants::BLOCK_SIZE + y;

                    if (imageX < width && imageY < height) {
                        block.data[y][x] = imageData[imageY][imageX];
                    } else {
                        block.data[y][x] = {0, 128, 128}; // Default YCbCr values
                    }
                }
            }

            blocks.push_back(block);
        }
    }

    return blocks;
}

void Image::compressImage(const std::string& outputFilename) {
    auto ycbcrImage = convertToYCbCr();
    auto blocks = divideIntoBlocks(ycbcrImage);

    for (auto& block : blocks) {
        DCTTransform::applyDCT(block);
        quantize(block);
    }

    auto huffmanCodes = huffmanEncode(blocks);

    for (auto& block : blocks) {
        dequantize(block);
        DCTTransform::applyIDCT(block);
    }

    auto rgbImage = convertToRGB(blocks);
    saveImage(outputFilename, rgbImage);
}

std::vector<std::vector<std::array<int, 3>>> Image::convertToRGB(const std::vector<PixelBlock>& blocks) {
    std::vector<std::vector<std::array<int, 3>>> rgbImage(height, std::vector<std::array<int, 3>>(width));

    int blockIndex = 0;
    for (int blockY = 0; blockY < height; blockY += JPEGConstants::BLOCK_SIZE) {
        for (int blockX = 0; blockX < width; blockX += JPEGConstants::BLOCK_SIZE) {
            const auto& block = blocks[blockIndex++];

            for (int y = 0; y < JPEGConstants::BLOCK_SIZE && (blockY + y) < height; ++y) {
                for (int x = 0; x < JPEGConstants::BLOCK_SIZE && (blockX + x) < width; ++x) {
                    int imageY = blockY + y;
                    int imageX = blockX + x;

                    double Y = block.data[y][x][0];
                    double Cb = block.data[y][x][1] - 128.0;
                    double Cr = block.data[y][x][2] - 128.0;

                    int R = std::round(Y + 1.402 * Cr);
                    int G = std::round(Y - 0.344136 * Cb - 0.714136 * Cr);
                    int B = std::round(Y + 1.772 * Cb);

                    rgbImage[imageY][imageX] = {
                        static_cast<int>(clamp(R, 0, 255)),
                        static_cast<int>(clamp(G, 0, 255)),
                        static_cast<int>(clamp(B, 0, 255))
                    };
                }
            }
        }
    }

    return rgbImage;
}

void Image::quantize(PixelBlock& block) {
    for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
        for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
            // Luminance quantization
            block.data[i][j][0] = std::round(block.data[i][j][0] /
                (JPEGConstants::LUMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));

            // Chrominance quantization
            for (int c = 1; c < 3; ++c) {
                block.data[i][j][c] = std::round(block.data[i][j][c] /
                    (JPEGConstants::CHROMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));
            }
        }
    }
}

void Image::dequantize(PixelBlock& block) {
    for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
        for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
            // Luminance dequantization
            block.data[i][j][0] = std::round(block.data[i][j][0] *
                (JPEGConstants::LUMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));

            // Chrominance dequantization
            for (int c = 1; c < 3; ++c) {
                block.data[i][j][c] = std::round(block.data[i][j][c] *
                    (JPEGConstants::CHROMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));
            }
        }
    }
}

std::map<int, std::string> Image::huffmanEncode(const std::vector<PixelBlock>& blocks) {
    std::vector<int> frequencies(JPEGConstants::COLOR_RANGE, 0);
    std::vector<int> values(JPEGConstants::COLOR_RANGE);

    for (int i = 0; i < JPEGConstants::COLOR_RANGE; ++i) {
        values[i] = i;
    }

    // Calculate frequencies
    for (const auto& block : blocks) {
        for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
            for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                for (int c = 0; c < 3; ++c) {
                    int value = clamp(block.data[i][j][c], 0, 255);
                    frequencies[value]++;
                }
            }
        }
    }

    HuffmanTree huffmanTree;
    return huffmanTree.encode(values, frequencies);
}

int Image::clamp(int value, int min, int max) {
    return value < min ? min : (value > max ? max : value);
}
