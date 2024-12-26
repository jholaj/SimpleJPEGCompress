#pragma once

#include <string>
#include <vector>
#include <array>
#include "pixel_block.hpp"
#include "huffman.hpp"

class Image {
public:
    explicit Image(const std::string& filename);
    static int clamp(int value, int min, int max);
    void compressImage(const std::string& outputFilename);

private:
    int width, height;
    std::vector<std::vector<std::array<int, 3>>> pixels;

    void loadImage(const std::string& filename);
    void saveImage(const std::string& filename,
                  const std::vector<std::vector<std::array<int, 3>>>& imageData);

    std::vector<std::vector<std::array<int, 3>>> convertToYCbCr();
    std::vector<std::vector<std::array<int, 3>>> convertToRGB(const std::vector<PixelBlock>& blocks);
    std::vector<PixelBlock> divideIntoBlocks(const std::vector<std::vector<std::array<int, 3>>>& imageData);

    void quantize(PixelBlock& block);
    void dequantize(PixelBlock& block);
    std::map<int, std::string> huffmanEncode(const std::vector<PixelBlock>& blocks);

};
