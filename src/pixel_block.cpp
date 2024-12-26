#include "pixel_block.hpp"
#include "constants.hpp"

PixelBlock::PixelBlock()
    : data(JPEGConstants::BLOCK_SIZE, std::vector<std::array<int, 3>>(JPEGConstants::BLOCK_SIZE)) {}
