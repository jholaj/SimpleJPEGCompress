#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <queue>
#include <map>
#include <memory>
#include <stdexcept>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

using namespace std;

// Forward declarations
class Image;
class HuffmanTree;

// Custom exceptions
class ImageException : public std::runtime_error {
public:
    explicit ImageException(const string& message) : std::runtime_error(message) {}
};

// Constants for JPEG compression
namespace JPEGConstants {
    constexpr int BLOCK_SIZE = 8;
    constexpr int COLOR_RANGE = 256;
    constexpr double QUALITY_FACTOR = 1.0; // (0.1 - 2.0)

    // Standard JPEG quantization matrices
    const vector<vector<int>> LUMINANCE_QUANT_MATRIX = {
        {16, 11, 10, 16, 24, 40, 51, 61},
        {12, 12, 14, 19, 26, 58, 60, 55},
        {14, 13, 16, 24, 40, 57, 69, 56},
        {14, 17, 22, 29, 51, 87, 80, 62},
        {18, 22, 37, 56, 68, 109, 103, 77},
        {24, 35, 55, 64, 81, 104, 113, 92},
        {49, 64, 78, 87, 103, 121, 120, 101},
        {72, 92, 95, 98, 112, 100, 103, 99}
    };

    const vector<vector<int>> CHROMINANCE_QUANT_MATRIX = {
        {17, 18, 24, 47, 99, 99, 99, 99},
        {18, 21, 26, 66, 99, 99, 99, 99},
        {24, 26, 56, 99, 99, 99, 99, 99},
        {47, 66, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99},
        {99, 99, 99, 99, 99, 99, 99, 99}
    };
}

// Struktura pro uchování bloku pixelů
struct PixelBlock {
    vector<vector<array<int, 3>>> data;

    PixelBlock() : data(JPEGConstants::BLOCK_SIZE,
                       vector<array<int, 3>>(JPEGConstants::BLOCK_SIZE)) {}
};

struct Node {
    int data;
    unsigned freq;
    unique_ptr<Node> left, right;

    Node(int data, unsigned freq) : data(data), freq(freq), left(nullptr), right(nullptr) {}
};

struct CompareNodes {
    bool operator()(const unique_ptr<Node>& l, const unique_ptr<Node>& r) const {
        return l->freq > r->freq;
    }
};

class HuffmanTree {
private:
    unique_ptr<Node> root;

    void generateCodes(Node* node, string code, map<int, string>& huffmanCodes) {
        if (!node) return;

        if (node->data != -1) {
            huffmanCodes[node->data] = code;
        }

        generateCodes(node->left.get(), code + "0", huffmanCodes);
        generateCodes(node->right.get(), code + "1", huffmanCodes);
    }

public:
    map<int, string> encode(const vector<int>& data, const vector<int>& freq) {
        priority_queue<unique_ptr<Node>, vector<unique_ptr<Node>>, CompareNodes> minHeap;

        // Create leaf nodes
        for (size_t i = 0; i < data.size(); ++i) {
            if (freq[i] > 0) {
                minHeap.push(make_unique<Node>(data[i], freq[i]));
            }
        }

        // Build Huffman tree
        while (minHeap.size() > 1) {
            unique_ptr<Node> left = std::move(const_cast<unique_ptr<Node>&>(minHeap.top()));
            minHeap.pop();
            unique_ptr<Node> right = std::move(const_cast<unique_ptr<Node>&>(minHeap.top()));
            minHeap.pop();

            auto parent = make_unique<Node>(-1, left->freq + right->freq);
            parent->left = move(left);
            parent->right = move(right);
            minHeap.push(move(parent));
        }

        root = std::move(const_cast<unique_ptr<Node>&>(minHeap.top()));
        map<int, string> huffmanCodes;
        generateCodes(root.get(), "", huffmanCodes);
        return huffmanCodes;
    }

    int decode(const string& encodedStr, const map<int, string>& huffmanCodes) {
        Node* current = root.get();
        for (char bit : encodedStr) {
            current = (bit == '0') ? current->left.get() : current->right.get();

            if (!current) {
                throw ImageException("Invalid Huffman code sequence");
            }

            if (!current->left && !current->right) {
                return current->data;
            }
        }
        throw ImageException("Incomplete Huffman code sequence");
    }
};

class Image {
private:
    int width, height;
    vector<vector<array<int, 3>>> pixels;
    static vector<vector<double>> dctMatrix;
    static vector<vector<double>> dctMatrixTranspose;

    // Precomputed DCT matrix for optimalization
    static void initializeDCTMatrices() {
        if (dctMatrix.empty()) {
            dctMatrix = createDCTMatrix();
            dctMatrixTranspose = transpose(dctMatrix);
        }
    }

    static vector<vector<double>> createDCTMatrix() {
        vector<vector<double>> matrix(JPEGConstants::BLOCK_SIZE,
                                    vector<double>(JPEGConstants::BLOCK_SIZE));

        for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
            for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                double alpha = (i == 0) ? 1.0 / sqrt(2.0) : 1.0;
                matrix[i][j] = alpha * cos((2.0 * j + 1.0) * i * M_PI / (2.0 * JPEGConstants::BLOCK_SIZE));
            }
        }

        double normalizationFactor = sqrt(2.0 / JPEGConstants::BLOCK_SIZE);
        for (auto& row : matrix) {
            for (double& val : row) {
                val *= normalizationFactor;
            }
        }

        return matrix;
    }

    static vector<vector<double>> transpose(const vector<vector<double>>& matrix) {
        vector<vector<double>> result(matrix[0].size(), vector<double>(matrix.size()));
        for (size_t i = 0; i < matrix.size(); ++i) {
            for (size_t j = 0; j < matrix[0].size(); ++j) {
                result[j][i] = matrix[i][j];
            }
        }
        return result;
    }

public:
    explicit Image(const string& filename) {
        initializeDCTMatrices();
        loadImage(filename);
    }

    void loadImage(const string& filename) {
        ifstream file(filename, ios::binary);
        if (!file) {
            throw ImageException("Cannot open input file: " + filename);
        }

        // Read BMP header
        vector<char> header(54);
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
        pixels.resize(height, vector<array<int, 3>>(width));

        // Read pixel data
        for (int y = height - 1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                array<unsigned char, 3> colors;
                if (!file.read(reinterpret_cast<char*>(colors.data()), 3)) {
                    throw ImageException("Failed to read pixel data");
                }
                pixels[y][x] = {colors[2], colors[1], colors[0]}; // BGR to RGB
            }
            file.seekg(padding, ios::cur);
        }
    }

    void saveImage(const string& filename, const vector<vector<array<int, 3>>>& imageData) {
        ofstream file(filename, ios::binary);
        if (!file) {
            throw ImageException("Cannot create output file: " + filename);
        }

        // Calculate padding and file size
        int padding = (4 - (width * 3) % 4) % 4;
        int fileSize = 54 + (3 * width + padding) * height;

        // Prepare BMP header
        vector<char> header(54, 0);
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
        vector<char> paddingBytes(padding, 0);
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

    vector<vector<array<int, 3>>> convertToYCbCr() {
        vector<vector<array<int, 3>>> ycbcrPixels(height, vector<array<int, 3>>(width));

        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                double R = pixels[y][x][0];
                double G = pixels[y][x][1];
                double B = pixels[y][x][2];

                // Přesnější konverzní koeficienty
                int Y  = round(0.299 * R + 0.587 * G + 0.114 * B);
                int Cb = round(128 - 0.168736 * R - 0.331264 * G + 0.5 * B);
                int Cr = round(128 + 0.5 * R - 0.418688 * G - 0.081312 * B);

                // Kontrola rozsahu hodnot
                Y = clamp(Y, 0, 255);
                Cb = clamp(Cb, 0, 255);
                Cr = clamp(Cr, 0, 255);

                ycbcrPixels[y][x] = {Y, Cb, Cr};
            }
        }

        return ycbcrPixels;
    }

    vector<PixelBlock> divideIntoBlocks(const vector<vector<array<int, 3>>>& imageData) {
        vector<PixelBlock> blocks;
        int blockCountX = (width + JPEGConstants::BLOCK_SIZE - 1) / JPEGConstants::BLOCK_SIZE;
        int blockCountY = (height + JPEGConstants::BLOCK_SIZE - 1) / JPEGConstants::BLOCK_SIZE;

        blocks.reserve(blockCountX * blockCountY);

        for (int blockY = 0; blockY < blockCountY; ++blockY) {
            for (int blockX = 0; blockX < blockCountX; ++blockX) {
                PixelBlock block;

                // Process each pixel in the block
                for (int y = 0; y < JPEGConstants::BLOCK_SIZE; ++y) {
                    for (int x = 0; x < JPEGConstants::BLOCK_SIZE; ++x) {
                        int imageX = blockX * JPEGConstants::BLOCK_SIZE + x;
                        int imageY = blockY * JPEGConstants::BLOCK_SIZE + y;

                        // Handle image boundaries
                        if (imageX < width && imageY < height) {
                            block.data[y][x] = imageData[imageY][imageX];
                        } else {
                            // Padding for incomplete blocks
                            block.data[y][x] = {0, 128, 128}; // Default YCbCr values
                        }
                    }
                }

                blocks.push_back(block);
            }
        }

        return blocks;
    }

    void applyDCT(PixelBlock& block) {
            // For each channel... (Y, Cb, Cr)
            for (int channel = 0; channel < 3; ++channel) {
                vector<vector<double>> channelData(JPEGConstants::BLOCK_SIZE,
                                                 vector<double>(JPEGConstants::BLOCK_SIZE));
                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        channelData[i][j] = block.data[i][j][channel] - (channel == 0 ? 128 : 0);
                    }
                }
                vector<vector<double>> temp(JPEGConstants::BLOCK_SIZE,
                                          vector<double>(JPEGConstants::BLOCK_SIZE, 0.0));
                vector<vector<double>> result(JPEGConstants::BLOCK_SIZE,
                                            vector<double>(JPEGConstants::BLOCK_SIZE, 0.0));

                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        for (int k = 0; k < JPEGConstants::BLOCK_SIZE; ++k) {
                            temp[i][j] += dctMatrix[i][k] * channelData[k][j];
                        }
                    }
                }

                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        for (int k = 0; k < JPEGConstants::BLOCK_SIZE; ++k) {
                            result[i][j] += temp[i][k] * dctMatrixTranspose[k][j];
                        }
                    }
                }
                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        block.data[i][j][channel] = round(result[i][j]);
                    }
                }
            }
        }

        void quantize(PixelBlock& block) {
            for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                    block.data[i][j][0] = round(block.data[i][j][0] /
                        (JPEGConstants::LUMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));

                    for (int c = 1; c < 3; ++c) {
                        block.data[i][j][c] = round(block.data[i][j][c] /
                            (JPEGConstants::CHROMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));
                    }
                }
            }
        }

        void dequantize(PixelBlock& block) {
            for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                    block.data[i][j][0] = round(block.data[i][j][0] *
                        (JPEGConstants::LUMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));

                    for (int c = 1; c < 3; ++c) {
                        block.data[i][j][c] = round(block.data[i][j][c] *
                            (JPEGConstants::CHROMINANCE_QUANT_MATRIX[i][j] * JPEGConstants::QUALITY_FACTOR));
                    }
                }
            }
        }

        void applyIDCT(PixelBlock& block) {
            for (int channel = 0; channel < 3; ++channel) {
                vector<vector<double>> channelData(JPEGConstants::BLOCK_SIZE,
                                                 vector<double>(JPEGConstants::BLOCK_SIZE));

                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        channelData[i][j] = block.data[i][j][channel];
                    }
                }

                vector<vector<double>> temp(JPEGConstants::BLOCK_SIZE,
                                          vector<double>(JPEGConstants::BLOCK_SIZE, 0.0));
                vector<vector<double>> result(JPEGConstants::BLOCK_SIZE,
                                            vector<double>(JPEGConstants::BLOCK_SIZE, 0.0));

                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        for (int k = 0; k < JPEGConstants::BLOCK_SIZE; ++k) {
                            temp[i][j] += dctMatrixTranspose[i][k] * channelData[k][j];
                        }
                    }
                }

                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        for (int k = 0; k < JPEGConstants::BLOCK_SIZE; ++k) {
                            result[i][j] += temp[i][k] * dctMatrix[k][j];
                        }
                    }
                }

                // Uložení výsledku zpět do bloku a přičtení 128 pro Y kanál
                for (int i = 0; i < JPEGConstants::BLOCK_SIZE; ++i) {
                    for (int j = 0; j < JPEGConstants::BLOCK_SIZE; ++j) {
                        double value = result[i][j] + (channel == 0 ? 128 : 0);
                        block.data[i][j][channel] = clamp(round(value), 0, 255);
                    }
                }
            }
        }

        map<int, string> huffmanEncode(const vector<PixelBlock>& blocks) {
            vector<int> frequencies(JPEGConstants::COLOR_RANGE, 0);
            vector<int> values(JPEGConstants::COLOR_RANGE);

            for (int i = 0; i < JPEGConstants::COLOR_RANGE; ++i) {
                values[i] = i;
            }

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

        vector<vector<array<int, 3>>> convertToRGB(const vector<PixelBlock>& blocks) {
            vector<vector<array<int, 3>>> rgbImage(height, vector<array<int, 3>>(width));

            int blockIndex = 0;
            for (int blockY = 0; blockY < height; blockY += JPEGConstants::BLOCK_SIZE) {
                for (int blockX = 0; blockX < width; blockX += JPEGConstants::BLOCK_SIZE) {
                    const auto& block = blocks[blockIndex++];

                    for (int y = 0; y < JPEGConstants::BLOCK_SIZE && (blockY + y) < height; ++y) {
                        for (int x = 0; x < JPEGConstants::BLOCK_SIZE && (blockX + x) < width; ++x) {
                            int imageY = blockY + y;
                            int imageX = blockX + x;

                            // YCbCr to RGB conversion
                            double Y = block.data[y][x][0];
                            double Cb = block.data[y][x][1] - 128.0;
                            double Cr = block.data[y][x][2] - 128.0;

                            int R = round(Y + 1.402 * Cr);
                            int G = round(Y - 0.344136 * Cb - 0.714136 * Cr);
                            int B = round(Y + 1.772 * Cb);

                            rgbImage[imageY][imageX] = {
                                clamp(R, 0, 255),
                                clamp(G, 0, 255),
                                clamp(B, 0, 255)
                            };
                        }
                    }
                }
            }

            return rgbImage;
        }

        static int clamp(int value, int min, int max) {
            return value < min ? min : (value > max ? max : value);
        }

        void compressImage(const string& outputFilename) {
            auto ycbcrImage = convertToYCbCr();

            auto blocks = divideIntoBlocks(ycbcrImage);

            for (auto& block : blocks) {
                applyDCT(block);
            }

            for (auto& block : blocks) {
                quantize(block);
            }

            auto huffmanCodes = huffmanEncode(blocks);

            for (auto& block : blocks) {
                dequantize(block);
                applyIDCT(block);
            }

            auto rgbImage = convertToRGB(blocks);

            saveImage(outputFilename, rgbImage);
        }
    };

    vector<vector<double>> Image::dctMatrix;
    vector<vector<double>> Image::dctMatrixTranspose;

    int main() {
        try {
            string inputFile = "../samples/sample_soho.bmp";
            string outputFile = "../results/compressed.bmp";

            Image img(inputFile);
            img.compressImage(outputFile);

            cout << "Image compression completed successfully!" << endl;
            return 0;
        }
        catch (const ImageException& e) {
            cerr << "Image processing error: " << e.what() << endl;
            return 1;
        }
        catch (const exception& e) {
            cerr << "Unexpected error: " << e.what() << endl;
            return 2;
        }
    }
