#include <iostream>
#include <vector>
#include <tuple>
#include <array>
#include <fstream>
#include <cmath>

using namespace std;

// https://arxiv.org/ftp/arxiv/papers/1405/1405.6147.pdf
// https://stackoverflow.com/questions/4931621/jpeg-encoding-technique

class Image {
public:
    int width, height;
    std::vector<std::vector<std::array<int, 3>>> pixels;

    Image(const std::string &filename) {
        loadImage(filename);
    }

    void loadImage(const std::string &filename) {
        // Loading image from file
        std::ifstream file(filename, std::ios::binary);
        cout << "LOADING  IMAGE DATA..." << endl;

        // BMP format
        // Skipping BMP header (54 bytes) 
        // https://www.ece.ualberta.ca/~elliott/ee552/studentAppNotes/2003_w/misc/bmp_file_format/bmp_file_format.htm
        file.seekg(54);

        // Width on 18. byte 
        // Height on 22. byte
        file.seekg(18);
        file.read(reinterpret_cast<char*>(&width), 4);  // retyped to char
        file.read(reinterpret_cast<char*>(&height), 4);

        cout << "WIDTH: " << width << endl;
        cout << "HEIGHT: " << height << endl;

        // Loading pixels from picture
        cout << "LOADING PIXELS FROM PICTURE..." << endl;
        pixels.resize(height, vector<array<int, 3>>(width));
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                char colors[3];
                file.read(colors, 3);
                pixels[y][x] = {colors[2], colors[1], colors[0]};  // BMP saves colors as BGR
            }
        }
    }


    // Dividing into 8x8 blocks
    vector<vector<vector<array<int, 3>>>> divideIntoBlocks() {
        
        vector<vector<vector<array<int, 3>>>> blocks; 
        cout << "DIVIDING INTO BLOCKS..." << endl;
        for (int i = 0; i < height; i += 8) {
            for (int j = 0; j < width; j += 8) {

                vector<vector<array<int, 3>>> block; // 8x8 pixels block
                
                // if y is not limited => overflow 
                for (int y = i; y < min(i + 8, height); ++y) {

                    vector<array<int, 3>> row;
                    // if x is not limited => overflow 
                    for (int x = j; x < min(j + 8, width); ++x) {
                        row.push_back(pixels[y][x]);
                    }
                    block.push_back(row);
                }
                blocks.push_back(block);
            }
        }
        
        cout << "NUMBER OF BLOCKS: " << blocks.size() << endl;

        return blocks;
    }

    // Converting to YCbCr
    void convertToYCbCr() {
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                int R = pixels[y][x][0];
                int G = pixels[y][x][1];
                int B = pixels[y][x][2];

                int Y = 0.299 * R + 0.587 * G + 0.114 * B;
                int Cb = 128 - 0.168736 * R - 0.331264 * G + 0.5 * B;
                int Cr = 128 + 0.5 * R - 0.418688 * G - 0.081312 * B;

                pixels[y][x] = {Y, Cb, Cr};  // Saving in YCbCr format
            }
        }
    }

    // Shifting pixel values from [0 - 255] to [-128 - 127]
    void shiftPixelValues() {
        std::cout << "SHIFTING PIXEL VALUES..." << std::endl;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                pixels[y][x][0] -= 128; // just needs substracting 128 from Y, Cb, Cr
                pixels[y][x][1] -= 128;
                pixels[y][x][2] -= 128;
            }
        }   
    }

    // Discrete Cosine Transform on each block 
    // From left to right, top to bottom
    vector<vector<vector<array<int, 3>>>> applyDCT(vector<vector<vector<array<int, 3>>>> &blocks) {
        cout << "APPLYING DCT..." << endl;

        vector<vector<vector<array<int, 3>>>> dctBlocks;
        dctBlocks.reserve(blocks.size());

        for (auto &block : blocks) {
            vector<vector<array<int, 3>>> dctBlock(8, vector<array<int, 3>>(8));
            // u frequency
            for (int u = 0; u < 8; ++u) {
                // v frequency
                for (int v = 0; v < 8; ++v) {

                    double sumY = 0.0, sumCb = 0.0, sumCr = 0.0;

                    for (int x = 0; x < 8; ++x) {
                        for (int y = 0; y < 8; ++y) {
                            int Y = block[y][x][0];
                            int Cb = block[y][x][1];
                            int Cr = block[y][x][2];
                            // cosin transformation for selected pixel and frequency
                            double cosineX = cos((2 * x + 1) * u * M_PI / 16.0);
                            double cosineY = cos((2 * y + 1) * v * M_PI / 16.0);

                            sumY += Y * cosineX * cosineY;
                            sumCb += Cb * cosineX * cosineY;
                            sumCr += Cr * cosineX * cosineY;
                        }
                    }

                    // constants
                    double Cu = (u == 0) ? sqrt(2.0) / 2.0 : 1.0;
                    double Cv = (v == 0) ? sqrt(2.0) / 2.0 : 1.0;

                    // dct value for Y/Cb/Cr
                    double dctValueY = 0.25 * Cu * Cv * sumY;
                    double dctValueCb = 0.25 * Cu * Cv * sumCb;
                    double dctValueCr = 0.25 * Cu * Cv * sumCr;
                    // saving computed values
                    dctBlock[u][v] = {static_cast<int>(dctValueY), static_cast<int>(dctValueCb), static_cast<int>(dctValueCr)};  // Saving in DCT format
                }
            }
            dctBlocks.push_back(dctBlock);
        }
        return dctBlocks;
    }

    // Each block compressed through quantization
    void quantize(vector<vector<vector<array<int, 3>>>> &blocks) {
        cout << "QUANTIZING..." << endl;

        // Quantization matrix for Y component 
        vector<vector<int>> quantizationMatrixY = {
            {16, 11, 10, 16, 24, 40, 51, 61},
            {12, 12, 14, 19, 26, 58, 60, 55},
            {14, 13, 16, 24, 40, 57, 69, 56},
            {14, 17, 22, 29, 51, 87, 80, 62},
            {18, 22, 37, 56, 68, 109, 103, 77},
            {24, 35, 55, 64, 81, 104, 113, 92},
            {49, 64, 78, 87, 103, 121, 120, 101},
            {72, 92, 95, 98, 112, 100, 103, 99}
        };

        // Quantization matrix for Cb and Cr components 
        vector<vector<int>> quantizationMatrixCbCr = {
            {17, 18, 24, 47, 99, 99, 99, 99},
            {18, 21, 26, 66, 99, 99, 99, 99},
            {24, 26, 56, 99, 99, 99, 99, 99},
            {47, 66, 99, 99, 99, 99, 99, 99},
            {99, 99, 99, 99, 99, 99, 99, 99},
            {99, 99, 99, 99, 99, 99, 99, 99},
            {99, 99, 99, 99, 99, 99, 99, 99},
            {99, 99, 99, 99, 99, 99, 99, 99}
        };


        for (auto &block : blocks) {
            for (int y = 0; y < 8; ++y) {
                for (int x = 0; x < 8; ++x) {
                    // extracting Y, Cb a Cr values from block
                    int Y = block[y][x][0];
                    int Cb = block[y][x][1];
                    int Cr = block[y][x][2];

                    // Quantizing Y, Cb, Cr values with matrices
                    Y = round(Y / quantizationMatrixY[y][x]);
                    Cb = round(Cb / quantizationMatrixCbCr[y][x]);
                    Cr = round(Cr / quantizationMatrixCbCr[y][x]);

                    // saving to block
                    block[y][x] = {Y, Cb, Cr};
                }
            }
        }

    }

    void huffmanEncoding() {
        
    }

    // Compressed image is reconstructed through reverse process (IDCT)
    void reconstructImage() {

    }
};

int main() {
    Image img("bmp_24.bmp");
    img.convertToYCbCr();
    img.shiftPixelValues();
    auto blocks = img.divideIntoBlocks();
    auto dctBlocks = img.applyDCT(blocks);
    img.quantize(dctBlocks);
    //img.huffmanEncoding();
    //img.reconstructImage();

    return 0;
}