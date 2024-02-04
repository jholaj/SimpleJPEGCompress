#include <iostream>
#include <vector>
#include <array>
#include <fstream>
#include <cmath>
#include <queue>
#include <map>

using namespace std;

// https://arxiv.org/ftp/arxiv/papers/1405/1405.6147.pdf
// https://stackoverflow.com/questions/4931621/jpeg-encoding-technique

struct Node {
    int data;
    unsigned freq;
    Node *left, *right;

    Node(int data, unsigned freq) {
        left = right = NULL;
        this->data = data;
        this->freq = freq;
    }
};

// For comparison of two heap nodes (needed in min heap)
struct compare {
    bool operator()(Node* l, Node* r) {
        return (l->freq > r->freq);
    }
};


struct EncodedPixel {
    string encodedY;
    string encodedCb;
    string encodedCr;
};

class HuffmanTree {
    private:
    Node* root;

    // recursive generating huffman tree
    void generateCodes(Node* root, string str, map<int, string> &huffmanCodes) {
        if (!root)
            return;

        if (root->data != -1)
            huffmanCodes[root->data] = str;

        generateCodes(root->left, str + "0", huffmanCodes);
        generateCodes(root->right, str + "1", huffmanCodes);
    }

public:
    map<int, string> encode(int data[], int freq[], int size) {
        Node *left, *right, *top;

        // min heap -> all symbols from data[]
        priority_queue<Node*, vector<Node*>, compare> minHeap;

        for (int i = 0; i < size; ++i)
            minHeap.push(new Node(data[i], freq[i]));

        // Iteration until heap size not 1
        while (minHeap.size() != 1) {
            left = minHeap.top();
            minHeap.pop();

            right = minHeap.top();
            minHeap.pop();

            top = new Node(-1, left->freq + right->freq);

            top->left = left;
            top->right = right;

            minHeap.push(top);
        }

        // Saving root node
        root = minHeap.top();

        // Creating codes for huffman tree
        map<int, string> huffmanCodes;
        generateCodes(root, "", huffmanCodes);

        cout << "huffmanCodes content:" << endl;
        for (const auto& pair : huffmanCodes) {
            cout << "Symbol: " << pair.first << ", Code: " << pair.second << endl;
        }

        return huffmanCodes;

    }

    int decode(const string &str, const map<int, string> &huffmanCodes) {
        Node* curr = root;
        int decodedValue = -1;

        for (char bit : str) {
            if (bit == '0') {
                curr = curr->left;
            } else {
                curr = curr->right;
            }

            if (curr->left == nullptr && curr->right == nullptr) {
                decodedValue = curr->data;
                break;
            }
        }

        return decodedValue;
    }
};

class Image {
public:
    int width, height, filesize;
    vector<vector<array<int, 3>>> pixels;

    Image(const string &filename) {
        loadImage(filename);
    }
    void loadImage(const string &filename) {
        // Loading image from file
        ifstream file(filename, ios::binary);
        cout << "LOADING IMAGE DATA..." << endl;

        // BMP format
        // Skipping BMP header (54 bytes) 
        // https://www.ece.ualberta.ca/~elliott/ee552/studentAppNotes/2003_w/misc/bmp_file_format/bmp_file_format.htm
        file.seekg(0, std::ios::end);
        filesize = file.tellg();

        // Width on 18. byte 
        // Height on 22. byte
        file.seekg(18);
        file.read(reinterpret_cast<char*>(&width), 4);  // retyped to char
        file.read(reinterpret_cast<char*>(&height), 4);

        cout << "WIDTH: " << width << endl;
        cout << "HEIGHT: " << height << endl;
        cout << "FILE SIZE: " << filesize << endl;

        // padding
        int padding = ((4 - (width * 3) % 4) % 4) %4;

        // Loading pixels from picture
        cout << "LOADING PIXELS FROM PICTURE..." << endl;
        pixels.resize(height, vector<array<int, 3>>(width));
        for (int y = height-1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                unsigned char colors[3];
                file.read(reinterpret_cast<char*>(colors), 3);
                pixels[y][x] = {colors[0], colors[2], colors[1]};
                //cout << static_cast<int>(colors[2]) << " " << static_cast<int>(colors[1]) << " " << static_cast<int>(colors[0]) << endl;
            }
            file.seekg(padding, std::ios::cur);
        }
        file.close();
        saveImage("01_TEST_LOADED_PICTURE.bmp", pixels);
    }

    void saveImage(const std::string &filename, const std::vector<std::vector<std::array<int, 3>>> &image) {
        std::ofstream file(filename, std::ios::binary);
        if (!file.is_open()) {
            std::cerr << "Failed to open file for writing: " << filename << std::endl;
            return;
        }

        const int width = image[0].size();
        const int height = image.size();

        int padding = (4 - (width * 3) % 4) % 4;

        // BMP header
        const int fileSize = 54 + (3 * width + padding) * height;
        cout << "COMPRESSED SIZE: " << fileSize << endl;
        const int dataOffset = 54;

        unsigned char header[54] = {
            'B', 'M', // Signature
            static_cast<unsigned char>(fileSize), static_cast<unsigned char>(fileSize >> 8), static_cast<unsigned char>(fileSize >> 16), static_cast<unsigned char>(fileSize >> 24), // File size
            0, 0, 0, 0, // Reserved
            static_cast<unsigned char>(dataOffset), static_cast<unsigned char>(dataOffset >> 8), static_cast<unsigned char>(dataOffset >> 16), static_cast<unsigned char>(dataOffset >> 24), // Data offset
            40, 0, 0, 0, // Header size
            static_cast<unsigned char>(width), static_cast<unsigned char>(width >> 8), static_cast<unsigned char>(width >> 16), static_cast<unsigned char>(width >> 24), // Image width
            static_cast<unsigned char>(height), static_cast<unsigned char>(height >> 8), static_cast<unsigned char>(height >> 16), static_cast<unsigned char>(height >> 24), // Image height
            1, 0, // Planes
            24, 0, // Bits per pixel
            0, 0, 0, 0, // Compression
            0, 0, 0, 0, // Image size
            0, 0, 0, 0, // X pixels per meter
            0, 0, 0, 0, // Y pixels per meter
            0, 0, 0, 0, // Colors in color table
            0, 0, 0, 0, // Important color count
        };

        // Write BMP header
        file.write(reinterpret_cast<const char*>(header), 54);

        // Write image data (BGR format)
        for (int y = height-1; y >= 0; --y) {
            for (int x = 0; x < width; ++x) {
                unsigned char pixel[3] = {
                    static_cast<unsigned char>(image[y][x][2]), // Blue
                    static_cast<unsigned char>(image[y][x][1]), // Green
                    static_cast<unsigned char>(image[y][x][0]), // Red
                };
                file.write(reinterpret_cast<const char*>(pixel), 3);
            }

            char padding_bytes[3] = {0};  // Fill padding with zeros
            file.write(padding_bytes, padding);
        }

        file.close();
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

                //cout << R << endl;

                int Y = 0.299 * R + 0.587 * G + 0.114 * B;
                int Cb = 128 - 0.168736 * R - 0.331264 * G + 0.5 * B;
                int Cr = 128 + 0.5 * R - 0.418688 * G - 0.081312 * B;

                pixels[y][x] = {Y, Cb, Cr};
            }
        }
        saveImage("02_TEST_TO_YCBCR.bmp", pixels);
    }

    void shiftPixelValues() {
        cout << "SHIFTING PIXEL VALUES..." << endl;
        for (int y = 0; y < height; ++y) {
            for (int x = 0; x < width; ++x) {
                pixels[y][x][0] -= 128; // just needs substracting 128 from R,G,B
                pixels[y][x][1] -= 128;
                pixels[y][x][2] -= 128;
            }
        }
        saveImage("03_TEST_SHIFTING.bmp", pixels); 
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
        auto image = combineBlocks(dctBlocks);
        saveImage("04_TEST_DCT.bmp", image);
        return dctBlocks;
    }

    // Each block compressed through quantization
    vector<vector<vector<array<int, 3>>>> quantize(vector<vector<vector<array<int, 3>>>> &blocks) {
        cout << "QUANTIZING..." << endl;

        vector<vector<vector<array<int, 3>>>> quantizedBlocks = blocks;

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


        for (auto &block : quantizedBlocks) {
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

        auto img = combineBlocks(quantizedBlocks);
        saveImage("05_TEST_QUANT.bmp", img);

        return quantizedBlocks;

    }


    map<int, string> huffmanEncoding(vector<vector<vector<array<int, 3>>>> &blocks) {
        cout << "HUFFMAN ENCODING..." << endl;

        HuffmanTree huffmanTree;

        // Prepare frequency array
        int freq[256] = {0};

        // Count frequency of each value in blocks
        for (auto &block : blocks) {
            for (int y = 0; y < 8; ++y) {
                for (int x = 0; x < 8; ++x) {
                    // extracting Y, Cb a Cr values from block
                    int Y = block[y][x][0];
                    int Cb = block[y][x][1];
                    int Cr = block[y][x][2];

                    //cout << "Y: " << Y << ", Cb: " << Cb << ", Cr: " << Cr << endl;

                    // Increase frequency count
                    freq[Y]++;
                    freq[Cb]++;
                    freq[Cr]++;
                }
            }
        }    

        int data[256];
        for (int i = 0; i < 256; i++) {
            data[i] = i;
        }

        // Build Huffman Tree and get the code map
        map<int, string> huffmanCodes = huffmanTree.encode(data, freq, 256);

        // Return the map of Huffman codes
        return huffmanCodes;
    }


    vector<vector<vector<array<int, 3>>>> applyIDCT(vector<vector<vector<array<int, 3>>>> &dctBlocks) {
        vector<vector<vector<array<int, 3>>>> blocks;
        blocks.reserve(dctBlocks.size());

        for (auto &dctBlock : dctBlocks) {
            vector<vector<array<int, 3>>> block(8, vector<array<int, 3>>(8));
            // x position in block
            for (int x = 0; x < 8; ++x) {
                // y position in block
                for (int y = 0; y < 8; ++y) {

                    double sumY = 0.0, sumCb = 0.0, sumCr = 0.0;

                    // u frequency
                    for (int u = 0; u < 8; ++u) {
                        // v frequency
                        for (int v = 0; v < 8; ++v) {
                            double Cu = (u == 0) ? sqrt(2.0) / 2.0 : 1.0;
                            double Cv = (v == 0) ? sqrt(2.0) / 2.0 : 1.0;

                            double cosineX = cos((2 * x + 1) * u * M_PI / 16.0);
                            double cosineY = cos((2 * y + 1) * v * M_PI / 16.0);

                            double dctValueY = dctBlock[u][v][0];
                            double dctValueCb = dctBlock[u][v][1];
                            double dctValueCr = dctBlock[u][v][2];

                            sumY += Cu * Cv * dctValueY * cosineX * cosineY;
                            sumCb += Cu * Cv * dctValueCb * cosineX * cosineY;
                            sumCr += Cu * Cv * dctValueCr * cosineX * cosineY;
                        }
                    }

                    // Discretization and rounding
                    int Y = static_cast<int>(sumY);
                    int Cb = static_cast<int>(sumCb);
                    int Cr = static_cast<int>(sumCr);

                    Y = max(0, min(255, Y));
                    Cb = max(0, min(255, Cb));
                    Cr = max(0, min(255, Cr));

                    // Saving computed YCbCr values to block
                    block[y][x] = {Y, Cb, Cr};
                }
            }
            blocks.push_back(block);
        }
        auto img = combineBlocks(blocks);
        saveImage("06_TEST_IDCT.bmp", img);
        return blocks;
    }


    vector<vector<vector<array<int, 3>>>> shiftPixelValuesBack(vector<vector<vector<array<int, 3>>>> &blocks) {
        for (auto &block : blocks) {
            for (auto &row : block) {
                for (auto &pixel : row) {
                    pixel[0] += 128; // Y component
                    //cout << " Y: " << pixel[0]<< " Cb: " << pixel[1] << " Cr: " << pixel[2];
                }
            }
        }

        auto img = combineBlocks(blocks);
        saveImage("07_TEST_SHIFTING_BACK.bmp", img);

        return blocks;
    }

    vector<vector<vector<array<int, 3>>>> convertYCbCrToRGB(const vector<vector<vector<array<int, 3>>>>& blocksYCbCr) {
        // Create blocksRGB with the same structure as blocksYCbCr
        vector<vector<vector<array<int, 3>>>> blocksRGB(blocksYCbCr.size());
        for (size_t i = 0; i < blocksYCbCr.size(); ++i) {
            blocksRGB[i].resize(blocksYCbCr[i].size());
            for (size_t j = 0; j < blocksYCbCr[i].size(); ++j) {
                blocksRGB[i][j].resize(blocksYCbCr[i][j].size());
            }
        }

        // Convert YCbCr to RGB
        for (size_t i = 0; i < blocksYCbCr.size(); ++i) {
            for (size_t j = 0; j < blocksYCbCr[i].size(); ++j) {
                for (size_t k = 0; k < blocksYCbCr[i][j].size(); ++k) {
                    int Y = blocksYCbCr[i][j][k][0];
                    int Cb = blocksYCbCr[i][j][k][1];
                    int Cr = blocksYCbCr[i][j][k][2];

                    int R = Y + 1.402 * (Cr - 128);
                    int G = Y - 0.344136 * (Cb - 128) - 0.714136 * (Cr - 128);
                    int B = Y + 1.772 * (Cb - 128);

                    blocksRGB[i][j][k] = {R, G, B};
                }
            }
        }

        auto img = combineBlocks(blocksRGB);
        saveImage("08_TEST_TO_RGB.bmp", img);

        return blocksRGB;
    }

    vector<vector<array<int, 3>>> combineBlocks(const vector<vector<vector<array<int, 3>>>>& blocks) {
        vector<vector<array<int, 3>>> combinedImage(height, vector<array<int, 3>>(width));

        int blockIndex = 0;
        for (int i = 0; i < height; i += 8) {
            for (int j = 0; j < width; j += 8) {
                const auto& block = blocks[blockIndex];
                for (size_t y = 0; y < block.size(); ++y) {
                    for (size_t x = 0; x < block[y].size(); ++x) {
                        combinedImage[i + y][j + x] = block[y][x];
                    }
                }
                ++blockIndex;
            }
        }

        return combinedImage;
    }



    // Compressed image is reconstructed through reverse process (IDCT)
    void reconstructImage(vector<vector<vector<array<int, 3>>>> &blocks, map<int, string> &huffmanCodes) {
        cout << "RECONSTRUCTING IMAGE..." << endl;

        // Create an instance of Huffman Tree
        HuffmanTree huffmanTree;
        cout << "DECODING HUFFMAN TREE..." << endl;
        // Decode Huffman codes back to quantized blocks
        for (auto &block : blocks) {
            for (int y = 0; y < block.size(); ++y) {
                for (int x = 0; x < block[y].size(); ++x) {
                    // Extracting Y, Cb a Cr values from block
                    int Y = block[y][x][0];
                    int Cb = block[y][x][1];
                    int Cr = block[y][x][2];

                    // Decode values using Huffman codes
                    int decodedY = huffmanTree.decode(huffmanCodes[Y], huffmanCodes);
                    int decodedCb = huffmanTree.decode(huffmanCodes[Cb], huffmanCodes);
                    int decodedCr = huffmanTree.decode(huffmanCodes[Cr], huffmanCodes);
                    //cout << decodedY << decodedCb << decodedCr << endl;

                    // Update block with decoded values
                    block[y][x] = {decodedY, decodedCb, decodedCr};
                }
            }
        }


        cout << "APPLYING INVERSE DCT..." << endl;
        // Apply inverse DCT on blocks
        auto idctBlocks = applyIDCT(blocks);

        cout << "SHIFTING PIXEL VALUES BACK..." << endl;
        // Shift pixel values back
        auto shiftedBlocks = shiftPixelValuesBack(idctBlocks);


        cout << "CONVERTING YCbCr BACK TO RGB..." << endl;
        // Convert YCbCr back to RGB
        auto rgbBlocks = convertYCbCrToRGB(shiftedBlocks);

        cout << "COMBINING BLOCKS BACK TO IMAGE..." << endl;
        // Combine blocks back to image
        auto image = combineBlocks(rgbBlocks);

        // same as to rgb image
        cout << "SAVING RECONSTRUCTED IMAGE..." << endl;
        saveImage("09_TEST_FINAL_IMAGE.bmp", image);
    }
};

int main() {
    Image img("sample_soho.bmp");
    img.convertToYCbCr();
    img.shiftPixelValues();
    auto blocks = img.divideIntoBlocks();
    auto dctBlocks = img.applyDCT(blocks);
    auto quantizedBlocks = img.quantize(dctBlocks);
    map<int, string> huffmanCodes = img.huffmanEncoding(quantizedBlocks);
    img.reconstructImage(quantizedBlocks,huffmanCodes);

    return 0;
}