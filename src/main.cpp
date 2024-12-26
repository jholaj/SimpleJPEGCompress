#include <iostream>
#include <filesystem>
#include <algorithm>
#include <string>
#include <cstdlib>
#include "image.hpp"
#include "exceptions.hpp"
#include "constants.hpp"

void printUsage(const char* programName) {
    std::cerr << "\nImage Compression Tool\n"
              << "=====================\n\n"
              << "Description:\n"
              << "  A simple command-line tool for compressing BMP images using a JPEG-like compression algorithm.\n\n"
              << "Usage:\n"
              << "  " << programName << " [options] <input_file> <output_file>\n\n"
              << "Arguments:\n"
              << "  <input_file>    Path to the input BMP image (uncompressed)\n"
              << "  <output_file>   Path to save the compressed image in BMP format\n\n"
              << "Options:\n"
              << "  -b, --block-size <size>         Block size for compression (default: 8)\n"
              << "  -q, --quality-factor <factor>   Quality factor for compression (default: 1.0)\n\n"
              << "Example:\n"
              << "  " << programName << " ./samples/input.bmp ./results/output.bmp -q 0.1 -b 8 \n\n";
}

bool parseArguments(int argc, char* argv[], int& blockSize, double& qualityFactor,
                    std::string& inputFile, std::string& outputFile) {
    blockSize = JPEGConstants::BLOCK_SIZE;
    qualityFactor = JPEGConstants::QUALITY_FACTOR;

    int positionalArgs = 0;

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-b" || arg == "--block-size") {
            if (i + 1 < argc) {
                blockSize = std::stoi(argv[++i]);
            } else {
                std::cerr << "Error: Missing value for block size.\n";
                return false;
            }
        } else if (arg == "-q" || arg == "--quality-factor") {
            if (i + 1 < argc) {
                qualityFactor = std::stod(argv[++i]);
            } else {
                std::cerr << "Error: Missing value for quality factor.\n";
                return false;
            }
        } else if (arg[0] != '-') {
            if (positionalArgs == 0) {
                inputFile = arg;
                positionalArgs++;
            } else if (positionalArgs == 1) {
                outputFile = arg;
                positionalArgs++;
            } else {
                std::cerr << "Error: Too many positional arguments.\n";
                return false;
            }
        } else {
            std::cerr << "Error: Unknown argument '" << arg << "'.\n";
            return false;
        }
    }

    if (positionalArgs < 2) {
        std::cerr << "Error: Missing input or output file.\n";
        return false;
    }

    return true;
}

int main(int argc, char* argv[]) {
    int blockSize;
    double qualityFactor;
    std::string inputFile;
    std::string outputFile;

    if (!parseArguments(argc, argv, blockSize, qualityFactor, inputFile, outputFile)) {
        printUsage(argv[0]);
        return 1;
    }

    std::cout << "Settings:\n"
              << "  Block Size: " << blockSize << "\n"
              << "  Quality Factor: " << qualityFactor << "\n";

    try {
        // Process the image
        std::cout << "Processing image '" << inputFile << "'...\n";
        Image img(inputFile, blockSize, qualityFactor);
        img.compressImage(outputFile);
        std::cout << "\nImage compression completed successfully!\n"
                  << "Output saved to: " << outputFile << "\n";
        return 0;
    } catch (const ImageException& e) {
        std::cerr << "\nImage processing error: " << e.what() << "\n";
        return 1;
    } catch (const std::exception& e) {
        std::cerr << "\nUnexpected error: " << e.what() << "\n";
        return 2;
    }
}
