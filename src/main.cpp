#include <iostream>
#include <filesystem>
#include <algorithm>
#include "image.hpp"
#include "exceptions.hpp"

void printUsage(const char* programName) {
    std::cerr << "\nImage Compression Tool\n"
              << "=====================\n\n"
              << "Description:\n"
              << "  This is a simple command-line tool for compressing BMP images using a JPEG-like compression algorithm. The tool\n"
              << "  is designed to demonstrate basic image compression techniques and is not intended for professional use.\n\n"
              << "Usage:\n"
              << "  " << programName << " <input_file> <output_file>\n\n"
              << "Arguments:\n"
              << "  <input_file>    Path to the input BMP image (uncompressed)\n"
              << "  <output_file>   Path to save the compressed image in BMP format\n\n"
              << "Example:\n"
              << "  " << programName << " ./samples/input.bmp ./results/output.bmp\n\n"
              << "Notes:\n"
              << "  - This tool only supports BMP images as input and output.\n"
              << "  - The output file will be saved in a compressed BMP format.\n"
              << "  - Make sure the specified directories for input and output exist or are writable.\n\n"
              << "Supported formats:\n"
              << "  - Input: BMP (Bitmap)\n"
              << "  - Output: Compressed BMP\n\n";
}

bool validateInputFile(const std::string& path) {
    namespace fs = std::filesystem;

    if (!fs::exists(path)) {
        std::cerr << "Error: Input file '" << path << "' does not exist.\n";
        return false;
    }

    if (fs::is_directory(path)) {
        std::cerr << "Error: '" << path << "' is a directory, not a file.\n";
        return false;
    }

    std::string extension = fs::path(path).extension().string();
    std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

    if (extension != ".bmp") {
        std::cerr << "Error: Input file must be a BMP image (got '" << extension << "').\n";
        return false;
    }

    return true;
}

bool validateOutputPath(const std::string& path) {
    namespace fs = std::filesystem;

    try {
        // Create output directory if it doesn't exist
        fs::path outputPath = fs::path(path).parent_path();
        if (!outputPath.empty() && !fs::exists(outputPath)) {
            fs::create_directories(outputPath);
        }

        // Check if output file extension is correct
        std::string extension = fs::path(path).extension().string();
        std::transform(extension.begin(), extension.end(), extension.begin(), ::tolower);

        if (extension != ".bmp") {
            std::cerr << "Error: Output file must have .bmp extension (got '" << extension << "').\n";
            return false;
        }

        return true;
    }
    catch (const fs::filesystem_error& e) {
        std::cerr << "Error: Cannot create output directory: " << e.what() << "\n";
        return false;
    }
}

int main(int argc, char* argv[]) {
    try {
        if (argc != 3) {
            std::cerr << "Error: Incorrect number of arguments.\n";
            printUsage(argv[0]);
            return 1;
        }

        std::string inputFile = argv[1];
        std::string outputFile = argv[2];

        // Validate input and output paths
        if (!validateInputFile(inputFile) || !validateOutputPath(outputFile)) {
            printUsage(argv[0]);
            return 1;
        }

        // Process the image
        std::cout << "Processing image '" << inputFile << "'...\n";
        Image img(inputFile);
        img.compressImage(outputFile);
        std::cout << "\nImage compression completed successfully!\n"
                  << "Output saved to: " << outputFile << "\n";
        return 0;
    }
    catch (const ImageException& e) {
        std::cerr << "\nImage processing error: " << e.what() << "\n";
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "\nUnexpected error: " << e.what() << "\n";
        return 2;
    }
}
