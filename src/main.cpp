#include <iostream>
#include "image.hpp"
#include "exceptions.hpp"

int main() {
    try {
        std::string inputFile = "./samples/sample_soho.bmp";
        std::string outputFile = "./results/compressed.bmp";

        Image img(inputFile);
        img.compressImage(outputFile);

        std::cout << "Image compression completed successfully!" << std::endl;
        return 0;
    }
    catch (const ImageException& e) {
        std::cerr << "Image processing error: " << e.what() << std::endl;
        return 1;
    }
    catch (const std::exception& e) {
        std::cerr << "Unexpected error: " << e.what() << std::endl;
        return 2;
    }
}
