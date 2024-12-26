# Image Compression using DCT and Huffman Encoding

This project is an image compression tool that uses *Discrete Cosine Transform (DCT) and Huffman Encoding* to compress images in a customizable way.

The goal is to show how compress algorhitms works.

## Features

- **DCT-based compression**: The image is first converted to the YCbCr color space, followed by applying the Discrete Cosine Transform (DCT) for frequency-domain representation.
- **Quantization**: The DCT coefficients are quantized using a user-defined quality factor and color range, which impacts the level of compression.
- **Huffman Encoding**: The quantized image data is encoded using Huffman coding to further reduce the file size.
---
- **Customizable Parameters**:
  - *Block Size*: Size of the blocks used for DCT.
  - *Quality Factor*: A value controlling the level of compression and quality.

## Compilation
1. Clone the repository:
```bash
git clone https://github.com/jholaj/SimpleJPEGCompress.git
cd SimpleJPEGCompress
```
2. Compile the project:
```bash
g++ -Iinclude -o compressor src/*.cpp
```
3. Run:
```bash
./compressor sample_images/sample.bmp results/compressed_img.bmp -q 2.0 -b 6
```

## Parameters:
- **input_image**: The path to the input image (e.g., a BMP file).
- **output_image**: The path where the compressed image will be saved (e.g., compressed_output.bmp).
- **-q <quality_factor>**: Quality factor used in quantization (default is 0.1). A lower value means less compression and better quality.
- **-b <block_size>**: Size of the blocks used for DCT (default is 8). Typical values are 8 or 16.

## Known Issues:
- Performance: The current implementation is not optimized for large images and may be slow for high-resolution inputs.
- Image Formats: The project currently supports only BMP image format for input and output.
- File Size: Current implementation does not reduce file size.

## Examples:
![No. 1 Example](https://imgur.com/XjNuLlY.png)
![No. 2 Example](https://imgur.com/H2fI8m9.png)
