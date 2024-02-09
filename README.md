
# Simple educational JPEG compress from scratch

- This code is not meant for production, It's only for educational reasons

## Progress

#### Loading Bitmap picture ✅
![loading bitmap](https://i.imgur.com/NzZltzY.png)
#### Converting to YCbCr ✅
![ycbcr](https://i.imgur.com/u2qtxQn.png)
#### Dividing Into blocks ✅
#### Shifting Y value ✅
![shift](https://i.imgur.com/CFjirZe.png)
#### Discrete Cosine Transform (DCT) ✅
![dct](https://i.imgur.com/tlFD87J.png)
#### Quantization ✅
![quant](https://i.imgur.com/xxbjQ6Z.png)
#### IDCT ✅
![idct](https://i.imgur.com/k26LIv1.png)
#### Shifting Y value back ❌
#### Converting YCbCr back to RGB ✅
#### Huffman encoding ❓
#### Combining blocks back to final image ✅

## Needs to repair/update

#### Check/Debug Huffman encoding
#### When shifting pixel values back, accidentally shifting RGB channel instead of Y value
#### Add argument instead of hardcoded path

