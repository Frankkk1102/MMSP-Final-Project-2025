# MMSP Final Project 

## 1. Project Overview

這次的project是在設計一套影像壓縮、解壓縮的系統，並用decoder/encoder來執行，並完成mode0到mode3的檢驗。
而mode0到mode3都會用diff/cmp來驗證內容是否正確
---

## 2. System Architecture and Implementation Progress

### Overall Pipeline

BMP Image (RGB)  
→ Reversible Color Transform (Y, Cb, Cr)  
→ 8×8 Block Processing  
→ DCT  
→ Quantization  
→ Zigzag Scan  
→ Run-Length Encoding (Mode 2)  
→ Huffman Encoding (Mode 3)  
→ Bitstream / Codebook  

而 Decoder 端則依照相反順序進行處理，將影像重建回 RGB BMP 格式。

---

#### Mode 0 – 原始RGB重建
- Separates RGB channels into text files
- Reconstructs BMP image from R/G/B channels
- Verified using `cmp` for exact match

#### Mode 1 – 量化、重建
- Applies DCT and quantization using provided quantization tables
- Supports decoding with and without error feedback
- Computes pixel-level SQNR (R, G, B channels)
- Reconstruction verified using `cmp`

#### Mode 2 – Run-Length Encoding (RLE)
- Implements RLE on zigzag-scanned quantized coefficients
- Supports both ASCII and binary formats
- Decoder reconstructs images identical to Mode 1 output
- Verified using `cmp`

#### Mode 3 – Huffman Coding
- Builds Huffman codebook from RLE symbols
- Supports ASCII and binary Huffman bitstreams
- Decoder reconstructs images identical to Mode 1 and Mode 2 outputs
- Verified using `cmp`
---

## 3. 驗證 and GitHub Workflow

A GitHub Actions workflow (`.github/workflows/main.yml`) is configured to run
automatically on every push. The workflow performs the following steps:

1. Compile encoder and decoder
2. Execute Mode 1 (Quantization) pipeline
3. Execute Mode 3 (Huffman) pipeline
4. Verify reconstructed images using `cmp`
5. Upload generated outputs as artifacts

---

## 4. Artifacts Description

The GitHub Actions workflow uploads the following artifacts:

- Reconstructed BMP images
- Quantized coefficient files (`*.raw`)
- Huffman codebook files (`*.txt`)
- Huffman bitstreams (`*.bin`)

These artifacts are used for grading and demonstrate the correctness of the
implementation.

---

## 5. 心得、感想
  這次的final project其實對我自己來說是蠻有難度的，其實老師在第三次作業開始就有讓我們練習用gitgub、熟悉一些功能的使用，但因為是分組作業，雖然有一起做，但我剛好都是負責程式的部分，所以有關workflow、README的內容我是比較不熟的，但這次的作業就是讓你從頭到尾自己完成，我也去問了組員他們當初在做這些部分的時候遇到的困難。我自己是要把我電腦執行的程式push到github的時候遇到了一些困難，因為檔名的問題，在我從github clone完後多了一個資料夾，我一直push我原本的資料夾，然後一直找不到檔案，最後問了同學也問了AI就有解決了。
  這次的project，也確實讓我學到很多內容，不管是DCT變換還是Huffman Coding我都有更深入的了解，但要說最熟的還是打在terminal上的指令，因為過程中一直出錯，每一次都重打基本上是背起來了，也希望自己能把這學期學到的內容好好的吸收、真的變成我的技能，在以後面對到相關的內容時，可以有自信的完成。
