# Software Project(Image Compression)
## Overview
This project highlights the use of Singular Value Decomposition of matrices to compress an image to required degree.

## Contents
- codes/
- figs/
- tables/
- report.tex
- report.pdf

## Roadmap
```
ee25btech11025/
|
SoftwareAssignment/
|_______codes/____c_libs/(empty)
|         |_______c_main/(empty)
|         |_______python_driver/(empty)
|         |_______hybrid_c_python/
|         |             |____c_backend/
|         |             |         |______svd_final_color.c
|         |             |         |______svd_final_greyscale.c
|         |             |         |______svd_final_color.so
|         |             |         |______svd_final_greyscale.so
|         |             |
|         |             |
|         |             |_____python_frontend/
|         |                       |______image_compressor.py
|         |              
|         |_________README_code.md        
|
|
|________figs/_____einstein.jpg
|           |______recent.png
|           |______flower.png
|           |______globe.jpg
|           |______greyscale.png
|           |______logo.png
|           |______sample1.png
|           |______sample1_stats.png
|           |______sample2.png
|           |______sample2_stats.png
|           |______sample3.png
|           |______sample3_stats.png
|           |______sample4.png
|           |______sample_stats.png
|
|
|_______tables/___frobenius.tex
|
|_______report.tex
|_______report.pdf
|_______report.aux
|_______report.log
|_______README.md
```

## How to compile report.tex
### Prerequisites
Install pdflatex packages in the debian environment of your termux app.
Use the following command:

```apt install texlive-full gnumeric```

Also git clone the repo

```git clone https://github.com/username/EE1030-2025-Project-Submissions```

### Procedure
cd to my folder. Then use the following command:

```pdflatex report.tex```

Then you will be able to see the compiled pdf.






