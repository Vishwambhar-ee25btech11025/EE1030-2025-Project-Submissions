# Overview
I have used hybrid c+python for this project. The following presents you the way in which you can compile files and view the image.

## C
It conatins my algorithm to find first k singular values.
There are two codes one for RGB(svd_final_color.c) and another for Greyscale(svd_final_greyscale).

### Prerequisites
Do these in the debian environment.
1. ```pkg update && pkg upgrade```
2. ```pkg install clang```
3. ```pkg install build-essential```

### Compiling those C files
First cd to the folder c_backend/.

To compile the greyscale file use the following command. 

```gcc -shared -o svd_final_greyscale.so -fPIC svd_final_greyscale.c -lm```

To compile the color file use the following command.

```gcc -shared -o svd_final_color.so -fPIC svd_final_color.c -lm```

These two command create a shared object (.so files) which are used by python to call C functions.

## Python
It contains command to I/O the image into/from matrix containing the information about each pixel.

### Prerequisites

Do these in the debian environment.

1. ```pkg install python```
2. ```pip install numpy scipy matplotlib pillow```

### Executing Python codes

cd to the python_frontend folder.
Then use the following command:

```python3 image_compressor.py```

The image will be saved to figs folder as 'recent.png'

To view the image first copy it to downloads folder in your file manager. Use the following command in figs folder:

```cp recent.png/sdcard/Download/```

## Conclusion 
You ac view the image in your file manager.
