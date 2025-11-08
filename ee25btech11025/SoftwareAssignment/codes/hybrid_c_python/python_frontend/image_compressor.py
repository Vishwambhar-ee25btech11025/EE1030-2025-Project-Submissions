import ctypes as ct
import numpy                                    
import matplotlib.pyplot as plot            #to plot the image
import time                                 #to calculate the runtime of my algorithm and np.linalg.svd
from PIL import Image                       #to read and convert images to matrix and vice-versa

image_name = input("Please enter which image do you want to compress(DONT FORGET TO ENTER ALONG WITH IMAGE TYPE! i.e., .jpg/.png/.jpeg): ")         #taking input of image name
type_of_image = input("Enter whether the given image is color or greyscale type: ")             #asking whether image type is color or greyscale

if type_of_image == "color":
    c_lib = ct.CDLL('../c_backend/./svd_final_color.so')                #if image is color
    dim = 3

elif type_of_image == "greyscale":
    c_lib = ct.CDLL('../c_backend/./svd_final_greyscale.so')            #if image is greyscale
    dim=2

#-------------Giving arguments type---------------------
if dim==3:
    c_lib.find_svd.argtypes = [
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=dim, flags='C_CONTIGUOUS'),
        ct.c_int,  ct.c_int,  ct.c_int,  ct.c_int,  
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=dim, flags='C_CONTIGUOUS')
    ]
    c_lib.Frobenius.argtypes = [
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=dim, flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=dim, flags='C_CONTIGUOUS'),
        ct.c_int,
        ct.c_int,
        ct.c_int
    ]
elif dim==2:
    c_lib.find_svd.argtypes = [
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=2, flags='C_CONTIGUOUS'),
        ct.c_int,  ct.c_int,  ct.c_int,
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=2, flags='C_CONTIGUOUS')
    ]
    c_lib.Frobenius.argtypes = [
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=2, flags='C_CONTIGUOUS'),
        numpy.ctypeslib.ndpointer(dtype=numpy.float64, ndim=2, flags='C_CONTIGUOUS'),
        ct.c_int,
        ct.c_int
    ]
c_lib.Frobenius.restype = ct.c_double               #giving return type

#-----------------Converting the image to numpy array and plotting the original image ---------------------
image = f'../../../figs/{image_name}'
if dim==3:
    main = numpy.array(Image.open(image).convert('RGB'), dtype=numpy.float64)
    r, c, ch = main.shape
elif dim==2:
    main = numpy.array(Image.open(image).convert('L'), dtype=numpy.float64)
    r, c = main.shape

kV = [10, 20, 50, 100]                      #these are the k values I am using

plot.figure(figsize=(18, 10))
plot.subplot(1, 5, 1)
if dim == 3:
    plot.imshow(main.astype(numpy.uint8))
elif dim==2:
    plot.imshow(main.astype(numpy.uint8), cmap='gray')
plot.title("Original")
plot.axis("off")

print("")
print("----The following are the frobenius errors for the chosen k values----")
for i, k in enumerate(kV, start=2):
    final = numpy.empty_like(main)
    if dim == 3:
        c_lib.find_svd(main, r, c, k, ch, final)
        e = c_lib.Frobenius(main, final, r, c, ch)
    elif dim==2:
        c_lib.find_svd(main, r, c, k, final)
        e = c_lib.Frobenius(main, final, r, c)
    print(f"Frobenius Error of the image for k = {k} is {e}")                   #printing the frobenius error
    final = numpy.clip(final, 0, 255).astype(numpy.uint8)

    #------------Plotting the compressed images--------------
    plot.subplot(1, 5, i)
    if dim == 3:
        plot.imshow(final)
    elif dim==2:
        plot.imshow(final, cmap='gray')
    plot.title(f"k = {k}")
    plot.axis("off")

#------------Calculating the runtime------------
k_cal = 50
final_c = numpy.empty_like(main)
time_sc = time.perf_counter()
if dim==3:
    c_lib.find_svd(main, r, c, k_cal, ch, final_c)
elif dim==2:
    c_lib.find_svd(main, r, c, k_cal, final_c)
time_ec = time.perf_counter() - time_sc 

final_np = numpy.empty_like(main)
time_snp = time.perf_counter()
if dim==3:
    for p in range(ch):
        A = main[:, :, p]
        U, S, Vt = numpy.linalg.svd(A, full_matrices=False)
        U_k = U[:, :k_cal]
        S_k = numpy.diag(S[:k_cal])
        V_k = Vt[:k_cal, :]
        final_np[:, :, p] = U_k @ S_k @ V_k
elif dim==2:
    U, S, Vt = numpy.linalg.svd(main, full_matrices=False)
    U_k = U[:, :k_cal]
    S_k = numpy.diag(S[:k_cal])
    V_k = Vt[:k_cal, :]
    final_np = U_k @ S_k @ V_k
time_enp = time.perf_counter() - time_snp

#-----------Printing the runtime and errors---------------
mse = numpy.mean((final_c - final_np) ** 2)                 #mean squared error
rel_error = numpy.linalg.norm(final_c - final_np) / numpy.linalg.norm(final_np)             #relative error
print("")
print("----The following are the stats to compare between my algorithm and python standard command----")
print(f"C Power Method runtime: {time_ec:.4f} s")
print(f"NumPy SVD runtime: {time_enp:.4f} s")
print(f"Mean Squared Error: {mse:.6f}")
print(f"Relative Error: {rel_error:.6e}")

plot.tight_layout()
plot.show()
