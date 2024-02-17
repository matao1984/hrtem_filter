import numpy as np
from scipy.fft import fft2, fftshift, ifftshift, ifft2
import pandas as pd
from scipy.signal import medfilt2d

# Crop the input images if they are not square
def crop_to_square(img):
    x, y = img.shape
    if x > y:
        new_start = int((x - y) / 2)
        new_end = new_start + y 
        img_crop = img[new_start:new_end,:]
    else:
        new_start = int((y - x) / 2)
        new_end = new_start + x 
        img_crop = img[:,new_start:new_end]
    return img_crop



# For radial integration, convert image indices to polar coordinates
def img_to_polar(img):
    # Feed an image array, generate a polar indices array 
    y, x = np.indices(img.shape)
    center = (img.shape[0] + 1) / 2, (img.shape[1] +1 ) / 2
    x = x - center[0]
    y = y - center[1]
    rho = np.hypot(y, x) # calculate sqrt(x**2 + y**2)
    #rho = rho.astype(int) # Convert rho to int for simplicity
    #phi = np.arctan2(y, x) # Don't need this    
    return rho # rho is the radial distance

# Gaussian low pass filter
def gaussian_lowpass(img, cutoff_ratio):
    """
    img: image array to be filtered, must be square
    cutoff_ratio: cutoff ratio in frequency domain
    """
    r = img_to_polar(img)
    
    # Compute the FFT to find the frequency transform
    fshift = fftshift(fft2(img))

    # Calculate the cutoff frequency
    cutoff = img.shape[0] * cutoff_ratio
    
    # Create Gaussian mask
    gaussian_filter = np.exp(- (r**2) / (2 * (cutoff**2)))
    
    # Apply the filter to the frequency domain representation of the image
    filtered_fshift = fshift * gaussian_filter

    # Apply the inverse FFT to return to the spatial domain
    img_glp = ifft2(ifftshift(filtered_fshift))
    img_glp = np.real(img_glp)
    return img_glp

# Butterworth lowpass filter
def bw_lowpass(img, order, cutoff_ratio):
    """
    img: image array to be filtered, must be square
    order: Butterworth order
    cutoff_ratio: cutoff ratio in frequency domain
    """
    r = img_to_polar(img) # Convert to polar indices
    bw = 1/(1+0.414*(r/(cutoff_ratio * img.shape[0]))**(2*order))
    
    # Compute the FFT to find the frequency transform
    fshift = fftshift(fft2(img))
    
    # Apply the filter to the frequency domain representation of the image
    filtered_fshift = fshift * bw

    # Apply the inverse FFT to return to the spatial domain
    img_bw = ifft2(ifftshift(filtered_fshift))
    img_bw = np.real(img_bw)
    return img_bw

    

# Function to get an averaged background from a real-space HR image
def get_avg_background(img, delta=10):
    """
    img: 2D array of real-space HR image data
    delta: a threashold for background averaging
    """
    # Get the polar indices array
    r = img_to_polar(img)
    r = r.astype(int)
    
    # Get a Butterworth filter on image to remove the edge effect
    noedgebw = 1/(1+0.414*(r/(0.4 * r.shape[0]))**(2*12))
    noedgeimg = img * noedgebw
    f_noedge = fftshift(fft2(noedgeimg))
    # Light filter the FFT for processing
    f_img = medfilt2d(np.abs(f_noedge),kernel_size=5)

    # Convert the FFT data into flattened real magnitude data
    f_mag = f_img.flatten()
    # Flatten the r array
    r_mag = r.flatten()

    # Sort the polar indices and separate the FFT data by r 
    df = pd.DataFrame({'r_mag': r_mag, 'f_mag': f_mag})    
    grouped = df.groupby('r_mag')['f_mag'].apply(list)   
    r_bin = np.array(grouped.index)
    f_bin = np.array(grouped.values) 

    # Get a mean plot for the FFT
    f_mean0 = [np.mean(bin) for bin in f_bin]


    # For each r bin, replace the pixels > mean with mean0, and take the new mean1.
    # Compare the mean1 and mean0, if the difference is < threashold%, stop
    # This needs to be done through all the bins
    # Make a copy of the original mean list
    f_mean = f_mean0[:]
    for i in range(len(f_bin)):
        bin = np.array(f_bin[i])
        mean = np.mean(bin)
        # Initialize difference as infinity
        diff_pc = np.inf
        # While the percentage difference between mean values is greater than the threshold
        while diff_pc > delta and mean > 0:
            # Replace any data > mean 
            bin[bin > mean] = mean   
            # Recalculate mean
            new_mean = np.mean(bin)   
            # Calculate the percentage difference
            diff_pc = np.abs((new_mean - mean) / mean) * 100    
            # Update the mean value
            mean = new_mean
        if mean < 0:
            mean = 0
        # Update the mean data
        f_mean[i] = mean

    # Construct the background array
    f_avg_r = dict(zip(r_bin, f_mean))
    f_avg = np.zeros(f_img.shape)
    for i in range(r.shape[0]):
        for j in range(r.shape[1]):
            f_avg[i, j] = f_avg_r[r[i, j]]
                                    
    return f_avg


# Wiener filter function
def wiener_filter(img, delta=5, lowpass=True, lowpass_cutoff=0.3, lowpass_order=2):
    """
    Wiener filter for HRTEM images
    img: the image data array
    delta: a threashold for background averaging
    lowpass: also apply a lowpass filter after filtering
    lowpass_cutoff: a cutoff ratio in frequency domain for the lowpass
    lowpass_order: order for the Butterworth filter; smaller int retults more tapered cutoff
    Return: filtered image array and difference
    """
    f_img = fftshift(fft2(img))
    fu = np.abs(f_img)
    fa = get_avg_background(img, delta=delta)
    wf = (fu**2 - fa**2)/fu**2
    wf[wf<0] = 0
    f_img_wf = f_img * wf
    img_wf = ifft2(ifftshift(f_img_wf))
    img_wf = np.real(img_wf)
    if lowpass:
        img_wf = bw_lowpass(img_wf, lowpass_order, lowpass_cutoff)
    img_wf = np.single(img_wf)
    img_diff = img - img_wf
    return img_wf, img_diff
    
# Average background subtraction filter function
def abs_filter(img, delta=5, lowpass=True, lowpass_cutoff=0.3, lowpass_order=2):
    """
    ABS filter for HRTEM images
    img: the image data array
    delta: a threashold for background averaging
    lowpass: also apply a lowpass filter after filtering
    lowpass_cutoff: a cutoff ratio in frequency domain for the lowpass
    lowpass_order: order for the Butterworth filter; smaller int retults more tapered cutoff
    Return: filtered image array and difference
    """
    f_img = fftshift(fft2(img))
    fu = np.abs(f_img)
    fa = get_avg_background(img, delta=delta)
    absf = (fu - fa)/fu
    absf[absf<0] = 0
    f_img_absf = f_img * absf
    img_absf = ifft2(ifftshift(f_img_absf))
    img_absf = np.real(img_absf)
    if lowpass:
        img_absf = bw_lowpass(img_absf, lowpass_order, lowpass_cutoff)
    img_absf = np.single(img_absf)
    img_diff = img - img_absf
    return img_absf, img_diff


# Nonlinear filter function            
def nlfilter(x_in, N, delta=10, lowpass_cutoff=0.3, final_lowpass = False, lowpass_order=2):
    """
    Non-linear filter
    img: img 2D-array
    N: number of iterations
    lowpass_cutoff: cutoff of the low pass filter
    final_lowpass: apply a Butterworth lowpass filter on the final image
    The Butterworth filter will use lowpass_order and lowpass_cutoff
    Return: filtered image array and difference
    """
    i=0
    while i < N:
        x_lp = gaussian_lowpass(x_in, lowpass_cutoff)
        x_diff = x_in - x_lp
        x_diff_wf, _ = wiener_filter(x_diff,
                                  delta=delta,
                                 lowpass=False)
        x_in = x_lp + x_diff_wf
        i = i+1

    if final_lowpass:
        x_in = bw_lowpass(x_in, lowpass_order, lowpass_cutoff)
    img_filtered = np.single(x_in) # Convert to 32 bit float
    img_diff = x_in - img_filtered
    return img_filtered, img_diff
