# HR(S)TEM filter

## Introduction
`hrtem_filter` provides a set of python functions to denoise HR(S)TEM images. A Wiener filter and an average background subtraction filter were designed based on __R. Kilaas J. Microscopy, 1998, 190, 45-51__. Some steps are adopted from D. R. G. Mitchell's script for GMS. A nonlinear filter was adopted by the "non-linear filter plugin from GMS", originally developed by Dr. Hongchu Du. Refer to the original paper: __Hongchu Du, A nonlinear filtering algorithm for denoising HR(S)TEM micrographs, Ultramicroscopy 2015, 151, 62-67__

## Installation
Use pip:

`pip install hrtem_filter`

## Usage
The package contains three main filter functions:

`wiener_filter(img, delta=5, lowpass=True, lowpass_cutoff=0.3, lowpass_order=2)`

`abs_filter(img, delta=5, lowpass=True, lowpass_cutoff=0.3, lowpass_order=2)`

These are for Wiener filter or average background subtraction filter. It takes `img`, an image array as an input, and returns the filtered image array and a difference image array. Parameters are the following:

delta: a threashold for background averaging. Smaller number results in more iterations in refining the averaged background and hense longer time. 

lowpass: also apply a lowpass filter after filtering to remove the residual high frequency noise.

lowpass_cutoff: a cutoff ratio in frequency domain for the lowpass. 1 means no filtering and vice versa.

lowpass_order: order for the Butterworth filter; smaller int retults more tapered cutoff

`nlfilter(img, N, delta=10, lowpass_cutoff=0.3, final_lowpass = False, lowpass_order=2)`

This function carries out the nonlinear filter, a combination of Gaussian lowpass and Wiener filters. It takes `img`, an image array as an input, and returns the filtered image array and a difference image array. Parameters are the following:

N: number of iterations for the lowpass + wiener filtering. More iterations give better noise reduction and takes more time.

final_lowpass: if True, also applys a Butterworth lowpass filter after the nonlinear filtering. This helps to remove residual high frequency noise.

All other parameters are the same as in the `wiener_filter`.

In addition, `hrtem_filter` also provides two lowpass filters: Butterworth filter and Gaussian filter.

`bw_lowpass(img, order, cutoff_ratio)`

`gaussian_lowpass(img, cutoff_ratio)`

## Contact
Send your questions and suggestions to Dr. Tao Ma at matao1984@gmail.com



