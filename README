###############################################
# IRAC photometry for WISP program
# Ivano Baronchelli 2017
###############################################

Requirements
------------

- Python (2.7 or later), and its modules numpy, scipy, astropy, matplotlib
  (get them using:
  $ sudo apt-get install python-numpy
  $ sudo pip install astropy
  and similar commands)

- SExtractor software (2.19.5 or higher) 
  See https://www.astromatic.net/software/sextractor

- SWarp (2.19.1 or higher)
  See: https://www.astromatic.net/software/swarp
  FOR UNIX SYSTEMS:
  - Download the pre-compiled rpm version of the code.
  - Then, in a bash shell run the following commands:
  $ ./configure
  $ sudo make
  $ sudo make install
  For a Linux-Ubuntu systems, an older (but stable) version of SWarp can be obtained from Ubuntu software center.

- python scikit package (http://scikit-image.org/docs/dev/auto_examples/filters/plot_deconvolution.html)

- photutils library (0.3.2)

- CFITSIO library
  (download the tarball from website http://heasarc.nasa.gov/fitsio/fitsio.html
  and follow instructions)

- FFTW3 library
  (download the tarball from website http://www.fftw.org/download.html and follow instructions).


How to run IRAC_photometry.
---------------------------
1) Set the paths in the IRAC_photometry/paths.py file. 
   The paths that need to be changed are highlighted and are:
   self.IRAC_DATA_PATH_I1='set_your_path/' # path to the original IRAC 1 data images
   self.IRAC_DATA_PATH_I2='set_your_path/' # path to the original IRAC 2 data images
   self.IRAC_PHOT_FOLDER='set_your_path/'  # location of IRAC_photometry.py
   self.TPHOT_LOCATION='set_your_path'     # PATH to the TPHOT executable




NOTE on IRAC_photometry.py (how it works)
------------------
############################
IRAC RECENTERING
############################
- IRAC images are recentered to the HST images during the process. This is needed in order to make tphot able to work using the HST source positions and shapes as priors for the IRAC photometry. HST catalogs of sources were created when running the WISP pipeline. The recentering is performed using an iterative procedure, where IRAC positions are determined running SExtractor. The corresponding HST counterparts are searched using an initial searching radius of 2 arcseconds. Then, the RA,DEC positions in the HST final catalog and in the IRAC catalog just created are compared and a shift correction is computed. Finally, the shift correction is applied to the IRAC images. After this initial correction, the entire process is repeated using, this time, a shorter searching radius (0.75 arcseconds). This allows us to use only those sources for which a secure counterpart is found, making the final correction more precise. At hte end of the process, the WCS of IRAC images is recentered to the HST reference (original images are saved in a separate folder called "IRAC_original".

############################
IRAC image SWarping
############################
- TPHOT requires input images with identical (or integer multiple) pixel-scale and identical pixel orientation. In order to make the IRAC images consistent with the HST reference, the "SWarp" software is used. SWarp intermediate products are images with similar characteristics (pixel scale and orientation). This SWarp intermediate products are used (while we do not need SWarp final outputs).

############################
HST SEGMENTATION MAP
############################
- Tphot requires both a catalog of sources extracted in the high-resolution image (HST), and the correspondent segmentation image. The IDs in the input catalog must correspond to the pixel values of the correspondent sources in the segmentation map. For this reason, the segmentation maps created by the WISP pipeline can not be used as they are. While SExtractor assigns IDs and creates segmentation maps with values based on positional bases, in the final wisp catalog, the IDs are re-ordered as a function of the source magnitudes. Additinally, some sources are removed from the catalog, mainly because located in the image borders or for other problematic issues. However, their correspondent counterparts in the segmentation map are still there. Finally, there are some sources detected in the J or in the H band, but not in the final combined (J+H) image. These sources are indicated in the catalog using IDs higher than 10000, 20000 and 30000. These sources do not have a counterpart in the J+H segmentation map, computed when SExtracting the combined J+H image, for the reasons explained. For these reasons, there is not an immediate correspondence between IDs in the catalog and in the segmentation map.
This program (IRAC_photometry.py) creates an appropriate segmentation map called "HST_segmentation.fits". In this segmentation map:
 i) there is an univocal correspondence between the pixel values in the segmentation map and the IDs of the corresponding sources in the catalog; 
 ii) sources that are not included in the catalog (border or othe problems) are removed from the segmentation map;
 iii) sources detected in J or H but not in the J+H image are included in the segmentation map with the correct final value (10000+, 20000+, 30000+). The shape of these sources in the segmentation map is computed from the original segmentation maps created when extracting  J and H images singularly.
Creating this segmentation map requires the complete output of the WISP pipeline. In particular, the "SEX" folder, that is not included in the delivered version of the data, is required.

############################
Tphot input catalog
############################
- While creating the segmentation map, the program also writes an output catalog that is used as input by TPHOT. This catalog includes the IDs (HST), X and Y positions in the HST reference image, X and Y borders for each source (min and max taken from the segmentation map), a background (set to zero, since it is already subtracted from the images before running TPHOT) and the reference Flux (HST J or H). 

############################
Convolution Kernel:
############################
Another input required by tphot is a convolution kernel image, that allows to "transform" the HST PSF to the IRAC one.
 PSF(LRI) = Kernel * PSF(HRI)
with LRI meaning "Low Resolution Image" (IRAC) and HRI "High Resolution Image" (HST J or H).
We compute the convolution kernel deconvolving the IRAC PSF using the HST PSF as deconvolution kernel.
IRAC and HST PSFs are obtained averaging the images of the point-like and not saturated sources in the field.
This approach is preferred to the use of synthetic PSFs because IRAC PSF is not rotationally symmetric (it looks more like a triangle than a circle). The orientation of this triangle is different from field to field and it is also the result of a SWarping (re-drizzling) process on a smaller pixel-scale.


############################
Tphot run and output catalog
############################
Tphot is run on the IRAC swarped images using the HST images as a reference. Tphot is run in two separate passages. 
The tphot output is not delivered as it is, but a series of changes are made by IRAC_photometry.py:
- Fluxes (in counts) are converted to AB magnitudes
- An ID correspondence is created with the HST catalogs (same number of sources, same order.
- The 3-sigma depth is computed and reported for IRAC undetected sources.
- Coordinates in the HST reference and in the ORIGINAL (not recentered) IRAC image are included. 
  These IRAC not-recentered coordinates should be the reference RA and DEC when observing targets in the sky.


############################
MAIN OUTPUTS 
############################
NOTE: main outputs are saved in the folder:
wisp_path/axe/Parxx/IRAC/Parxx_IRAC_V6.2.tar.gz

- CATALOG: fin_I1.cat (fin_I2.cat for IRAC 2). This catalog contains the IRAC photometry. 
    IMPORTANT NOTES: 
      1) the photometry refers to HST detected sources. 
          Outside the HST image borders, photometry can't be computed using tphot. 
      2) IRAC photometry for sources outside (and inside) the HST image borders can be
          taken from the catalogs saved in IRAC_SEx_cat.
	  In particular, "tmp_cat_CH1_it0.fits" should be used to have correct IRAC original coordinates.
	  A copy of this file can be found in the "SEx_complete_cat_CH1.fits", in the main folder.
- HST SEGMENTATION MAP: "HST_segmentation.fits" is a segmentation map where IDs correspond in an univocal way
     to the IDs in the HST and IRAC catalogs.
- RESAMPLED IRAC IMAGES: IRAC resampled "sci" and "rms" images are created. Their names are (for par 68, channel 1):
      "WISPS68_ch1_mosaic.resamp_bck_sub.fits" and "WISPS68_ch1_unc.resamp.fits"
     These images have teh following characteristics:
       1) Same pixel size as HST J and/or H images.
       2) Pixel orientation identical to HST J and/or H images.
       3) Same total size as J and/or H images.
       4) Sci image is background subtracted.
- PSFS: an image of the IRAC and HST aproximated PSF is created ("PSF_HST.fits","PSF_IRAC1.fits"). 
  Note that these PSFs could be noisy for some purposes. They work well for what is needed n tphot.
- CONVOLUTION KERNEL: "Kernel_HST_I1.fits" that can be used, for example, to smooth an HST image making it similar to an IRAC one.
