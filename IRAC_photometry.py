###############################################
# IRAC photometry for WISP program
# Ivano Baronchelli 2017
#----------------------------------------------
# - For instructions, run:
#     python IRAC_photometry.py
# - This program can be run using run_IRAC_photometry.sh:
#     source run_IRAC_photometry.sh 
###############################################

###############################################
###############################################
# IMPORT MODULES
###############################################
###############################################
import os
import astropy
import numpy as np
from numpy import where
from numpy import array
from astropy.io import fits # read fits files
import pdb # Python debugger (to stop and interact with the variables during the run, use pdb.set_trace())
from pdb import set_trace as stop # in this way you can just write stop() to stop the program.
# import matplotlib.pyplot as plt # to plot
import math
import sys # if arguments are needed (example: >python progr.py 'argument'
import matplotlib.pyplot as plt # to plot
from scipy.signal import convolve2d as conv2    # Used for deconvolution purposes
from scipy import stats
from skimage import restoration #, color, data # Used for deconvolution purposes
import shutil # to remove directories using remove_tree

# Functions written by Ivano for general purposes (not strictly related to this program:)
from classes.Ivano_functions_1_1 import mk_regionfile
from classes.Ivano_functions_1_1 import image_shift
from classes.Ivano_functions_1_1 import image_depth
# Classes written by Ivano for general purposes (not strictly related to this program:)
from classes.Ivano_classes_1_1 import Deg2sex
from classes.Ivano_classes_1_1 import Match
from classes.Ivano_classes_1_1 import Match_cat
from classes.Ivano_classes_1_1 import Readcol

from paths import Get_paths
###############################################
###############################################



###############################################
###############################################
# CLASSES
###############################################
###############################################

class Data_def(object):
################################################
# 1) Determine the kind of data contained in path0
# Options for HST : J+H , J only, H only
# For IRAC: IRAC-1,  IRAC2
# 2) Determines the naming convenction of the 
# original IRAC2 image
################################################
    def __init__(self,paths,PARN):
        self.PARN=str(PARN)
        self.J=os.path.isfile(paths.p0+'F110W_sci.fits')
        self.H1=os.path.isfile(paths.p0+'F140W_sci.fits')
        self.H2=os.path.isfile(paths.p0+'F160W_sci.fits')
        self.I1=os.path.isfile(paths.p1+'WISPS'+self.PARN+'_ch1_mosaic.fits')
        # CH2 special bahaviour #########################
        ADD_CH2_STR_tmp='_chan2' # NOTE: some fields (all?) have an additional string on their name
        self.I2A=os.path.isfile(paths.p2+'WISPS'+self.PARN+'_ch2_mosaic.fits')
        self.I2B=os.path.isfile(paths.p2+'WISPS'+self.PARN+ADD_CH2_STR_tmp+'_ch2_mosaic.fits')
        
        self.ADD_CH2_STR=''
        if self.I2A==True:
            self.I2=self.I2A
            self.ADD_CH2_STR=''
        if self.I2B==True:
            self.I2=self.I2B
            self.ADD_CH2_STR=ADD_CH2_STR_tmp
        if (self.I2A!=True) and (self.I2B!=True):
            self.I2=False
            self.ADD_CH2_STR=''

        
    def show(self):
    ################################################
    # Shows the kind of data contained in the folder
    ################################################
        print "####################"
        print " This field has"
        string=''
        if self.J==True:
            string=string+'F110'
        if self.H1==True:
            if self.J==True:
                string=string+' and '
            string=string+'F140'
        if self.H2==True:
            if self.J==True:
                string=string+' and '
            string=string+'F160'
        print string
        print '-> Grism data not tested'
        # print "####################"
        if self.I1==True:
            print 'IRAC-1 IMAGE PRESENT'
        else:
            print 'IRAC-1 IMAGE NOT PRESENT'
        if self.I2==True: print 'IRAC-2 IMAGE PRESENT'
        else:
            print 'IRAC-2 IMAGE NOT PRESENT'
        print "####################"


class Get_Filenames(object):
    #---------------------------------------------------
    # files definitions
    #---------------------------------------------------
    def __init__(self,PARN,DS):
        self.PARN=str(PARN)
        # Suffix for resampled and background subtracted images
        self.resamp_suffix=".resamp"                                      # resampled 
        self.bck_sub_suffix='_bck_sub'                                    # background subtracted
        self.resamp_bck_sub_suffix=self.resamp_suffix+self.bck_sub_suffix # resampled and background subtracted

        #####################################################################################
        ###############################   MAIN OUTPUTS   ####################################
        #####################################################################################
        # The following catalogs are the main outputs of this program. output format is ascii.
        self.IRAC_cat_tphot_1='fin_I1.cat' # IRAC-1 T-phot cat. Recentered high-res. image
        self.IRAC_cat_tphot_2='fin_I2.cat' # IRAC-2 T-phot cat. Recentered high-res. image
        # Sextractor catalogs (fits files). They include all the sources SExtracted in the IRAC area,
        #  also outside the HST image borders
        self.SEx_complete_cat_I1="SEx_complete_cat_I1.fits" # IRAC SExtractor cat (includes sources outside HST area)
        self.SEx_complete_cat_I2="SEx_complete_cat_I2.fits" # IRAC SExtractor cat (includes sources outside HST area)
        # output HST Segmentation map. IDs correspond to the IDs in the WISPs catalogs.
        self.out_seg='HST_segmentation.fits'   # output Segmentation map (input for TPHOT)
        # IRAC resampled (SWarped) and background subtracted images
        self.IRAC_resamp1='WISPS'+self.PARN+'_ch1_mosaic'+self.resamp_bck_sub_suffix+'.fits' # IRAC 1 resampled (SWarped) image
        self.IRAC_resamp2='WISPS'+self.PARN+'_ch2_mosaic'+self.resamp_bck_sub_suffix+'.fits' # IRAC 2 resampled (SWarped) image
        self.IRAC_unc_resamp1='WISPS'+self.PARN+'_ch1_unc'+self.resamp_suffix+'.fits' #  IRAC 1 resampled (SWarped) uncert. img
        self.IRAC_unc_resamp2='WISPS'+self.PARN+'_ch2_unc'+self.resamp_suffix+'.fits' #  IRAC 2 resampled (SWarped) uncert. img
        # Output folder name:
        self.output_folder='Par'+self.PARN+'_IRAC_V6.2'  # in this folder (.tar.gz), the main outputs will be saved
        #####################################################################################
        #####################################################################################
        #####################################################################################
        
        # IRAC input images #################################
        # CH1 
        self.IRAC_image1='WISPS'+self.PARN+'_ch1_mosaic.fits' # IRAC 1 original image name 
        self.IRAC_unc1='WISPS'+self.PARN+'_ch1_unc.fits'      # IRAC 1 original uncertainty image
        self.IRAC_cov1='WISPS'+self.PARN+'_ch1_cov.fits'      # IRAC 1 original coverage image
        # CH2 BAD ORIGINAL NAMES ############################
        # NOTE: a bad original naming convenction (_chan2_ch2_) is changed to  _ch2_" inside the working directory
        self.IRAC_image2_BAD='WISPS'+self.PARN+DS.ADD_CH2_STR+'_ch2_mosaic.fits' # **BAD** IRAC 2 original image 
        self.IRAC_unc2_BAD='WISPS'+self.PARN+DS.ADD_CH2_STR+'_ch2_unc.fits'      # **BAD** IRAC 2 original uncertainty image
        self.IRAC_cov2_BAD='WISPS'+self.PARN+DS.ADD_CH2_STR+'_ch2_cov.fits'      # **BAD** IRAC 2 original coverage image
        # CH2 GOOD NAMING CONVENCTION (similar to CH1) ######
        self.IRAC_image2='WISPS'+self.PARN+'_ch2_mosaic.fits' # IRAC 2 new name into the working directory)
        self.IRAC_unc2='WISPS'+self.PARN+'_ch2_unc.fits'      # IRAC 2 original uncertainty image
        self.IRAC_cov2='WISPS'+self.PARN+'_ch2_cov.fits'      # IRAC 2 original coverage image
        #####################################################

        # IRAC temporary SExtractor catalogs. IRAC images are recentered to HST positions, at this point.
        # The following catalogs are used for internal temporary purposes and are in fits format
        self.IRAC_cat_SEx_LRES_1='cat_CH1_SExtractor_LRES.fits' # IRAC-1 SExtractor cat. Recentered image (low resolution)
        self.IRAC_cat_SEx_LRES_2='cat_CH2_SExtractor_LRES.fits' # IRAC-2 SExtractor cat. Recentered image (low resolution)
        self.IRAC_cat_SEx_HRES_1='cat_CH1_SExtractor_HRES.fits' # IRAC-1 SExtractor cat. Recentered image (high resolution)
        self.IRAC_cat_SEx_HRES_2='cat_CH2_SExtractor_HRES.fits' # IRAC-2 SExtractor cat. Recentered image (high resolution)
        #-------------------------------------------------------------------------------------------------------
        # HST sci images (input) and background subtracted sci images (output) names.
        self.J_image="F110W_sci.fits"
        self.J_image_bck_sub='F110W_sci'+self.bck_sub_suffix+'.fits'
        self.J_bck='F110W_bck.fits'
        if DS.H1==True and DS.H2==False:
            self.H_image="F140W_sci.fits"
            self.H_image_bck_sub='F140W_sci'+self.bck_sub_suffix+'.fits'
            self.H_bck='F140W_bck.fits'
        if DS.H1==False and DS.H2==True:
            self.H_image="F160W_sci.fits"
            self.H_image_bck_sub='F160W_sci'+self.bck_sub_suffix+'.fits'
            self.H_bck='F160W_bck.fits'
        if DS.H1==False and DS.H2==False:
            self.H_image="undefined"
            self.H_image_bck_sub="undefined"
        self.JH_image='JH_combined_sci.fits'
        self.JH_image_bck_sub='JH_combined_sci'+self.bck_sub_suffix+'.fits'
        self.JH_bck='JH_combined_bck.fits'
        #-------------------------------------------------------------------------------------------------------
        # Reference HST catalog, image and segmentation image for HRI (High Resolution Image)
        if DS.J==True and (DS.H1==True or DS.H2==True): 
            self.HST_image=self.JH_image           # HST Reference image
            self.HST_bck_sub=self.JH_image_bck_sub # HST Reference image, background subtracted
            self.HST_catalog='fin_F110.cat'       # HST reference catalog
            self.HST_seg='JH_combined_seg.fits'   # HST reference segmentation
        if DS.J==True and (DS.H1==False and DS.H2==False): 
            self.HST_image='F110W_sci.fits'       # HST Reference image:
            self.HST_bck_sub=self.J_image_bck_sub # HST Reference image, background subtracted
            self.HST_catalog='fin_F110.cat'       # HST reference catalog 
            self.HST_seg='F110W_seg.fits'         # HST reference segmentation
        if DS.J==False and (DS.H1==True):
            self.HST_image='F140W_sci.fits'       # HST Reference image:
            self.HST_bck_sub=self.H_image_bck_sub # HST Reference image, background subtracted
            self.HST_catalog='fin_F140.cat'       # HST reference catalog 
            self.HST_seg='F140W_seg.fits'         # HST reference segmentation
        if DS.J==False and (DS.H2==True): 
            self.HST_image='F160W_sci.fits'       # HST Reference image:
            self.HST_bck_sub=self.H_image_bck_sub # HST Reference image, background subtracted
            self.HST_catalog='fin_F160.cat'       # HST reference catalog 
            self.HST_seg='F160W_seg.fits'         # HST reference segmentation
        #-------------------------------------------------------------------------------------------------------
        # WISP input Segmentation maps names
        self.J_seg='F110W_seg.fits'
        if DS.H1==True: self.H_seg='F140W_seg.fits'
        if DS.H2==True: self.H_seg='F160W_seg.fits'
        self.JH_seg='JH_combined_seg.fits'
        #-------------------------------------------------------------------------------------------------------
        self.HRI_catalog="HRI_catalog.cat"     # output High Resolution Image catalog. This is a TPHOT input.
        #-------------------------------------------------------------------------------------------------------
        self.outname1='SWarped_1'              # Do not include extenction! Name of output image CH-1 
                                               # and corresponding orientation header file (used by SWarp).
        self.outname2='SWarped_2'              # Do not include extenction! Name of output image CH2
                                               #  and corresponding orientation header file (used by SWarp).
        self.headname1=self.outname1+'.head'   # SWarp additional header for orientation - CH1
        self.headname2=self.outname2+'.head'   # SWarp additional header for orientation - CH2
        #-------------------------------------------------------------------------------------------------------
        self.PSF_HST='PSF_HST.fits'            # Average HST PSF computed for sources in this field (F110, 140 or 160)
        self.PSF_IRAC1='PSF_IRAC1.fits'        # Average IRAC-1 PSF computed for sources in this field
        self.PSF_IRAC2='PSF_IRAC2.fits'        # Average IRAC-2 PSF computed for sources in this field
        self.KERNEL_1='Kernel_HST_I1.fits'     # deconvolution Kernel for: PSF(IRAC1)=K*PSF(HST)
        self.KERNEL_2='Kernel_HST_I2.fits'     # deconvolution Kernel for: PSF(IRAC2)=K*PSF(HST)
        #-------------------------------------------------------------------------------------------------------
        # TPHOT output catalogs. These are raw outputs (fluxes in counts, no IRAC original coordinates included)
        self.tphot_outcat_CH1='RAW_tphot_outcat_CH1.cat'
        self.tphot_outcat_CH2='RAW_tphot_outcat_CH2.cat'
        


        
###############################################
###############################################
# FUNCTIONS
###############################################
###############################################

def check_arg():
    # Check if the correct number of arguments are given in input
    if len(sys.argv) <= 1:
        print "-----------------------------------------"
        print "Par number not specified."
        print "Please specify the Par number as in the following examples:"
        print "> python progr.py 322"
        print " OR "
        print "> python IRAC_photometry.py 322 I ## (for interactive mode)"
        print "> python IRAC_photometry.py 322 A ## (for Automatic mode=default)"
        print "> python IRAC_photometry.py 322 T ## (for Test mode1: interactive. Uses test paths)"
        print "> python IRAC_photometry.py 322 P ## (for Test mode2: automatic. Uses test paths)"
        print "This program is stopped here"
        print "-----------------------------------------"
        exit()
    return len(sys.argv)



def clean_old(paths,PARN,TS):

    # Clean outputs of previous runs. 
    #    ONLY FOR TEST: copy original files in the test folder
    PARN=str(PARN)
    print '---------------------------------------------------------------------'
    print 'Directories, created by previous runs, are eliminated:   '
    print '---------------------------------------------------------------------'
    # The directories eliminated here are only those created in the "create_tree() function"
    if os.path.isdir(paths.p3) == True: 
        print "removing "+ paths.p3+' directory'
        shutil.rmtree(paths.p3, ignore_errors=True)
    if os.path.isdir(paths.pwd) == True: 
        print "removing "+ paths.pwd+' directory'
        shutil.rmtree(paths.pwd, ignore_errors=True)
        
    
    if TS ==1:
        print '---------------------------------------------------------------------'
        print 'The next copy/paste commands are run only for the test phase'
        print '---------------------------------------------------------------------'
        if os.path.isdir(paths.path_test) == False: # SAFETY TEST !
            print "You are trying to eliminate original files for which you don't have a backup copy "
            print "Original files will not be eliminated!"
        if os.path.isdir(paths.path_test) == True:
            print '---------------------------------------------------------------------'
            print "TEST-ONLY: removing old fake-original tree"
            print '---------------------------------------------------------------------'
            if os.path.isdir(paths.path_test+'DIRECT_GRISM/') == True: 
                print "removing "+paths.path_test+'DIRECT_GRISM/'+' directory'
                shutil.rmtree(paths.path_test+'DIRECT_GRISM/', ignore_errors=True)

            if os.path.isdir(paths.path_test+'SEX/') == True: 
                print "removing "+paths.path_test+'SEX/'+' directory'
                shutil.rmtree(paths.path_test+'SEX/', ignore_errors=True)

            if os.path.isdir(paths.path_test+'IRAC_images/') == True: 
                print "removing "+paths.path_test+'IRAC_images/'+' directory'
                shutil.rmtree(paths.path_test+'IRAC_images/', ignore_errors=True)
            print '---------------------------------------------------------------------'
            print "TEST-ONLY: creating a new fake-original tree"
            print '---------------------------------------------------------------------'
            print "creating"+paths.path_test+'DIRECT_GRISM/' +' directory'
            os.system('mkdir '+paths.path_test+'DIRECT_GRISM/' )
            print "creating"+paths.path_test+'SEX/' +' directory'
            os.system('mkdir '+paths.path_test+'SEX/' )
            print "creating"+paths.path_test+'IRAC_images/' +' directory'
            os.system('mkdir '+paths.path_test+'IRAC_images/' )
            print "creating"+paths.path_test+'IRAC_images/CH1/' +' directory'
            os.system('mkdir '+paths.path_test+'IRAC_images/CH1' )
            print "creating"+paths.path_test+'IRAC_images/CH2/' +' directory'
            os.system('mkdir '+paths.path_test+'IRAC_images/CH2/' )
            print '---------------------------------------------------------------------'
            print 'Copying the files in the test-original directories'
            print '---------------------------------------------------------------------'
            print "copying "+paths.path_test +'ORIGINAL_DATA/*sci* and *wht* images to '+paths.p0
            print " AND "+paths.path_test +'ORIGINAL_DATA/*seg* images to '+paths.p5
#            if DS.J==True: 
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F110W_sci* '+paths.p0)
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F110W_wht* '+paths.p0)
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F110W_*_seg* '+paths.p5)
#            if DS.H1==True:
#                os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F140W_sci* '+paths.p0)
#                os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F140W_wht* '+paths.p0)
#                os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F140W_*_seg* '+paths.p5)
#            if DS.H2==True:
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F160W_sci* '+paths.p0)
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F160W_wht* '+paths.p0)
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/F160W_*_seg* '+paths.p5)
            
            print "copying "+paths.path_test +'ORIGINAL_DATA/fin_* to '+paths.p0
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/fin_* '+paths.p0)
#            if DS.J==True and (DS.H1==True or DS.H2==True):
            print "copying "+paths.path_test +'ORIGINAL_DATA/JH_combined_sci to '+paths.p0
            print " AND "+paths.path_test +'ORIGINAL_DATA/JH_combined_seg to '+paths.p5
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/JH_combined_sci* '+paths.p0)
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/JH_combined_seg* '+paths.p5)
            print "copying "+paths.path_test +'ORIGINAL_DATA/WISPS'+PARN+'_ch1* to '+paths.p1
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/WISPS'+PARN+'_ch1* '+paths.p1)
            # IRAC 2 : two possiblities:
            print "copying "+paths.path_test +'ORIGINAL_DATA/WISPS'+PARN+'_ch2* to '+paths.p2
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/WISPS'+PARN+'_ch2* '+paths.p2)
            print "copying "+paths.path_test +'ORIGINAL_DATA/WISPS'+PARN+'_chan2* to '+paths.p2
            os.system('cp '+paths.path_test  +'ORIGINAL_DATA/WISPS'+PARN+'_chan2* '+paths.p2)

def create_tree(paths):
    # creates directory three
    os.system('mkdir '+ paths.pwd) # Create working directory
    os.system('mkdir '+ paths.p3)  # Directory for IRAC catalogs

def copy_files(paths,files,DS):
    #---------------------------------------------------
    # Copy original data into the working directory
    #---------------------------------------------------
    print '--------------------------------------------------------------'
    print ' Copying original files into the '+paths.pwd +' directory'
    print ' For IRAC, a backup folder (named "original") is created'
    print ' in the working directory. The HST reference image is copied too.'
    print '--------------------------------------------------------------'
    # IRAC IMAGES
    print 'Copying IRAC images and uncertainties'
    if DS.I1==True:
        os.system('cp '+paths.p1+files.IRAC_image1+' '+paths.pwd+files.IRAC_image1)
        os.system('cp '+paths.p1+files.IRAC_unc1+' '+paths.pwd+files.IRAC_unc1)
        os.system('cp '+paths.p1+files.IRAC_cov1+' '+paths.pwd+files.IRAC_cov1)
        os.system('mkdir '+ paths.pwd+'IRAC_original/')
        os.system('cp '+paths.p1+files.IRAC_image1+' '+paths.pwd+'IRAC_original/'+files.IRAC_image1)
        os.system('cp '+paths.p1+files.IRAC_unc1+' '+paths.pwd+'IRAC_original/'+files.IRAC_unc1)
        os.system('cp '+paths.p1+files.IRAC_cov1+' '+paths.pwd+'IRAC_original/'+files.IRAC_cov1)
    if DS.I2==True:
        os.system('cp '+paths.p2+files.IRAC_image2_BAD+' '+paths.pwd+files.IRAC_image2)
        os.system('cp '+paths.p2+files.IRAC_unc2_BAD+' '+paths.pwd+files.IRAC_unc2)
        os.system('cp '+paths.p2+files.IRAC_cov2_BAD+' '+paths.pwd+files.IRAC_cov2)
        I_ORIG=os.path.isdir(paths.pwd+'IRAC_original/')
        if I_ORIG==False: os.system('mkdir '+ paths.pwd+'IRAC_original/')
        os.system('cp '+paths.p2+files.IRAC_image2_BAD+' '+paths.pwd+'IRAC_original/'+files.IRAC_image2_BAD)
        os.system('cp '+paths.p2+files.IRAC_unc2_BAD+' '+paths.pwd+'IRAC_original/'+files.IRAC_unc2_BAD)
        os.system('cp '+paths.p2+files.IRAC_cov2_BAD+' '+paths.pwd+'IRAC_original/'+files.IRAC_cov2_BAD)
    # HST reference image (F110, F40 or F160). Needed for SWarp
    print 'Copying HST reference image.'
    os.system('cp '+paths.p0+files.HST_image+' '+paths.pwd)
    print 'Copying HST reference catalog.'
    os.system('cp '+paths.p0+files.HST_catalog+' '+paths.pwd)
    # Copying HST segmentation maps
    print 'Copying HST images and segmentation maps.'
    if DS.J==True: 
        os.system('cp '+paths.p5+'F110W_*_seg* '+paths.pwd+files.J_seg)
        os.system('cp '+paths.p0+'F110W_sci.fits '+paths.pwd+files.J_image)
    if DS.H1==True:
        os.system('cp '+paths.p5+'F140W_*_seg* '+paths.pwd+files.H_seg)
        os.system('cp '+paths.p0+'F140W_sci.fits '+paths.pwd+files.H_image)
    if DS.H2==True:
        os.system('cp '+paths.p5+'F160W_*_seg* '+paths.pwd+files.H_seg)
        os.system('cp '+paths.p0+'F160W_sci.fits '+paths.pwd+files.H_image)
    if DS.J==True and (DS.H1==True or DS.H2==True):
        os.system('cp '+paths.p5+files.JH_seg+' '+paths.pwd+files.JH_seg)
    
    ## # Copying SWarp configuration files to working directory
    ## print 'Copying SWarp configuration files (needed in the same location where SWarp runs)'
    ## os.system('mkdir '+ paths.pwd+paths.p6)
    ## os.system('cp '+paths.p6+'* '+paths.pwd+'SWarp_config/')
    os.system('cp '+paths.p6+'config_1.swarp '+paths.pwd)
    
    # Copying TPHOT configuration files to working directory
    print 'Copying TPHOT configuration files to working directory'
    os.system('cp '+paths.p8+'*.param '+paths.pwd)


def errors(files,err):
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    print "ERROR: "
    print "SWarp needs an header file defining the orientation"
    print "of the resampled IRAC image."
    if err==1: print files.headname1 + '  NOT FOUND !'
    if err==1: print files.headname2 + '  NOT FOUND !'
    print "The program is stopped here"
    print "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
    exit()


def SExtractor_gain(CHANNEL,coverage_img):
    # Determine SExtractor gain using the IRAC parameters (WARM PHASE)
    # and using the coverage image to compute the average coverage
    # CHANNEL = 'CH1' or 'CH2'. IRAC Channel considered
    # coverage_img =filename of the IRAC coverage image for the input channel
    #
    #### GAIN 
    # For reference, see SExtractor manual and, for example:
    # Timlin+2016 (equation 9 in the arxive version)
    #
    # SExtractor GAIN must be set as follows. Given:
    # g=instrumental gain [ADU/e-]
    # T= exposure time single passage (coverage=1)
    # N=Average coverage along the field
    # K=conversion factor from ADU/s --> MJy/str (map units)
    # Then GAIN= (gxTxN)/K
    #---------------------------
    # g_CH1=3.7 [e-/adu] (wharm phase)
    # g_CH2=3.71 [e-/adu]
    # K_CH1=0.1257 [MJy/str]/[ADU/s] from Flux Conversions table (column FLUXCONV) at
    #             http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/warmimgcharacteristics/
    # K_CH2=0.1447 [MJy/str]/[ADU/s] 
    # N_CH1=15 Average coverage in the test image
    # T=30 [s/Coverage]
    # GAIN for ch1 ~ 3.7 x 30 x 15 / 0.1257= 13245 [e-/(MJy/str)]
    # GAIN for ch2 ~ 3.71 x 30 x 15 / 0.1447= 11537 [e-/(MJy/str)]
    ####
    expt=30 # exposure time for a single passage [s/Coverage]
    if str(CHANNEL) == 'CH1':
        instr_gain=3.7 # instrumental gain
        KK=0.1257 #conversion factorfrom ADU/s --> MJy/str =map units
    if str(CHANNEL) == 'CH2':
        instr_gain=3.71 # instrumental gain
        KK=0.1447 #conversion factorfrom ADU/s --> MJy/str =map units
    #-----------------------------------------------------------------------------------------
    # Read Coverage-image to compute average coverage and appropriate GAIN for SExtractor
    fits_file=fits.open(coverage_img)
    fits_file.info()
    img_cov = fits_file[0].data # Image data are stored in this np array
    fits_file.close()
    #-------------------------------------
    # Mean coverage. NOTE: coverages = 1 are not considered (borders)
    Average_cov=np.nanmean(img_cov[img_cov>1]) # mean coverage (NaN excluded) 
    #-------------------------------------
    gain=(instr_gain*expt*Average_cov)/KK
    return gain


def SExtract(CHANNEL,catname,files,paths, PIXSCALE):
    # CHANNEL= 'CH1' or 'CH2'. IRAC Channel considered
    # catname=output catalog name
    # files=filenames
    # paths=path names
    # PIXSCALE= Pixel scale: 0.6 for original image, 0.08 for SWarped image
    ######## MAG_ZEROPOINT NOTES 
    # A) from IRAC handbook:
    # http://irsa.ipac.caltech.edu/data/SPITZER/docs/irac/iracinstrumenthandbook/24/
    # Zeropoint CH1= 18.80+(ch1_vega_to_AB_conversion)= 18.80+2,79=21.59
    # Zeropoint CH2= 18.32+(ch2_vega_to_AB_conversion)= 18.32+3.26=21,58 (identical to CH1, as expected)
    # B) My method to compute zeropoint magnitude (just to be sure)
    # 1) Handbook conversion factor (images in [MJy/str]): 
    #    Flux[Jy]/pixel = (8.461595x10^-6) Flux[MJy/str]/pixel --> 
    #    --> delta_mag=-2.5log(8.461595x10^-6)=-2.5x(-5,07255)=12,68137
    # 2) Given that for tdefinition AB=-2.5log(Flux[uJy])+8.9 -->
    # 3) AB=-2.5log(Flux[MJy/str])+12,68137+8.9 =-2.5log(Flux[MJy/str]) + 21,58
    # --> Since the first term is computed by SExtractor, the second one must 
    # be the zeropoint: 21,58
    ######## GAIN NOTES -- > see SExtractor_gain function

    PIXSCALE=float(PIXSCALE)
    pxs=str(PIXSCALE)
    
    REF_DETECT_THRESH=1.0
    REF_BACKSIZE=32.0
    REF_BACKTHICK=20.0
    
    if PIXSCALE==0.6: # STANDARD CASE (original IRAC IMAGE):
        print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        print 'SExtractor applied to IRAC images  with original resolution'
        print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        I1_img=paths.pwd+files.IRAC_image1
        I1_unc=paths.pwd+files.IRAC_unc1
        I2_img=paths.pwd+files.IRAC_image2
        I2_unc=paths.pwd+files.IRAC_unc2
        filt="gauss_2.5_5x5.conv"
        backsize=str(REF_BACKSIZE)     #"32"
        backthick=str(REF_BACKTHICK)   #"20"
        DETECT_THRESH=str(REF_DETECT_THRESH)   #"1.0"
        ANALYSIS_THRESH=str(REF_DETECT_THRESH) #"1.0"
        gain_mult=1.
        
    if PIXSCALE==0.08: # SWarped IRAC image CASE:
        print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        print 'SExtractor applied to SWarped IRAC images (increased resolution)'
        print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        I1_img=paths.pwd+files.IRAC_resamp1
        I1_unc=paths.pwd+files.IRAC_unc_resamp1
        I2_img=paths.pwd+files.IRAC_resamp2
        I2_unc=paths.pwd+files.IRAC_unc_resamp2
        filt="tophat_5.0_5x5.conv"
        backsize=str(int(round(REF_BACKSIZE*0.6/0.08)))
        backthick=str(int(round(REF_BACKTHICK*0.6/0.08)))
        DETECT_THRESH=str(REF_DETECT_THRESH*((0.08/0.6)**0.5))  #"0.365" # ("2.75" #using a resampled not-scaled unc-map)
        ANALYSIS_THRESH=str(REF_DETECT_THRESH*((0.08/0.6)**0.5))#"0.365" # ("2.75" #using a resampled not-scaled unc-map)
        gain_mult=1. # (0.6/0.08)**2 #using a resampled not-scaled unc-map (however, no difference when using an unc. map!)

    if CHANNEL == 'CH1':
        # Compute actual gain for CH1 image
        # ---------------------------------
        gain_CH1=gain_mult*SExtractor_gain(CHANNEL,paths.pwd+files.IRAC_cov1)
        print 'SExtractor Gain adopted for CH1 channel: '+str(gain_CH1)
        # ---------------------------------
        # SExtract image CH1
        #stop()
        os.system('sex '+I1_img+' -c '+paths.p4+'config_A.txt -CATALOG_NAME '+catname+' -PARAMETERS_NAME '+paths.p4+'default.param -FILTER_NAME '+paths.p4+filt+' -STARNNW_NAME '+paths.p4+'default.nnw -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+I1_unc+' -PIXEL_SCALE '+pxs+' -SEEING_FWHM 1.66 -DETECT_THRESH '+DETECT_THRESH+' -ANALYSIS_THRESH '+ANALYSIS_THRESH+' -DEBLEND_MINCONT 0.0001 -BACK_FILTERSIZE 3 -BACK_SIZE '+backsize+' -BACKPHOTO_THICK '+backthick+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+paths.pwd+'CH1_test.fits -MAG_ZEROPOINT 21.58 -GAIN '+str(gain_CH1))


    if CHANNEL == 'CH2':
        # Compute actual gain for CH2 image
        # ---------------------------------
        gain_CH2=gain_mult*SExtractor_gain(CHANNEL,paths.pwd+files.IRAC_cov2)
        print 'SExtractor Gain adopted for CH2 channel: '+str(gain_CH2)
        # ---------------------------------
        # SExtract image CH2
        os.system('sex '+I2_img+' -c '+paths.p4+'config_A.txt -CATALOG_NAME '+catname+' -PARAMETERS_NAME '+paths.p4+'default.param -FILTER_NAME '+paths.p4+filt+' -STARNNW_NAME '+paths.p4+'default.nnw -WEIGHT_TYPE MAP_RMS -WEIGHT_IMAGE '+I2_unc+' -PIXEL_SCALE '+pxs+' -SEEING_FWHM 1.66 -DETECT_THRESH '+DETECT_THRESH+' -ANALYSIS_THRESH '+ANALYSIS_THRESH+' -DEBLEND_MINCONT 0.0001 -BACK_FILTERSIZE 3 -BACK_SIZE '+backsize+' -BACKPHOTO_THICK '+backthick+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+paths.pwd+'CH2_test.fits -MAG_ZEROPOINT 21.58 -GAIN '+str(gain_CH2) )


def size_orient_from_ref(ref_image,header_CH1,header_CH2):
    # Creates an header file containing the parameters defining
    # the size and orientation of the reference image. This header
    # can be used by SWarp as input file to set the orientation 
    # of the outputs.
    # It returns a string  array containing the values of center 
    # and size of the reference image.
    # Note: two header files are created: one for IRAC channel 1
    #       and one for channel 2. This is required because the 
    #       name of the headerfile must be consistent with the 
    #       image name.
    
    # A) read ref Image (HST) 
    HST_img = fits.open(ref_image)
    HST_DEC_center=HST_img[0].header['CRVAL2']
    HST_RA_center=HST_img[0].header['CRVAL1']
    # A2) DEG --> SEXAGESIMAL coordinates Conversion 
    #            (input type required by SWarp)
    HST_coord=Deg2sex(HST_RA_center,HST_DEC_center)
    string_center=HST_coord.RA(1,2,ndec=2)+' '+HST_coord.DEC(1,2,ndec=2)
    # B) read HST Image size
    HST_X_size=HST_img[0].header['NAXIS1']
    HST_Y_size=HST_img[0].header['NAXIS2']
    string_size=str(HST_X_size)+' '+str(HST_Y_size)
    # C) read HST orientation
    CD1_1=HST_img[0].header['CD1_1']
    CD1_2=HST_img[0].header['CD1_2']
    CD2_1=HST_img[0].header['CD2_1']
    CD2_2=HST_img[0].header['CD2_2']
    # C2) Write HST image orientation in a header file (usde by SWarp)
    # The header MUST be saved where SWarp will be run!
    # CH-1
    f1 = open(header_CH1, 'w')
    f1.writelines('CD1_1   = '+str(CD1_1)+'\n')
    f1.writelines('CD1_2   = '+str(CD1_2)+'\n')
    f1.writelines('CD2_1   = '+str(CD2_1)+'\n')
    f1.writelines('CD2_2   = '+str(CD2_2)+'\n')
    f1.writelines('END'+'\n')
    f1.close()
    # CH-2
    f2 = open(header_CH2, 'w')
    f2.writelines('CD1_1   = '+str(CD1_1)+'\n')
    f2.writelines('CD1_2   = '+str(CD1_2)+'\n')
    f2.writelines('CD2_1   = '+str(CD2_1)+'\n')
    f2.writelines('CD2_2   = '+str(CD2_2)+'\n')
    f2.writelines('END'+'\n')
    f2.close()
    
    HST_img.close()
    return (string_center,string_size)


def resample(CHANNEL,files,paths,string_center,string_size):

    # Reasmple IRAC images to the same resolution of the HST reference.
    # SWarp is run in the working directory 
    #
    pix_AREA_ratio=(0.6/0.08)**2
    # SWarp on IRAC "sci" and "unc" maps ######### and convert units from [MJy/str] to [uJy/pixel]
    # NOTE:
    # 1) files.outname1 is not created. We are interested in the intermediate products. THE NAME OF THIS FILE, HOWEVER, IS NECESSARY: the orientation must be specified in a file "filename.head" where filename must have the same name as the output image
    # 2) output uncertainty maps are multiplied for the square root of the pixel areal ratio. 
    if CHANNEL == 'CH1':
        # sci IMAGE 
        os.system('SWarp '+files.HST_image+' '+files.IRAC_image1+' -c config_1.swarp -IMAGEOUT_NAME '+files.outname1+'.fits -PIXELSCALE_TYPE MIN -WRITE_XML N -DELETE_TMPFILES N -CELESTIAL_TYPE NATIVE -COMBINE N -CENTER_TYPE MANUAL -CENTER "'+string_center+'" -IMAGE_SIZE "'+string_size+'" -RESAMPLE_DIR ./ -RESAMPLE_SUFFIX '+files.resamp_bck_sub_suffix+'.fits -SUBTRACT_BACK Y -BACK_TYPE AUTO -BACK_SIZE 32 -BACK_FILTERSIZE 3')
        # -----------------------------------------------------------------------------
        # unc (rms) IMAGE (no bck subtraction + multiply for [pixel_scale ratio]^2 )
        os.system('SWarp '+files.HST_image+' '+files.IRAC_unc1+' -c config_1.swarp -IMAGEOUT_NAME '+files.outname1+'.fits -PIXELSCALE_TYPE MIN -WRITE_XML N -DELETE_TMPFILES N -CELESTIAL_TYPE NATIVE -COMBINE N -CENTER_TYPE MANUAL -CENTER "'+string_center+'" -IMAGE_SIZE "'+string_size+'" -RESAMPLE_DIR ./ -RESAMPLE_SUFFIX '+files.resamp_suffix+'.fits -SUBTRACT_BACK N')
        # Load image just created
        I1_unc_resamp = fits.open(files.IRAC_unc_resamp1)
        # Multiply for (pixel ratio)^2 : needed to get the same SExtractor output uncertainties using High or low-res. imgs 
        I1_unc_resamp[0].data=I1_unc_resamp[0].data*(pix_AREA_ratio**0.5)
        # -----------------------------------------------------------------------------
        # Set negative values (created by the resampling) to the highest possible value
        # This passage is needed by tphot
        neg=np.where(I1_unc_resamp[0].data <= 0.0)
        neg=np.array(neg)
        YYY=neg[0]
        XXX=neg[1]
        I1_unc_resamp[0].data[YYY,XXX]=np.max(I1_unc_resamp[0].data)
        # -----------------------------------------------------------------------------
        # Remove saved unc-image
        os.system('rm '+files.IRAC_unc_resamp1)
        # Save new unc-image with same name
        I1_unc_resamp.writeto(files.IRAC_unc_resamp1)
        # close
        I1_unc_resamp.close()

    if CHANNEL == 'CH2':
        # sci IMAGE 
        os.system('SWarp '+files.HST_image+' '+files.IRAC_image2+' -c config_1.swarp -IMAGEOUT_NAME '+files.outname2+'.fits -PIXELSCALE_TYPE MIN -WRITE_XML N -DELETE_TMPFILES N -CELESTIAL_TYPE NATIVE -COMBINE N -CENTER_TYPE MANUAL -CENTER "'+string_center+'" -IMAGE_SIZE "'+string_size+'" -RESAMPLE_DIR ./ -RESAMPLE_SUFFIX '+files.resamp_bck_sub_suffix+'.fits -SUBTRACT_BACK Y -BACK_TYPE AUTO -BACK_SIZE 32 -BACK_FILTERSIZE 3')
        # -----------------------------------------------------------------------------
        # unc (rms) IMAGE (no bck subtraction)
        os.system('SWarp '+files.HST_image+' '+files.IRAC_unc2+' -c config_1.swarp -IMAGEOUT_NAME '+files.outname2+'.fits -PIXELSCALE_TYPE MIN -WRITE_XML N -DELETE_TMPFILES N -CELESTIAL_TYPE NATIVE -COMBINE N -CENTER_TYPE MANUAL -CENTER "'+string_center+'" -IMAGE_SIZE "'+string_size+'" -RESAMPLE_DIR ./ -RESAMPLE_SUFFIX '+files.resamp_suffix+'.fits -SUBTRACT_BACK N')
        # Load image just created
        I2_unc_resamp = fits.open(files.IRAC_unc_resamp2)
        # Multiply for (pixel ratio)^1/2 : needed to get the same SExtractor output uncertainties using High or low-res. imgs 
        I2_unc_resamp[0].data=I2_unc_resamp[0].data*(pix_AREA_ratio**0.5)
        # -----------------------------------------------------------------------------
        # Set negative values (created by the resampling) to the highest possible value
        # This passage is needed by tphot
        neg=np.where(I2_unc_resamp[0].data <= 0.0)
        neg=np.array(neg)
        YYY=neg[0]
        XXX=neg[1]
        I2_unc_resamp[0].data[YYY,XXX]=np.max(I2_unc_resamp[0].data)
        # -----------------------------------------------------------------------------
        # Remove saved unc-image
        os.system('rm '+files.IRAC_unc_resamp2)
        # Save new unc-image with same name
        I2_unc_resamp.writeto(files.IRAC_unc_resamp2)
        # close
        I2_unc_resamp.close()
        
def clean_after(DS,paths,files,PARN):
    # CLean resampled images not used in the successive steps
    if (DS.J == True and (DS.H1== True or DS.H2== True)): 
        os.system('rm '+paths.pwd+'JH_combined_sci'+files.resamp_suffix+'*') # remove resampled HST image (not used)
    if (DS.J == True and (DS.H1== False and DS.H2== False)): 
        os.system('rm '+paths.pwd+'F110W_sci'+files.resamp_suffix+'*' ) # remove resampled HST image (not used)
    if (DS.J == False and DS.H1== True): 
        os.system('rm '+paths.pwd+'F140W_sci'+files.resamp_suffix+'*' ) # remove resampled HST image (not used)
    if (DS.J == False and DS.H2== True): 
        os.system('rm '+paths.pwd+'F160W_sci'+files.resamp_suffix+'*' ) # remove resampled HST image (not used)
    # Remove IRAC recentered imgs with original pixel-scale 
    #   (originals are already saved in a different folder) 
    if (DS.I1 == True): 
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch1_mosaic.fits' ) # remove IRAC-sci image with original pixel-scale
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch1_mosaic'+files.resamp_bck_sub_suffix+'.weight.fits' ) # remove IRAC-sci wht image made by swarp (wrong and not used)
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch1_unc.fits' ) # remove IRAC-unc image with original pixel-scale
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch1_unc'+files.resamp_suffix+'.weight.fits' ) # remove IRAC-unc wht image made by swarp (wrong and not used)
    if (DS.I2 == True): 
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch2_mosaic.fits' ) # remove IRAC-sci image with original pixel-scale
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch2_mosaic'+files.resamp_bck_sub_suffix+'.weight.fits' ) # remove IRAC-sci wht img made by swarp (wrong and not used)
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch2_unc.fits' ) # remove IRAC-unc image with original pixel-scale
        os.system('rm '+paths.pwd+'WISPS'+PARN+'_ch2_unc'+files.resamp_suffix+'.weight.fits' ) # remove IRAC-unc wht img made by swarp (wrong and not used)


def background_HST(DS,paths,files):
    # Creates HST background subtracted images
    # and the correspondent background images.
    # ------------------------------------------
    catname=paths.pwd+'temporary_cat.fits' # Temporary catalog. Removed at the end of the process!
    # NOTE: many of the SExtractor parameters are not important at all, in this context.
    # Important paramenters:
    backsize="45"
    backthick="20"
    back_filtersize="3"
    # ----------------- 
    print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    print 'SExtractor run - ONLY FOR BACKGROUND SUBTRACTION-'
    print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
    
    doit='no'
    cy=0
    while cy<3:
        if cy==0 and DS.J==True: 
            imagename=paths.pwd+files.J_image
            Background_image=paths.pwd+files.J_bck
            Bck_subtracted_image=paths.pwd+files.J_image_bck_sub
            doit='yes'
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
            print ' F110 image'
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        if cy==1 and DS.H1==True:
            imagename=paths.pwd+files.H_image
            Background_image=paths.pwd+files.H_bck
            Bck_subtracted_image=paths.pwd+files.H_image_bck_sub
            doit='yes'
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
            print ' F140 image '
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        if cy==1 and DS.H2==True:
            imagename=paths.pwd+files.H_image
            Background_image=paths.pwd+files.H_bck
            Bck_subtracted_image=paths.pwd+files.H_image_bck_sub
            doit='yes'
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
            print ' F160 image'
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        if cy==2 and DS.J==True and (DS.H1==True or DS.H2==True):
            imagename=paths.pwd+files.HST_image
            Background_image=paths.pwd+files.JH_bck
            Bck_subtracted_image=paths.pwd+files.JH_image_bck_sub
            doit='yes'
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
            print ' J & H combined image '
            print 'xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx'
        
        # -----------------     
        if doit=='yes':
            os.system('sex '+imagename+' -c '+paths.p4+'config_HST.sex -CATALOG_NAME '+catname+' -PARAMETERS_NAME '+paths.p4+'default.param -FILTER_NAME '+paths.p4+'gauss_2.0_5x5.conv -STARNNW_NAME '+paths.p4+'default.nnw -WEIGHT_TYPE NONE -PIXEL_SCALE 0.08 -DETECT_THRESH 2.0 -ANALYSIS_THRESH 2.0 -DEBLEND_MINCONT 0.005 -BACK_FILTERSIZE '+back_filtersize+' -BACKPHOTO_TYPE LOCAL -BACK_SIZE '+backsize+' -BACKPHOTO_THICK '+backthick+' -CHECKIMAGE_TYPE BACKGROUND -CHECKIMAGE_NAME '+Background_image+' -MAG_ZEROPOINT 0 -GAIN 1')
            
            # Remove catalog - not used -
            os.system('rm '+catname)
            # load sci and bck images, create a bck subtracted image
            bck=fits.open(Background_image)
            img=fits.open(imagename)
            img_bck_sub=fits.open(imagename)
            # find bck border pixels
            border=np.where(img[0].data == 0)
            border=np.array(border)
            YYY=border[0]
            XXX=border[1]
            # set bck border to proper value (zero)
            bck[0].data[YYY,XXX]=0
            # overwrite 
            bck.writeto(Background_image,clobber='True')
            # Compute img-bck
            img_bck_sub[0].data=img[0].data[:]-bck[0].data[:]
            # Write to a file
            img_bck_sub.writeto(Bck_subtracted_image)
            # Close images
            bck.close()
            img.close()
            
        # -----------------     
        doit='no'
        cy=cy+1

def Segment_and_cat(paths,files,HST_ID,X_IMAGE,Y_IMAGE,FLUX,INTERACT):
    # This function creates a segmentation map and a catalog to be used
    # as tphot input (as High-resolution reference). 
    # Catalog and segmentation map are created in the same cycles 
    # to save computational time. 
    #
    # NOTE: 
    # the number in the detection map DOES NOT CORRESPOND to the 
    # catalog ID Number. Originally the ID is assigned by SExtractor
    # in an image-positional order. Instead, in a successive phase 
    # of the WISP data reduction process, the ID is ordered as a 
    # function of the magnitude, with brighter sources having lower IDs.
    # Moreover, some sources in the border are eliminated from the catalog
    # For this reason, a segmentation map consistent with the catalog 
    # must be created.
    #
    #
    # A) Creates arrays needed for output catalog (input catalog for TPHOT) 
    LOC_BCK=np.zeros(len(HST_ID)) # FLUX*0.
    X_MIN=np.zeros(len(HST_ID)) # FLUX*0.
    X_MAX=np.zeros(len(HST_ID)) # FLUX*0.
    Y_MIN=np.zeros(len(HST_ID)) # FLUX*0.
    Y_MAX=np.zeros(len(HST_ID)) # FLUX*0.
    # B) Creates an array used as a base for the output segmentation map:
    HST_seg_base=fits.open(paths.pwd+files.HST_seg)
    HST_seg_base[0].data[:]=0
    OP1=0
    OP2=0
    OP3=0
    # Cycle on IDs to make catalog and set correspondent value in the segmentation map
    iii=0
    while iii < len(HST_ID):
        if HST_ID[iii]<10000 and OP1==0:
            print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
            print 'Opening reference segmentation map'
            print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
            HST_seg=fits.open(paths.pwd+files.HST_seg) # J or H or JH
            HST_seg[0].data=-HST_seg[0].data
            if INTERACT == 'I': 
                plt.imshow(HST_seg[0].data, cmap='gray')
                plt.title('Initial Segmentation map')
                plt.colorbar()
                plt.show()
            OP1=1
        if HST_ID[iii]>=10000 and HST_ID[iii]<20000 and OP2==0:
            print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
            print 'Opening J segmentation map' 
            print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
            # If IDs like this does exist, it means that both the J and H image are present
            HST_seg=fits.open(paths.pwd+files.J_seg)
            HST_seg[0].data=-HST_seg[0].data
            OP2=1
        if HST_ID[iii]>=20000 and HST_ID[iii]<30000 and OP3==0:
            print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
            print 'Opening H segmentation map' 
            print 'HHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHHH'
            # If IDs like this does exist, it means that both the J and H image are present
            HST_seg=fits.open(paths.pwd+files.H_seg)
            HST_seg[0].data=-HST_seg[0].data
            OP3=1

        # IDX assigned in the original segmentation map
        SEG_VAL=HST_seg[0].data[int(round(Y_IMAGE[iii]-1)),int(round(X_IMAGE[iii]-1))]
        # IF NO ASSOCIATION IS FOUND, AN AUTOMATIC VALUE IS ASSIGNED 
        # to 5 pixels sorrounding the source center. This should never happen,
        # but it is set to prevent possible crashes. A Warning is issued in this case.
        
        
        if SEG_VAL==0:
            # Unexpected behaviour: --------------------------------------------------
            # This should not happen!
            print 'WARNING: ID '+str(HST_ID[iii])+' x='+str(X_IMAGE[iii])+' y='+str(Y_IMAGE[iii])
            print '         This program did not find a correspondence '
            print '         in the segmentation map for ID '+str(HST_ID[iii])
            print '         This is probably due to sources with weird shapes in' 
            print '         the segmentation map. The central region of the source'
            print '         is probably considered outside the source borders.'
            print '         The value will be found among the closest pixels.'
            print '         If no valid value is found, a 5 pixels source is'
            print '         automatically created in the segmentation map, in this position...'
            values1=HST_seg[0].data[int(round(Y_IMAGE[iii]-3)):int(round(Y_IMAGE[iii]+2)),int(round(X_IMAGE[iii]-3)):int(round(X_IMAGE[iii]+2))]
            if len(values1[values1<0])>0:
                MODE_values=stats.mode(values1[values1<0])
                SEG_VAL=MODE_values[0]
                print '---------------------------------'
                print 'Neighbor value found:',-1*SEG_VAL
                print '---------------------------------'
            if len(values1[values1<0])<=0:
                print '-------------------------------------------------------------'
                print 'No Neighbors found. A 5 pixel source is automatically created'
                print '-------------------------------------------------------------'
        
        if SEG_VAL==0:
            HST_seg_base[0].data[int(round(Y_IMAGE[iii]-1)):int(round(Y_IMAGE[iii]+2)),int(round(X_IMAGE[iii]-1)):int(round(X_IMAGE[iii]+2))]=HST_ID[iii]
            minx=X_IMAGE[iii]-1 # NOTE: X_IMAGE and Y_IMAGE start from 1 (not from 0)
            miny=Y_IMAGE[iii]-1 # NOTE: X_IMAGE and Y_IMAGE start from 1 (not from 0)
            maxx=X_IMAGE[iii]+1 # NOTE: X_IMAGE and Y_IMAGE start from 1 (not from 0)
            maxy=Y_IMAGE[iii]+1 # NOTE: X_IMAGE and Y_IMAGE start from 1 (not from 0)
            # Unexpected behaviour - end ----------------------------------------------
        if SEG_VAL!=0:
            # Now associate the new value, given to this source in the catalog 
            DDD=[np.where(HST_seg[0].data==SEG_VAL)] # one dimensional index array
            EEE=np.where(HST_seg[0].data==SEG_VAL) # multi-dimensional index array
            EEE=np.array(EEE)
            YYY=EEE[0]
            XXX=EEE[1]
            #
            DDD=np.array(DDD)
            DDD=DDD[0]
            HST_seg_base[0].data[np.where(HST_seg[0].data==SEG_VAL)]=HST_ID[iii]
            minx=np.min(XXX)+1 # NOTE: +1 is required since python indexes start from 0 (not 1)
            miny=np.min(YYY)+1 # NOTE: +1 is required since python indexes start from 0 (not 1)
            maxx=np.max(XXX)+1 # NOTE: +1 is required since python indexes start from 0 (not 1)
            maxy=np.max(YYY)+1 # NOTE: +1 is required since python indexes start from 0 (not 1)
    
        X_MIN[iii]=minx
        Y_MIN[iii]=miny
        X_MAX[iii]=maxx
        Y_MAX[iii]=maxy
        iii=iii+1

    if INTERACT == 'I': 
        plt.imshow(np.log(HST_seg_base[0].data), cmap='gray')
        plt.title('Final Segmentation map (Intensity reflects brightness)')
        plt.colorbar()
        plt.show()

    # super-Segmentation map saved into output file
    HST_seg_base.writeto(paths.pwd+files.out_seg)
    # Save output catalog (tphot input)
    ALLDATA = np.zeros(len(HST_ID), dtype=[('HST_ID', int), ('X_IMAGE', float), ('Y_IMAGE', float), ('X_MIN', int), ('Y_MIN', int), ('X_MAX', int), ('Y_MAX', int), ('LOC_BCK', float), ('FLUX', float)  ])
    ALLDATA['HST_ID']=HST_ID
    ALLDATA['X_IMAGE']=X_IMAGE
    ALLDATA['Y_IMAGE']=Y_IMAGE
    ALLDATA['X_MIN']=X_MIN
    ALLDATA['Y_MIN']=Y_MIN
    ALLDATA['X_MAX']=X_MAX
    ALLDATA['Y_MAX']=Y_MAX
    ALLDATA['LOC_BCK']=LOC_BCK # null value.
    ALLDATA['FLUX']=FLUX
    ALLDATA=ALLDATA.T # TRANSPOSE!!!
    STRING=" %6i %11.4f %11.4f %6i %6i %6i %6i %15.6f %15.6f "
    np.savetxt(paths.pwd+files.HRI_catalog,ALLDATA, fmt=STRING)


def make_output_dir(DS,paths,files,PARN):
    # creates an output directory, copies the main outputs in it and then tar.gz
    #############################################################################
    # 1) make output directory
    os.system('mkdir '+paths.pwd+files.output_folder)

    # 2) copy main outputs in it

    # HST segmentation map
    os.system('cp '+paths.pwd+files.out_seg+' '+paths.pwd+files.output_folder+'/')
    # PSF image HST (reference filter)
    os.system('cp '+paths.pwd+files.PSF_HST+' '+paths.pwd+files.output_folder+'/')
    
    if DS.I1==True:
        # Tphot catalog
        os.system('cp '+paths.pwd+files.IRAC_cat_tphot_1+' '+paths.pwd+files.output_folder+'/')
        # SExtractor catalog (original image)
        os.system('cp '+paths.pwd+files.SEx_complete_cat_I1+' '+paths.pwd+files.output_folder+'/')
        # Irac resampled image (sci)
        os.system('cp '+paths.pwd+files.IRAC_resamp1+' '+paths.pwd+files.output_folder+'/')
        # Irac resampled image (unc)
        os.system('cp '+paths.pwd+files.IRAC_unc_resamp1+' '+paths.pwd+files.output_folder+'/')
        # Residuals map
        os.system('cp '+paths.pwd+'WISPS'+str(PARN)+'_ch1_mosaic.resamp_bck_sub_SingleFit_CLIP/resid_lores_collage_CH1_pass2.fits '+paths.pwd+files.output_folder+'/')
        # PSF image IRAC
        os.system('cp '+paths.pwd+files.PSF_IRAC1+' '+paths.pwd+files.output_folder+'/')
        # Kernel image IRAC1-HST
        os.system('cp '+paths.pwd+files.KERNEL_1+' '+paths.pwd+files.output_folder+'/')
        
        
    if DS.I2==True:
        # Tphot catalog
        os.system('cp '+paths.pwd+files.IRAC_cat_tphot_2+' '+paths.pwd+files.output_folder+'/')
        # SExtractor catalog (original image)
        os.system('cp '+paths.pwd+files.SEx_complete_cat_I2+' '+paths.pwd+files.output_folder+'/')
        # Irac resampled image (sci)
        os.system('cp '+paths.pwd+files.IRAC_resamp2+' '+paths.pwd+files.output_folder+'/')
        # Irac resampled image (unc)
        os.system('cp '+paths.pwd+files.IRAC_unc_resamp2+' '+paths.pwd+files.output_folder+'/')
        # Residuals map
        os.system('cp '+paths.pwd+'WISPS'+str(PARN)+'_ch2_mosaic.resamp_bck_sub_SingleFit_CLIP/resid_lores_collage_CH2_pass2.fits '+paths.pwd+files.output_folder+'/')
        # PSF image IRAC
        os.system('cp '+paths.pwd+files.PSF_IRAC2+' '+paths.pwd+files.output_folder+'/')
        # Kernel image IRAC2-HST
        os.system('cp '+paths.pwd+files.KERNEL_2+' '+paths.pwd+files.output_folder+'/')

    # 3) tar.gz this folder (and eliminate original)

    PATH_NOW=os.getcwd() # get current directory
    os.chdir(paths.pwd) # go to the working directory where the folder has been created
    os.system('tar -cvf '+files.output_folder+'.tar '+files.output_folder)
    os.system('gzip '+files.output_folder+'.tar')
    shutil.rmtree(files.output_folder, ignore_errors=True)
    os.chdir(PATH_NOW) #return to current directory
    
    


def clean_after2(DS,paths,files,PARN):
    # Remove files already saved elsewere.
    #############################################################################
    if DS.J==True:
        os.system('rm '+paths.pwd+'F110W_sci.fits')
        os.system('rm '+paths.pwd+'F110W_seg.fits')
    if DS.H1==True:
        os.system('rm '+paths.pwd+'F140W_sci.fits')
        os.system('rm '+paths.pwd+'F140W_seg.fits')
    if DS.H2==True:
        os.system('rm '+paths.pwd+'F160W_sci.fits')
        os.system('rm '+paths.pwd+'F160W_seg.fits')
    if (DS.J==True) and ((DS.H1==True) or (DS.H2==True)):
        os.system('rm '+paths.pwd+'JH_combined_sci.fits')
        os.system('rm '+paths.pwd+'JH_combined_seg.fits')

    if DS.I1==True:
        os.system('rm '+paths.pwd+files.IRAC_cov1)
    if DS.I2==True:
        os.system('rm '+paths.pwd+files.IRAC_cov2)
        
    

#######################################################
####################################################
################################################
# MAIN PROGRAM STARTS HERE
################################################
####################################################
#######################################################

# Run the program giving the Par number as an input.
# Example:
# > python IRAC_photometry.py 322


def main():

    ##############################################################
    # Initialization of the working directory
    ##############################################################
    wherearewe=os.getcwd() # At the end, return to this directory
    TS=0                # Test or real run? Test--> TS=1, Normal run TS=0 
    Narg=check_arg()-1  # Check that the correct number of arguments is given in input (Par#)
    PARN=sys.argv[1]    # Set the Par number to the first input argument
    if Narg==1:INTERACT='A' # Default input argument for run MODE 
    if Narg>=2:INTERACT=sys.argv[len(sys.argv)-1] # last argument for MODE ('I'nteractive or 'A'utomatic, 'T'est)
                            # 'I' interactive
                            # 'A' automatic, to see all the plots,
                            # 'T' test. Uses different paths and is interactive
    if INTERACT=='T':
        INTERACT='I' # set INTEARCT keyword to Interactive mode.
        TS=1         # activate test mode
    if INTERACT=='P':
        INTERACT='A' # set INTEARCT keyword to Interactive mode.
        TS=1         # activate test mode
        
    paths=Get_paths(TS) # get path definitions. Argument=1 for test, 0 for real run
    clean_old(paths,PARN,TS) # Clean old files made by previous (possibly failed?) runs
    create_tree(paths)       # creates directory three
    DS=Data_def(paths,PARN)         # Define Data-Set(check which kind of data are present)
    DS.show()                       # Shows kind of dataset used
    files=Get_Filenames(PARN,DS)    # Set files names
    copy_files(paths,files,DS)      # Copy original data into the working directory and
                                    # creates an additional backup copy (folder: original)

    ##############################################################
    # Read REFERENCE J+H catalog (catalog made by WISP pipeline)
    ##############################################################
    cat_HST=Readcol(paths.pwd+files.HST_catalog,'a,i,f,f,f,f,a,f,f,a,a,a,f',skipline=17)
    # ---- needed for TPHOT-input catalog: ----
    WISP_ID=cat_HST.col(0)
    HST_ID=cat_HST.col(1)
    X_IMAGE=cat_HST.col(2)
    Y_IMAGE=cat_HST.col(3)
    A_IMAGE=cat_HST.col(4)
    B_IMAGE=cat_HST.col(5)
    # ---- needed for recentering: ----
    X_WORLD=cat_HST.col(7)
    Y_WORLD=cat_HST.col(8)
    MAG=cat_HST.col(12)
    # compute flux [uJy] (HST J or H)
    FLUX=10**((23.9-MAG)/2.5)
    ##############################################################
    # ITERATE: SExtract IRAC img - MATCH with HST - RECENTER IRAC imgs
    ##############################################################
    total_RA_shift_I1=0.
    total_DEC_shift_I1=0.
    total_RA_shift_I2=0.
    total_DEC_shift_I2=0.

    ITER=0
    TOT_ITER=2 # Total number of recentering itarations
    # An additional iteration is required to extract the final catalog
    while ITER < TOT_ITER+1:
        print '---------------------'
        print 'ITERATION:', str(ITER)
        print '---------------------'
        # catalogs and regions code: 
        # tmp="temporary"
        # it#= iteration number
        if ITER<TOT_ITER:
            catname1=paths.p3+'tmp_cat_CH1_it'+str(ITER)+'.fits'
            catname2=paths.p3+'tmp_cat_CH2_it'+str(ITER)+'.fits'
            region_name1='tmp_reg_CH1_it'+str(ITER)+'.reg'
            region_name2='tmp_reg_CH2_it'+str(ITER)+'.reg'
        if ITER==TOT_ITER:
            catname1=paths.p3+files.IRAC_cat_SEx_LRES_1 # Final CH1 catalog
            catname2=paths.p3+files.IRAC_cat_SEx_LRES_2 # Final CH2 catalog
            region_name1="CH1_recentered.reg"
            region_name2="CH2_recentered.reg"
        
        RADIUS=5 # pixels radius for region files

        ##############################################################
        # SExtract IRAC images and write region files
        ##############################################################
        if DS.I1 == True:
            # Compute GAIN using coverage images then run sextractor
            SExtract('CH1',catname1,files,paths,0.6)
            # Read output catalog created by SExtractor
            cat_I1 = fits.getdata(catname1)
            #print cat_I1.columns # print columns in the catalog
            # write region files
            mk_regionfile(cat_I1.ALPHA_J2000,cat_I1.DELTA_J2000,RADIUS,file=paths.pwd+region_name1)
            if ITER==0:
                # opy IRAC SExtractor catalogs referring to the original images and NOT RECENTERED to HST YET!!,
                # on the working directory. This is one of the outputs of this program, since tphot can't compute
                # IRAC photometry on sources outside the HST area otherwise. The IRAC original coordinates are
                # more precise than HST and the final TPHOT catalog will contain IRAC original coordinates too.
                os.system('cp '+catname1+' '+paths.pwd+files.SEx_complete_cat_I1)

        if DS.I2 == True:
            # Compute GAIN using coverage images then run sextractor
            SExtract('CH2',catname2,files,paths,0.6)
            # Read output catalog created by SExtractor
            cat_I2 = fits.getdata(catname2)
            print cat_I2.columns # print columns in the catalog
            # Region files from SExtracted catalogs
            # write region files
            mk_regionfile(cat_I2.ALPHA_J2000,cat_I2.DELTA_J2000,RADIUS,file=paths.pwd+region_name2)
            if ITER==0:
                # opy IRAC SExtractor catalogs referring to the original images and NOT RECENTERED to HST YET!!,
                # on the working directory. This is one of the outputs of this program, since tphot can't compute
                # IRAC photometry on sources outside the HST area otherwise. The IRAC original coordinates are
                # more precise than HST and the final TPHOT catalog will contain IRAC original coordinates too.
                os.system('cp '+catname2+' '+paths.pwd+files.SEx_complete_cat_I2)

        ##############################################################
        # Match coordinates in the catalogs (IRAC Vs HST)
        ##############################################################
        if ITER==0: max_DIST=2.
        if ITER!=0: max_DIST=0.75
        
        
        # IRAC Channel 1
        if DS.I1 == True: 
            N1=len(cat_I1.ALPHA_J2000)
            Nnot1=N1*5/100 # 5% brightest sources that should not be considered (probably saturated)
            Nnot1B=N1*(100-10)/100 # 10% faintest sources that should not be considered (to faint)
            idx1=cat_I1.MAG_AUTO.argsort() # Sort magnitudes in the catalog
            # Match after removing the saturated and faint sources
            MATCH1=Match_cat(cat_I1.ALPHA_J2000[idx1[Nnot1:Nnot1B]],cat_I1.DELTA_J2000[idx1[Nnot1:Nnot1B]],X_WORLD,Y_WORLD,max_DIST)
            mk_regionfile(cat_I1.ALPHA_J2000[idx1[Nnot1:Nnot1B]],cat_I1.DELTA_J2000[idx1[Nnot1:Nnot1B]],RADIUS-1,color='green',file=paths.pwd+'srcs_used_recenter_CH1.reg')
            if INTERACT == 'I': 
                if ITER<TOT_ITER: ptit1='IRAC CH1 recentering - Iteration '+str(ITER)
                if ITER==TOT_ITER: ptit1='IRAC CH1 recentering - FINAL RESULTS'
                MATCH1.plot_RA_DEC(title=ptit1)        # plot Ra Vs Dec
                MATCH1.plot_delta_RA_DEC(title=ptit1)  # plot Delta Ra Vs Delta Dec
                MATCH1.plot_delta2_RA_DEC(title=ptit1) # plot Delta Ra and Delta Dec as a function of Ra and dec
        # IRAC Channel 2
        if DS.I2 == True:
            N2=len(cat_I2.ALPHA_J2000)
            Nnot2=N2*5/100 # 10% brightest sources that should not be considered (probably saturated)
            Nnot2B=N2*(100-10)/100 # 10% faintest sources that should not be considered (to faint)
            idx2=cat_I2.MAG_AUTO.argsort() # Sort magnitudes in the catalog
            # Match after removing the saturated and faint sources
            MATCH2=Match_cat(cat_I2.ALPHA_J2000[idx2[Nnot2:Nnot2B]],cat_I2.DELTA_J2000[idx2[Nnot2:Nnot2B]],X_WORLD,Y_WORLD,max_DIST)
            mk_regionfile(cat_I2.ALPHA_J2000[idx2[Nnot2:Nnot2B]],cat_I2.DELTA_J2000[idx2[Nnot2:Nnot2B]],RADIUS-1,color='green',file=paths.pwd+'srcs_used_recenter_CH2.reg')
            if INTERACT == 'I': 
                if ITER<TOT_ITER: ptit2='IRAC CH2 recentering - Iteration '+str(ITER)
                if ITER==TOT_ITER: ptit2='IRAC CH2 recentering - FINAL RESULTS'
                MATCH2.plot_RA_DEC(title=ptit2)        # plot Ra Vs Dec
                MATCH2.plot_delta_RA_DEC(title=ptit2)  # plot Delta Ra Vs Delta Dec
                MATCH2.plot_delta2_RA_DEC(title=ptit2) # plot Delta Ra and Delta Dec as a function of Ra and dec

        if ITER<TOT_ITER: # RECENTER only TWO TIMES
            ##############################################################
            # Recenter IRAC images, uncertainty images and coverage images 
            ##############################################################
            # Recenter IRAC images using a simple shift correction.
            #  To this purpose, the reference pixel is shifted for the 
            #  corresponding delta_ra, delta_dec computed before.
            if DS.I1 == True: 
                image_shift(paths.pwd+files.IRAC_image1,dist_RA=np.median(MATCH1.DIST_RA),dist_DEC=np.median(MATCH1.DIST_DEC))
                image_shift(paths.pwd+files.IRAC_unc1,dist_RA=np.median(MATCH1.DIST_RA),dist_DEC=np.median(MATCH1.DIST_DEC))
                image_shift(paths.pwd+files.IRAC_cov1,dist_RA=np.median(MATCH1.DIST_RA),dist_DEC=np.median(MATCH1.DIST_DEC))
                total_RA_shift_I1=total_RA_shift_I1+np.median(MATCH1.DIST_RA)
                total_DEC_shift_I1=total_DEC_shift_I1+np.median(MATCH1.DIST_DEC)
            if DS.I2 == True: 
                image_shift(paths.pwd+files.IRAC_image2,dist_RA=np.median(MATCH2.DIST_RA),dist_DEC=np.median(MATCH2.DIST_DEC))
                image_shift(paths.pwd+files.IRAC_unc2,dist_RA=np.median(MATCH2.DIST_RA),dist_DEC=np.median(MATCH2.DIST_DEC))
                image_shift(paths.pwd+files.IRAC_cov2,dist_RA=np.median(MATCH2.DIST_RA),dist_DEC=np.median(MATCH2.DIST_DEC))
                total_RA_shift_I2=total_RA_shift_I2+np.median(MATCH2.DIST_RA)
                total_DEC_shift_I2=total_DEC_shift_I2+np.median(MATCH2.DIST_DEC)

        ITER=ITER+1


    ##############################################################
    # Compute JH image reference values for SWarp, from headers
    ##############################################################
    # get center, size and orientation of the reference HST image
    # and write two header files (one for IRAC ch-1 and one for ch-2).
    # In particular, the output header files contain the orientation
    # parameters of the reference image. They will be used by SWarp
    # to set the pixel orientation of the output images.
    AAA=size_orient_from_ref(paths.pwd+files.HST_image,paths.pwd+files.headname1,paths.pwd+files.headname2)
    string_center=AAA[0]
    string_size=AAA[1]

    ##############################################################
    # Swarp
    ##############################################################
    # Run SWARP into the working directory. 
    # ONLY on IRAC images (HST catalogs, images and segmentation maps 
    # are those created during the run of the wisp pipeline.)
    # Note: paths have now to be considered with respect to this dir.
    os.chdir(paths.pwd)
    if DS.I1 == True: 
        # The orientation must be specified in a file "filename.head" where filename must have the same name as the output image
        filetest1=os.path.isfile(files.headname1)
        if filetest1 == False: errors(files,1)
        if filetest1 == True: resample('CH1',files,paths,string_center,string_size) # sci and unc resample + bck subtraction

    if DS.I2 == True: 
        # The orientation must be specified in a file "filename.head" where filename must have the same name as the output image
        filetest2=os.path.isfile(files.headname2)
        if filetest2 == False: errors(files,2)
        if filetest2 == True:resample('CH2',files,paths,string_center,string_size) # sci and unc resample + bck subtraction

    os.chdir(wherearewe) # Return back to pyton program directory

    ##############################################################
    # Clean & rename remove files to save space + rename SWARP out
    ##############################################################
    clean_after(DS,paths,files,PARN)
    
    ##############################################################
    # Subtract background from HST images
    ##############################################################
    # NOTE: we are not using the SWarped HST images (that could be 
    # background subtracted by SWarp). This choise is made because 
    # there must be consistence among the catalog created during
    # the run of the wisp pipeline and the correspondent segmentation
    # map and image. If we use the swarped image, this correspondence
    # is lost at the level of single pixels. Example: the x,y center 
    # of a source in the swarped image could shift of a pixel, losing
    # correspondence with the sementation map.
    #
    background_HST(DS,paths,files)

    ##############################################################
    # Creates HST segmentation map and catalog with TPHOT requirements  
    ##############################################################
    # Creates a segmentation map and a consistent tphot input catalog
    # Note that for this passage, the segmentation maps created during 
    # the WISP reduction process of the HST images are needed.  
    Segment_and_cat(paths,files,HST_ID,X_IMAGE,Y_IMAGE,FLUX,INTERACT)

    ##############################################################
    # Kernel from deconvolution:  PSF(LRI) = Kernel * PSF(HRI) 
    ##############################################################
    # We compute the convolution kernel from the averaged PSFs of HST and IRAC,
    # obtained directly from the images of the field, for point-like sources.
    # This approach is preferred to the use of synthetic PSFs for the
    # following reason: IRAC PSF is not rotationally symmetric
    # (it looks more like a triangle than a circle). The orientation
    # of this triangle is different from field to field and it is also
    # the result of a SWarping (re-drizzling) process 
    # on a smaller pixel-scale.
    #
    #
    # The kernel image shoud be obtained using a Deconvolution algorithm
    # see:
    # https://stackoverflow.com/questions/17473917/is-there-a-equivalent-of-scipy-signal-deconvolve-for-2d-arrays
    # Input needed: 

    RADIUS_HRES=15  # For IRAC H-res catalog test (C)
    max_DIST_TEST=1 # For IRAC H-res catalog test (C)

    if DS.I1==True:
        # A) SExtract IRAC SWarped Image (High resolution) to get sources
        #    positions in the high resolution version of the IRAC image
        catname1B=paths.p3+files.IRAC_cat_SEx_HRES_1 #+'tmp_cat_CH1_HRES.fits'
        SExtract('CH1',catname1B,files,paths,0.08)
        # B) Load recentered catalog
        cat_I1_HRES = fits.getdata(catname1B)
        #---------------------------------------------------------------
        # C) TEST the new IRAC High resolution catalog
        # C1) Write region files
        region_name1B="CH1_recentered_HRES.reg"
        mk_regionfile(cat_I1_HRES.ALPHA_J2000,cat_I1_HRES.DELTA_J2000,RADIUS_HRES,file=paths.pwd+region_name1B,color='green')
        # C2) Match with low_res and compare magnitudes and associated uncertainties.
        MATCH_1_TEST=Match_cat(cat_I1_HRES.ALPHA_J2000,cat_I1_HRES.DELTA_J2000,cat_I1.ALPHA_J2000,cat_I1.DELTA_J2000,max_DIST_TEST)
        if INTERACT == 'I': 
            MATCH_1_TEST.plot_RA_DEC(title='CH1 Match between High and Low res IRAC images')        # plot Ra Vs Dec
            MATCH_1_TEST.plot_delta_RA_DEC(title='CH1 Match between High and Low res IRAC images')  # plot Delta Ra Vs Delta Dec
            MATCH_1_TEST.plot_delta2_RA_DEC(title='CH1 Match between High and Low res IRAC images') # plot Delta Ra and Delta Dec as a function of Ra and dec
        # MAGNITUDES -
            plt.plot(np.array([0,100]),np.array([0,100]))
            plt.plot(cat_I1.MAG_AUTO[MATCH_1_TEST.IDX2],cat_I1_HRES.MAG_AUTO[MATCH_1_TEST.IDX1], color='red',marker='.',linestyle='none')
            axes = plt.gca()
            axes.set_xlim([20,26])
            axes.set_ylim([20,26])
            plt.ylabel('I1 mag (HRES)')
            plt.xlabel('I1 mag (LRES)')
            plt.title('SExtractor CH1 magnitude (High Vs Low resolution)')
            plt.show()
        # MAGNITUDE UNCERTAINTIES -  
            plt.plot(np.array([0,100]),np.array([0,100]))
            plt.plot(cat_I1.MAGERR_AUTO[MATCH_1_TEST.IDX2],cat_I1_HRES.MAGERR_AUTO[MATCH_1_TEST.IDX1], color='red',marker='.',linestyle='none')
            axes = plt.gca()
            axes.set_xlim([0.,1.])
            axes.set_ylim([0.,1.])
            plt.ylabel('I1 mag unc. (HRES)')
            plt.xlabel('I1 mag unc. (LRES)')
            plt.title('SExtractor CH1 magnitude uncertainty (High Vs Low resolution)')
            plt.show()
        #---------------------------------------------------------------

        # D) Match IRAC-Hres with HST (Hres)
        max_DIST1_B=0.5 # ONLY SECURE MATCHES USED!
        MATCH1B=Match_cat(cat_I1_HRES.ALPHA_J2000,cat_I1_HRES.DELTA_J2000,X_WORLD,Y_WORLD,max_DIST1_B)
        if INTERACT == 'I': 
            MATCH1B.plot_RA_DEC(title='Match IRAC-1 SWarped with HST')        # plot Ra Vs Dec
            MATCH1B.plot_delta_RA_DEC(title='Match IRAC-1 SWarped with HST')  # plot Delta Ra Vs Delta Dec
            MATCH1B.plot_delta2_RA_DEC(title='Match IRAC-1 SWarped with HST') # plot Delta Ra and Delta Dec as a function of Ra and dec

    if DS.I2==True:
        # A) SExtract IRAC SWarped Image (High resolution) to get sources
        #    positions in the high resolution version of the IRAC image
        catname2B=paths.p3+files.IRAC_cat_SEx_HRES_2 #+'tmp_cat_CH2_HRES.fits'
        SExtract('CH2',catname2B,files,paths,0.08)
        # B) Load recentered catalog
        cat_I2_HRES = fits.getdata(catname2B)
        #---------------------------------------------------------------
        # C) TEST the new IRAC High resolution catalog
        # C1) Write region files
        region_name2B="CH2_recentered_HRES.reg"
        mk_regionfile(cat_I2_HRES.ALPHA_J2000,cat_I2_HRES.DELTA_J2000,RADIUS_HRES,file=paths.pwd+region_name2B,color='green')
        # C2) Match with low_res and compare magnitudes and associated uncertainties.
        MATCH_2_TEST=Match_cat(cat_I2_HRES.ALPHA_J2000,cat_I2_HRES.DELTA_J2000,cat_I2.ALPHA_J2000,cat_I2.DELTA_J2000,max_DIST_TEST)
        if INTERACT == 'I': 
            MATCH_2_TEST.plot_RA_DEC(title='CH2 Match between High and Low res IRAC images')        # plot Ra Vs Dec
            MATCH_2_TEST.plot_delta_RA_DEC(title='CH2 Match between High and Low res IRAC images')  # plot Delta Ra Vs Delta Dec
            MATCH_2_TEST.plot_delta2_RA_DEC(title='CH2 Match between High and Low res IRAC images') # plot Delta Ra and Delta Dec as a function of Ra and dec
        # MAGNITUDES -
            plt.plot(np.array([0,100]),np.array([0,100]))
            plt.plot(cat_I2.MAG_AUTO[MATCH_2_TEST.IDX2],cat_I2_HRES.MAG_AUTO[MATCH_2_TEST.IDX1], color='red',marker='.',linestyle='none')
            axes = plt.gca()
            axes.set_xlim([20,26])
            axes.set_ylim([20,26])
            plt.ylabel('I2 mag (HRES)')
            plt.xlabel('I2 mag (LRES)')
            plt.title('SExtractor CH2 magnitude (High Vs Low resolution)')
            plt.show()
        # MAGNITUDE UNCERTAINTIES -  
            plt.plot(np.array([0,100]),np.array([0,100]))
            plt.plot(cat_I2.MAGERR_AUTO[MATCH_2_TEST.IDX2],cat_I2_HRES.MAGERR_AUTO[MATCH_2_TEST.IDX1], color='red',marker='.',linestyle='none')
            axes = plt.gca()
            axes.set_xlim([0.,1.])
            axes.set_ylim([0.,1.])
            plt.ylabel('I2 mag unc. (HRES)')
            plt.xlabel('I2 mag unc. (LRES)')
            plt.title('SExtractor CH2 magnitude uncertainty  (High Vs Low resolution)')
            plt.show()
        #---------------------------------------------------------------

        # D) Match IRAC-Hres with HST (Hres)
        max_DIST2_B=0.5 # ONLY SECURE MATCHES USED!
        MATCH2B=Match_cat(cat_I2_HRES.ALPHA_J2000,cat_I2_HRES.DELTA_J2000,X_WORLD,Y_WORLD,max_DIST2_B)
        if INTERACT == 'I': 
            MATCH2B.plot_RA_DEC(title='Match IRAC-2 SWarped with HST')        # plot Ra Vs Dec
            MATCH2B.plot_delta_RA_DEC(title='Match IRAC-2 SWarped with HST')  # plot Delta Ra Vs Delta Dec
            MATCH2B.plot_delta2_RA_DEC(title='Match IRAC-2 SWarped with HST') # plot Delta Ra and Delta Dec as a function of Ra and dec





    # E)  Select not saturated - not too faint and approximately point-like sources from the HST catalog
    median_HST_mag=np.median(MAG[MAG<np.max(MAG)]) # median HST Magnitude
    # minimum acceptable magnitude
    min_HST_MAG=np.max([median_HST_mag-3,np.min(MAG),MAG[int(round(len(MAG)*5/100.))] ]) # Don't change it. It is correct!
    # maximum acceptable magnitude
    max_HST_MAG=np.min([median_HST_mag+1,np.max(MAG[MAG<np.max(MAG)])]) # Don't change it! It is correct!

    SEL_OK_1=np.array([-1])
    SEL_OK_2=np.array([-1])

    if DS.I1==True:
        # SELECTION in the IRAC catalog using HST properties (magnitude, shape, dimension)
        # A/B<1.5
        # MAG< too faint
        # MAG> too bright
        # A<15 pixels
        SEL_OK_1=np.where((A_IMAGE[MATCH1B.IDX2]/B_IMAGE[MATCH1B.IDX2]<1.5) & (MAG[MATCH1B.IDX2]>min_HST_MAG) & (MAG[MATCH1B.IDX2]<max_HST_MAG) & (A_IMAGE[MATCH1B.IDX2]<15)) # NOTE: for sets, "&", not "and" must be used! 
        SEL_OK_1=SEL_OK_1[0]
        #   --> Load IRAC SWarped (and bck subtracted) Image (High resolution) 
        I1_resamp_img = fits.open(paths.pwd+files.IRAC_resamp1)

    if DS.I2==True:
        # SELECTION in the IRAC catalog using HST properties (magnitude, shape, dimension)
        # A/B<1.5
        # MAG< too faint
        # MAG> too bright
        # A<15 pixels
        SEL_OK_2=np.where((A_IMAGE[MATCH2B.IDX2]/B_IMAGE[MATCH2B.IDX2]<1.5) & (MAG[MATCH2B.IDX2]>min_HST_MAG) & (MAG[MATCH2B.IDX2]<max_HST_MAG) & (A_IMAGE[MATCH2B.IDX2]<15)) # NOTE: for sets, "&", not "and" must be used! 
        SEL_OK_2=SEL_OK_2[0]
        #   --> Load IRAC SWarped (and bck subtracted) Image (High resolution) 
        I2_resamp_img = fits.open(paths.pwd+files.IRAC_resamp2)



    #   --> Load HST bck subtracted Image
    HST_img_bck_sub = fits.open(paths.pwd+files.HST_bck_sub)
    
    # Compute number of cycles to perform to compute average stamps (PSFs)
    if (DS.I1==True) & (DS.I2==False): N_stamps=len(SEL_OK_1)
    if (DS.I1==False) & (DS.I2==True): N_stamps=len(SEL_OK_2)
    if (DS.I1==True) & (DS.I2==True):  N_stamps=len(SEL_OK_1)+len(SEL_OK_2)

    if DS.I1==True:
        # Create void stamp for IRAC PSFs
        PSF_STAMP_I1=np.empty([101, 101,len(SEL_OK_1)])
        PSF_STAMP_I1[:,:,:]=0
    if DS.I2==True:
        # Create void stamp for IRAC PSFs
        PSF_STAMP_I2=np.empty([101, 101,len(SEL_OK_2)])
        PSF_STAMP_I2[:,:,:]=0
        
    # Create void stamp for HST PSFs
    if DS.I1==True: PSF_STAMP_HST=np.empty([101, 101,len(SEL_OK_1)])
    if DS.I2==True & DS.I1==False:  PSF_STAMP_HST=np.empty([101, 101,len(SEL_OK_2)])
    PSF_STAMP_HST[:,:,:]=0


    # Compute average PSF, from images stamps, for IRAC and HST images
    WW=0
    while WW < N_stamps:
        # IRAC stamps 
        if (DS.I1==True) & (WW < len(SEL_OK_1)):
            STAMP_I1=I1_resamp_img[0].data[int(round(cat_I1_HRES.Y_IMAGE[MATCH1B.IDX1[SEL_OK_1[WW]]]-51)):int(round(cat_I1_HRES.Y_IMAGE[MATCH1B.IDX1[SEL_OK_1[WW]]]+50)),int(round(cat_I1_HRES.X_IMAGE[MATCH1B.IDX1[SEL_OK_1[WW]]]-51)):int(round(cat_I1_HRES.X_IMAGE[MATCH1B.IDX1[SEL_OK_1[WW]]]+50))]
            # plt.imshow(STAMP_I1, cmap='gray')
            # plt.colorbar()
            # plt.show()
            PSF_STAMP_I1[:,:,WW]=STAMP_I1/np.sum(abs(STAMP_I1))#*(10**cat_I1_HRES.MAG_AUTO[MATCH1B.IDX1[SEL_OK_1[WW]]])
       
        if (DS.I2==True) & (WW < len(SEL_OK_2)):
            STAMP_I2=I2_resamp_img[0].data[int(round(cat_I2_HRES.Y_IMAGE[MATCH2B.IDX1[SEL_OK_2[WW]]]-51)):int(round(cat_I2_HRES.Y_IMAGE[MATCH2B.IDX1[SEL_OK_2[WW]]]+50)),int(round(cat_I2_HRES.X_IMAGE[MATCH2B.IDX1[SEL_OK_2[WW]]]-51)):int(round(cat_I2_HRES.X_IMAGE[MATCH2B.IDX1[SEL_OK_2[WW]]]+50))]
            # plt.imshow(STAMP_I1, cmap='gray')
            # plt.colorbar()
            # plt.show()
            PSF_STAMP_I2[:,:,WW]=STAMP_I2/np.sum(abs(STAMP_I2))#*(10**cat_I1_HRES.MAG_AUTO[MATCH1B.IDX1[SEL_OK_1[WW]]])
        
        # HST stamps (same sources as IRAC!)
        if (DS.I1==True) & (WW < len(SEL_OK_1)):
            STAMP_HST=HST_img_bck_sub[0].data[int(round(Y_IMAGE[MATCH1B.IDX2[SEL_OK_1[WW]]]-51)):int(round(Y_IMAGE[MATCH1B.IDX2[SEL_OK_1[WW]]]+50)),int(round(X_IMAGE[MATCH1B.IDX2[SEL_OK_1[WW]]]-51)):int(round(X_IMAGE[MATCH1B.IDX2[SEL_OK_1[WW]]]+50))]
            PSF_STAMP_HST[:,:,WW]=STAMP_HST/np.sum(abs(STAMP_HST))
        # Only if IRAC-1 is not present, match with IRAC 2 is taken as a reference
        if (DS.I2==True) & (DS.I1==False): 
            STAMP_HST=HST_img_bck_sub[0].data[int(round(Y_IMAGE[MATCH2B.IDX2[SEL_OK_2[WW]]]-50)):int(round(Y_IMAGE[MATCH2B.IDX2[SEL_OK_2[WW]]]+51)),int(round(X_IMAGE[MATCH2B.IDX2[SEL_OK_2[WW]]]-50)):int(round(X_IMAGE[MATCH2B.IDX2[SEL_OK_2[WW]]]+51))]
            # plt.imshow(STAMP_HST, cmap='gray')
            # plt.colorbar()
            # plt.show()
            PSF_STAMP_HST[:,:,WW]=STAMP_HST/np.sum(abs(STAMP_HST))

        WW=WW+1


    # Average STAMP with average point-like source
    if DS.I1==True: AV_PSF_STAMP_I1=np.median(PSF_STAMP_I1[:,:,:],axis=2)
    if DS.I2==True: AV_PSF_STAMP_I2=np.median(PSF_STAMP_I2[:,:,:],axis=2)
    AV_PSF_STAMP_HST=np.median(PSF_STAMP_HST[:,:,:],axis=2)
    # Normalize STAMP with average point-like source
    if DS.I1==True: AV_PSF_STAMP_I1=AV_PSF_STAMP_I1/np.sum(AV_PSF_STAMP_I1)
    if DS.I2==True: AV_PSF_STAMP_I2=AV_PSF_STAMP_I2/np.sum(AV_PSF_STAMP_I2)
    AV_PSF_STAMP_HST=AV_PSF_STAMP_HST/np.sum(AV_PSF_STAMP_HST)
    # Write to a fits file. 
    # NOTE: fro a np.ndarray, writeto function is used differently than for an HDU!!
    if DS.I1==True: fits.writeto(paths.pwd+files.PSF_IRAC1,AV_PSF_STAMP_I1,clobber='True')
    if DS.I2==True: fits.writeto(paths.pwd+files.PSF_IRAC2,AV_PSF_STAMP_I2,clobber='True')
    fits.writeto(paths.pwd+files.PSF_HST,AV_PSF_STAMP_HST,clobber='True')
    # To use an hdu.writeto() function, first, convert np.ndarray to HDU
    # > I1_PSF_hdu = fits.PrimaryHDU(AV_PSF_STAMP_I1)
    # Then, Save PSF to a FITS file
    # > I1_PSF_hdu.writeto(paths.pwd+files.PSF_IRAC1,clobber='True')

    if INTERACT == 'I':
        plt.imshow(AV_PSF_STAMP_HST)#, cmap='gray')
        plt.plot([0,100],[50,50],color='red')
        plt.plot([50,50],[0,100],color='red')
        plt.title('Average HST PSF')
        plt.colorbar()
        plt.show()
        if DS.I1 == True:
            plt.imshow(AV_PSF_STAMP_I1)#), cmap='gray')
            plt.plot([0,100],[50,50],color='red')
            plt.plot([50,50],[0,100],color='red')
            plt.title('Average IRAC-1 PSF')
            plt.colorbar()
            plt.show()
        if DS.I2 == True:
            plt.imshow(AV_PSF_STAMP_I2)#), cmap='gray')
            plt.plot([0,100],[50,50],color='red')
            plt.plot([50,50],[0,100],color='red')
            plt.title('Average IRAC-2 PSF')
            plt.colorbar()
            plt.show()
    
    # http://scikit-image.org/docs/dev/auto_examples/filters/plot_deconvolution.html
    # Linux and OSX
    # scikit-image can be installed from the Python packaging index using
    # pip install -U scikit-image
    # Richardson-Lucy approach XXXXXXXXXXXXX 
    # The number of iterations must be manually calibrated
    if DS.I1==True: 
        KERNEL_I1  = restoration.richardson_lucy(AV_PSF_STAMP_I1, AV_PSF_STAMP_HST, iterations=5)
        KERNEL_I1=KERNEL_I1/np.sum(KERNEL_I1)
        fits.writeto(paths.pwd+files.KERNEL_1,KERNEL_I1,clobber='True')
        if INTERACT == 'I':
            plt.imshow(KERNEL_I1)
            plt.title('Kernel (PSF_IRAC1=Kernel*PSF_IRAC_1)')
            plt.colorbar()
            plt.show()
    if DS.I2==True: 
        KERNEL_I2  = restoration.richardson_lucy(AV_PSF_STAMP_I1, AV_PSF_STAMP_HST, iterations=5)
        KERNEL_I2=KERNEL_I2/np.sum(KERNEL_I2)
        fits.writeto(paths.pwd+files.KERNEL_2,KERNEL_I2,clobber='True')
        if INTERACT == 'I':
            plt.imshow(KERNEL_I2)
            plt.title('Kernel (PSF_IRAC2=Kernel*PSF_IRAC_2)')
            plt.colorbar()
            plt.show()

    
    ##############################################################
    # TPHOT first and second passages
    ##############################################################
    # CH1
    os.chdir(paths.pwd) # go to working directory
    
    tphot_out_cat_CH1="tphot_CH1"
    tphot_out_cat_CH2="tphot_CH2"
    pass1_str="_pass1.cat"
    pass2_str="_pass2.cat"

    if DS.I1==True:
        # Pass 1
        os.system('python '+paths.p9+'tphot par1.param -hiresfile '+files.HST_bck_sub+' -hirescat '+files.HRI_catalog+' -hiresseg '+files.out_seg+' -loresfile '+files.IRAC_resamp1+' -loreserr '+files.IRAC_unc_resamp1+' -kernelfile '+files.KERNEL_1+' -cutoutdir cutouts_CH1 -cutoutcat cutouts_CH1/_cutouts_CH1.cat -modelsdir models_CH1 -modelscat models_CH1/models_CH1.cat -templatedir templates_CH1_pass1 -templatecat templates_CH1_pass1/_templates_CH1_pass1.cat -fitpars tpipe_tphot_CH1_pass1.param -tphotcat '+tphot_out_cat_CH1+pass1_str+' -tphotcell lores_tphot_CH1_pass1.cell -tphotcovar lores_tphot_CH1_pass1.covar -modelfile lores_collage_CH1_pass1.fits -ddiagfile ddiags_CH1_passes_1to2.txt -dlogfile dlog_CH1.txt > out_CH1_pass1.log')
        # Pass 2
        os.system('python '+paths.p9+'tphot par2.param -hiresfile '+files.HST_bck_sub+' -hirescat '+files.HRI_catalog+' -hiresseg '+files.out_seg+' -loresfile '+files.IRAC_resamp1+' -loreserr '+files.IRAC_unc_resamp1+' -kernelfile '+files.KERNEL_1+' -cutoutdir cutouts_CH1 -cutoutcat cutouts_CH1/_cutouts_CH1.cat -modelsdir models_CH1 -modelscat models_CH1/models_CH1.cat -templatedir templates_CH1_pass2 -templatecat templates_CH1_pass1/_templates_CH1_pass1.cat -fitpars tpipe_tphot_CH1_pass2.param -tphotcat '+tphot_out_cat_CH1+pass2_str+' -tphotcell lores_tphot_CH1_pass2.cell -tphotcovar lores_tphot_CH1_pass2.covar -modelfile lores_collage_CH1_pass2.fits -ddiagfile ddiags_CH1_passes_1to2.txt -dlogfile dlog_CH1.txt > out_CH1_pass2.log')
        # copy output catalog in the parent directory (the working directory), with a different name
        os.system('cp '+'WISPS'+DS.PARN+'_ch1_mosaic'+files.resamp_bck_sub_suffix+'_SingleFit_CLIP/'+tphot_out_cat_CH1+pass2_str+' '+files.tphot_outcat_CH1)

    if DS.I2==True:
        # Pass 1
        os.system('python '+paths.p9+'tphot par1.param -hiresfile '+files.HST_bck_sub+' -hirescat '+files.HRI_catalog+' -hiresseg '+files.out_seg+' -loresfile '+files.IRAC_resamp2+' -loreserr '+files.IRAC_unc_resamp2+' -kernelfile '+files.KERNEL_2+' -cutoutdir cutouts_CH2 -cutoutcat cutouts_CH2/_cutouts_CH2.cat -modelsdir models_CH2 -modelscat models_CH2/models_CH2.cat -templatedir templates_CH2_pass1 -templatecat templates_CH2_pass1/_templates_CH2_pass1.cat -fitpars tpipe_tphot_CH2_pass1.param -tphotcat '+tphot_out_cat_CH2+pass1_str+' -tphotcell lores_tphot_CH2_pass1.cell -tphotcovar lores_tphot_CH2_pass1.covar -modelfile lores_collage_CH2_pass1.fits -ddiagfile ddiags_CH2_pass_1to2.txt -dlogfile dlog_CH2.txt > out_CH2_pass1.log')
        # Pass 2
        os.system('python '+paths.p9+'tphot par2.param -hiresfile '+files.HST_bck_sub+' -hirescat '+files.HRI_catalog+' -hiresseg '+files.out_seg+' -loresfile '+files.IRAC_resamp2+' -loreserr '+files.IRAC_unc_resamp2+' -kernelfile '+files.KERNEL_2+' -cutoutdir cutouts_CH2 -cutoutcat cutouts_CH2/_cutouts_CH2.cat -modelsdir models_CH2 -modelscat models_CH2/models_CH2.cat -templatedir templates_CH2_pass2 -templatecat templates_CH2_pass1/_templates_CH2_pass1.cat -fitpars tpipe_tphot_CH2_pass2.param -tphotcat '+tphot_out_cat_CH2+pass2_str+' -tphotcell lores_tphot_CH2_pass2.cell -tphotcovar lores_tphot_CH2_pass2.covar -modelfile lores_collage_CH2_pass2.fits -ddiagfile ddiags_CH2_pass_1to2.txt -dlogfile dlog_CH2.txt > out_CH2_pass2.log')
        # copy output catalog in the parent directory (the working directory), with a different name
        os.system('cp '+'WISPS'+DS.PARN+'_ch2_mosaic'+files.resamp_bck_sub_suffix+'_SingleFit_CLIP/'+tphot_out_cat_CH2+pass2_str+' '+files.tphot_outcat_CH2)

    os.chdir(wherearewe) # Return back to pyton program directory

    ##############################################################
    # Modify tphot output catalog to convert output fluxes (to uJy and AB mag)
    ##############################################################
    
    # A) Open TPHOT output catalog (HSD-id order is not preserved and must be re-ordered)
    if DS.I1==True:
        cat_tphot_I1=Readcol(paths.pwd+files.tphot_outcat_CH1,'i,f,f,i,f,f,f,f,f,f,f,i',skipline=12)
        TPHOT_I1_ID_o=cat_tphot_I1.col(0) # _o stands for 'original'
        TPHOT_I1_X_o=cat_tphot_I1.col(1)
        TPHOT_I1_Y_o=cat_tphot_I1.col(2)
        TPHOT_I1_FLUX_o=cat_tphot_I1.col(7)#*cat_tphot_I1.col(10)
        TPHOT_I1_FLUXERR_o=cat_tphot_I1.col(8)
        TPHOT_I1_FLAG_o=cat_tphot_I1.col(11)
    if DS.I2==True:
        cat_tphot_I2=Readcol(paths.pwd+files.tphot_outcat_CH2,'i,f,f,i,f,f,f,f,f,f,f,i',skipline=12)
        TPHOT_I2_ID_o=cat_tphot_I2.col(0) # _o stands for 'original'
        TPHOT_I2_X_o=cat_tphot_I2.col(1)
        TPHOT_I2_Y_o=cat_tphot_I2.col(2)
        TPHOT_I2_FLUX_o=cat_tphot_I2.col(7)#*cat_tphot_I1.col(10)
        TPHOT_I2_FLUXERR_o=cat_tphot_I2.col(8)
        TPHOT_I2_FLAG_o=cat_tphot_I2.col(11)

    #........................................
    # IRAC-1 ARRAYS (for output catalogs)
    #........................................
    # New arrays with same dimension as HST input catalog (WISP pipeline main output)
    TPHOT_I1_ID=np.zeros(len(HST_ID), dtype=np.int) # will contain the re-orderer IDs
    TPHOT_I1_X=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I1_Y=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I1_FLUX=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I1_FLUXERR=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I1_FLAG=np.zeros(len(HST_ID), dtype=float)
    # ------------------------------------------------
    # Magnitudes (not in the tphot output!)
    TPHOT_I1_MAG=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I1_MAGERR=np.zeros(len(HST_ID), dtype=float)
    # ------------------------------------------------

    #........................................
    # IRAC-2 ARRAYS (for output catalogs)
    #........................................
    # New arrays with same dimension as HST input catalog (WISP pipeline main output)
    TPHOT_I2_ID=np.zeros(len(HST_ID), dtype=np.int) # will contain the re-orderer IDs
    TPHOT_I2_X=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I2_Y=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I2_FLUX=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I2_FLUXERR=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I2_FLAG=np.zeros(len(HST_ID), dtype=float)
    # ------------------------------------------------
    # Magnitudes (not in the tphot output!)
    TPHOT_I2_MAG=np.zeros(len(HST_ID), dtype=float)
    TPHOT_I2_MAGERR=np.zeros(len(HST_ID), dtype=float)
    # ------------------------------------------------

    #........................................
    # IRAC- 1 & 2 ARRAYS (for output catalogs)
    #........................................
    # COORDINATES - valid for both I1 and I2
    WISP_ID_IRAC=np.zeros(len(HST_ID), dtype=float)
    HST_RA=np.zeros(len(HST_ID), dtype=float)
    HST_DEC=np.zeros(len(HST_ID), dtype=float)
    IRAC_RA_I1=np.zeros(len(HST_ID), dtype=float)
    IRAC_DEC_I1=np.zeros(len(HST_ID), dtype=float)
    IRAC_RA_I2=np.zeros(len(HST_ID), dtype=float)
    IRAC_DEC_I2=np.zeros(len(HST_ID), dtype=float)

    WISP_ID_IRAC=WISP_ID
    HST_RA=X_WORLD
    HST_DEC=Y_WORLD
    IRAC_RA_I1=HST_RA-total_RA_shift_I1/3600.   # Coordinates in the original IRAC image (not the recentered one. Use HST for that)
    IRAC_DEC_I1=HST_DEC-total_DEC_shift_I1/3600.# Coordinates in the original IRAC image (not the recentered one. Use HST for that)
    IRAC_RA_I2=HST_RA-total_RA_shift_I2/3600.   # Coordinates in the original IRAC image (not the recentered one. Use HST for that)
    IRAC_DEC_I2=HST_DEC-total_DEC_shift_I2/3600.# Coordinates in the original IRAC image (not the recentered one. Use HST for that)

    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    # Flux Conversion factor: 
    # A) 8.461595 to pass from image units MJy/str to uJy/pixel (see IRAC handbook)
    # B) NOT NEEDED --> ((0.6/0.08)**2) to convert to the new pixel area. 
    #                    THIS IS NOT NECESSARY SINCE IT IS ALREADY KEPT INTO 
    #                    ACCOUNT DURING THE SWARPING PHASE 
    conv_fact=8.461595# /((0.6/0.08)**2.)
    #@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    
    if DS.I1==True:
        Mid1=Match(TPHOT_I1_ID_o,HST_ID) # ID MATCH - CH1

        TPHOT_I1_ID[:]=HST_ID[:]
        TPHOT_I1_X[Mid1.IDX_2_1]=TPHOT_I1_X_o[Mid1.IDX_1_2]
        TPHOT_I1_Y[Mid1.IDX_2_1]=TPHOT_I1_Y_o[Mid1.IDX_1_2]
        # ..........................................................
        # Fluxes and uncertainties
        TPHOT_I1_FLUX[Mid1.IDX_2_1]=TPHOT_I1_FLUX_o[Mid1.IDX_1_2] *conv_fact # Flux converted to uJy !!
        # PUT NEGATIVE FLUXES TO 0
        TPHOT_I1_FLUX[TPHOT_I1_FLUX<=0]=0.
        # Flux Uncertainties
        TPHOT_I1_FLUXERR[Mid1.IDX_2_1]=TPHOT_I1_FLUXERR_o[Mid1.IDX_1_2]*conv_fact # Flux-err converted to uJy !!
        # PUT uncertainties for NEGATIVE FLUXES TO -99
        TPHOT_I1_FLUXERR[TPHOT_I1_FLUX<=0]=-99.
        # ..........................................................
        # Magnitudes and uncertainties
        TPHOT_I1_MAG[Mid1.IDX_2_1]=-((2.5*np.log10(TPHOT_I1_FLUX[Mid1.IDX_2_1]))-23.9)
        # PUT mag. for NEGATIVE FLUXES TO 99
        TPHOT_I1_MAG[TPHOT_I1_FLUX<=0]=99.
        # Magnitude Uncertainties
        TPHOT_I1_MAGERR[Mid1.IDX_2_1]=TPHOT_I1_MAG[Mid1.IDX_2_1]+((2.5*np.log10(TPHOT_I1_FLUX[Mid1.IDX_2_1]+TPHOT_I1_FLUXERR[Mid1.IDX_2_1]))-23.9)
        # PUT uncertainties for NEGATIVE FLUXES TO -99
        TPHOT_I1_MAGERR[TPHOT_I1_FLUX<=0]=-99.
        # ..........................................................
        # Flags:
        # 1 for negative fluxes (here put to F=0 and AB=-99)
        # 2 for blended sources
        # 4 sources in the border
        TPHOT_I1_FLAG[Mid1.IDX_2_1]=TPHOT_I1_FLAG_o[Mid1.IDX_1_2]

        # #####################
        # Compute image depth
        # #####################
        # We compute the image depth using the SWarped image [0.08"/pixel]
        # The depth is computed in a 21 pixels radius aperture, corresponding
        # to 21*0.08"/pix=1.68" radius (Given the 
        # ~2" FWHM of both I1 and I2 PSFs, 21 pixels correspond to ~2 PSF's sigma
        # for both IRAC 1 and 2) 
        go_interactive='A'
        if INTERACT=='I':go_interactive='I'
        DEPTH_I1_A=image_depth(paths.pwd+files.IRAC_resamp1,radius=21,INTERACT=go_interactive,title='I1 depth computation')
        # convet to uJy Flux using the conversion factor:
        DEPTH_I1_F=3*DEPTH_I1_A*conv_fact                  # 3-sigma depth flux
        # convert to AB magnitude
        DEPTH_I1_M=-((2.5*np.log10(DEPTH_I1_F))-23.9)      # 3-sigma depth magnitude
        # Fill magnitude array with the depth when sources are undetected
        TPHOT_I1_MAG[TPHOT_I1_MAG==99]=DEPTH_I1_M
        
        # #####################
        # Save output catalog
        # #####################

        # HEADER:        
        Header_cat_I1=" 0) RA_DEC_NAME     # Name in the WISP HST reference catalog \n 1) HST_ID          # ID in the WISP HST reference catalog \n 2) HST_RA          # RA in the WISP HST reference image and catalog \n 3) HST_DEC         # DEC in the WISP HST reference image and catalog \n 4) IRAC_RA_I1      # RA in the original IRAC-1 image (not recentered to HST image and cat.) \n 5) IRAC_DEC_I1     # DEC in the original IRAC-1 image (not recentered to HST image and cat.) \n 6) X_IMAGE_HRES    # x position in the IRAC-1 SWarped image \n 7) Y_IMAGE_HRES    # y position in the IRAC-1 SWarped image \n 8) AB_MAG_I1       # IRAC-1 AB magnitude from tphot output \n 9) MAGERR_I1       # IRAC-1 magnitude uncertainty from tphot output (-99 for undetected sources) \n 10) TPHOT_FLAG_I1  # tphot flag: +1 for negative fluxes, +2 for blended srcs, 4 for sources in the border \n"

        OUT_CAT_CH1=np.zeros(len(HST_ID), dtype=[('WISP_ID', "S"+str(len(HST_ID))), ('HST_ID', int),('HST_RA', float),('HST_DEC', float), ('IRAC_RA', float), ('IRAC_DEC', float),  ('X_IMAGE_HRES', float), ('Y_IMAGE_HRES', float), ('AB_MAG', float), ('MAGERR', float), ('TPHOT_FLAG', int) ])

                             
        OUT_CAT_CH1['WISP_ID']=WISP_ID_IRAC
        OUT_CAT_CH1['HST_ID']=TPHOT_I1_ID
        OUT_CAT_CH1['HST_RA']=HST_RA
        OUT_CAT_CH1['HST_DEC']=HST_DEC

        OUT_CAT_CH1['IRAC_RA']=IRAC_RA_I1
        OUT_CAT_CH1['IRAC_DEC']=IRAC_DEC_I1

        OUT_CAT_CH1['X_IMAGE_HRES']=TPHOT_I1_X
        OUT_CAT_CH1['Y_IMAGE_HRES']=TPHOT_I1_Y
        
        OUT_CAT_CH1['AB_MAG']=TPHOT_I1_MAG
        OUT_CAT_CH1['MAGERR']=TPHOT_I1_MAGERR

        OUT_CAT_CH1['TPHOT_FLAG']=TPHOT_I1_FLAG

        OUT_CAT_CH1=OUT_CAT_CH1.T # TRANSPOSE!!!
        STRING_I1=" %24s %6i %12.6f %12.6f %12.6f %12.6f %11.4f %11.4f %11.4f %11.4f %6i"

        np.savetxt(paths.pwd+files.IRAC_cat_tphot_1,OUT_CAT_CH1,header=Header_cat_I1, fmt=STRING_I1)

 

    if DS.I2==True:
        Mid2=Match(TPHOT_I2_ID_o,HST_ID) # ID MATCH - CH2

        TPHOT_I2_ID[:]=HST_ID[:]
        TPHOT_I2_X[Mid2.IDX_2_1]=TPHOT_I2_X_o[Mid2.IDX_1_2]
        TPHOT_I2_Y[Mid2.IDX_2_1]=TPHOT_I2_Y_o[Mid2.IDX_1_2]
        # ..........................................................
        # Fluxes and uncertainties
        TPHOT_I2_FLUX[Mid2.IDX_2_1]=TPHOT_I2_FLUX_o[Mid2.IDX_1_2] *conv_fact # Flux converted to uJy !!
        # PUT NEGATIVE FLUXES TO 0
        TPHOT_I2_FLUX[TPHOT_I2_FLUX<=0]=0.
        # Flux Uncertainties
        TPHOT_I2_FLUXERR[Mid2.IDX_2_1]=TPHOT_I2_FLUXERR_o[Mid2.IDX_1_2]*conv_fact # Flux-err converted to uJy !!
        # PUT uncertainties for NEGATIVE FLUXES TO -99
        TPHOT_I2_FLUXERR[TPHOT_I2_FLUX<=0]=-99.
        # ..........................................................
        # Magnitudes and uncertainties
        TPHOT_I2_MAG[Mid2.IDX_2_1]=-((2.5*np.log10(TPHOT_I2_FLUX[Mid2.IDX_2_1]))-23.9)
        # PUT mag. for NEGATIVE FLUXES TO 99
        TPHOT_I2_MAG[TPHOT_I2_FLUX<=0]=99.
        # Magnitude Uncertainties
        TPHOT_I2_MAGERR[Mid2.IDX_2_1]=TPHOT_I2_MAG[Mid2.IDX_2_1]+((2.5*np.log10(TPHOT_I2_FLUX[Mid2.IDX_2_1]+TPHOT_I2_FLUXERR[Mid2.IDX_2_1]))-23.9)
        # PUT uncertainties for NEGATIVE FLUXES TO -99
        TPHOT_I2_MAGERR[TPHOT_I2_FLUX<=0]=-99.
        # ..........................................................
        # Flags:
        # 1 for negative fluxes (here put to F=0 and AB=-99)
        # 2 for blended sources
        # 4 sources in the border
        TPHOT_I2_FLAG[Mid2.IDX_2_1]=TPHOT_I2_FLAG_o[Mid2.IDX_1_2]

        # #####################
        # Compute image depth
        # #####################
        # We compute the image depth using the SWarped image [0.08"/pixel]
        # The depth is computed in a 21 pixels radius aperture, corresponding
        # to 21*0.08"/pix=1.68" radius (Given the 
        # ~2" FWHM of both I1 and I2 PSFs, 21 pixels correspond to ~2 PSF's sigma
        # for both IRAC 1 and 2) 
        go_interactive='A'
        if INTERACT=='I':go_interactive='I' 
        DEPTH_I2_A=image_depth(paths.pwd+files.IRAC_resamp2,radius=21,INTERACT=go_interactive,title='I2 depth computation')
        # convet to uJy Flux using the conversion factor:
        DEPTH_I2_F=3*DEPTH_I2_A*conv_fact                  # 3-sigma depth flux
        # convert to AB magnitude
        DEPTH_I2_M=-((2.5*np.log10(DEPTH_I2_F))-23.9)      # 3-sigma depth magnitude
        # Fill magnitude array with the depth when sources are undetected
        TPHOT_I2_MAG[TPHOT_I2_MAG==99]=DEPTH_I2_M
        
        # #####################
        # Save output catalog
        # #####################

        # HEADER:
        
        Header_cat_I2=" 0) RA_DEC_NAME    # Name in the WISP HST reference catalog \n 1) HST_ID         # ID in the WISP HST reference catalog \n 2) HST_RA         # RA in the WISP HST reference image and catalog \n 3) HST_DEC        # DEC in the WISP HST reference image and catalog \n 4) IRAC_RA        # RA in the original IRAC-2 image (not recentered to HST image and cat.) \n 5) IRAC_DEC       # DEC in the original IRAC-2 image (not recentered to HST image and cat.) \n 6) X_IMAGE_HRES   # x position in the IRAC-2 SWarped image \n 7) Y_IMAGE_HRES   # y position in the IRAC-2 SWarped image \n 8) AB_MAG         # IRAC-2 AB magnitude from tphot output \n 9) MAGERR         # IRAC-2 magnitude uncertainty from tphot output (-99 for undetected sources) \n 10) TPHOT_FLAG    # tphot flag: +1 for negative fluxes, +2 for blended srcs, 4 for sources in the border \n"
        
        OUT_CAT_CH2=np.zeros(len(HST_ID), dtype=[('WISP_ID', "S"+str(len(HST_ID))), ('HST_ID', int),('HST_RA', float),('HST_DEC', float), ('IRAC_RA', float), ('IRAC_DEC', float),  ('X_IMAGE_HRES', float), ('Y_IMAGE_HRES', float), ('AB_MAG', float), ('MAGERR', float), ('TPHOT_FLAG', int) ])

        OUT_CAT_CH2['WISP_ID']=WISP_ID_IRAC
        OUT_CAT_CH2['HST_ID']=TPHOT_I2_ID
        OUT_CAT_CH2['HST_RA']=HST_RA
        OUT_CAT_CH2['HST_DEC']=HST_DEC

        OUT_CAT_CH2['IRAC_RA']=IRAC_RA_I2
        OUT_CAT_CH2['IRAC_DEC']=IRAC_DEC_I2

        OUT_CAT_CH2['X_IMAGE_HRES']=TPHOT_I2_X
        OUT_CAT_CH2['Y_IMAGE_HRES']=TPHOT_I2_Y
        
        OUT_CAT_CH2['AB_MAG']=TPHOT_I2_MAG
        OUT_CAT_CH2['MAGERR']=TPHOT_I2_MAGERR
        OUT_CAT_CH2['TPHOT_FLAG']=TPHOT_I2_FLAG

        OUT_CAT_CH2=OUT_CAT_CH2.T # TRANSPOSE!!!
        STRING_I2=" %24s %6i %12.6f %12.6f %12.6f %12.6f %11.4f %11.4f %11.4f %11.4f %6i"


        np.savetxt(paths.pwd+files.IRAC_cat_tphot_2,OUT_CAT_CH2, header=Header_cat_I2, fmt=STRING_I2)

    

    ##############################################################
    # Compare results of closest-counterparts association with TPHOT output
    ##############################################################
    # Here we compare the results (fluxes and uncertainties) of the 
    # IRAC sources found as counterparts for the HST priors,
    # using the ii) closest countepart association method and 
    # ii) the TPHOT method.
    # Note for the fluxes comparison: tphot-output Vs SExtractor
    if DS.I1==True:
        # Match HST (=Tphot output) Vs IRAC low resolution 
        max_DIST1_C=.5
        MATCH1C=Match_cat(cat_I1.ALPHA_J2000,cat_I1.DELTA_J2000,X_WORLD,Y_WORLD,max_DIST1_C)
        # Match HST (=Tphot output) Vs IRAC high resolution 
        max_DIST1_D=.5
        MATCH1D=Match_cat(cat_I1_HRES.ALPHA_J2000,cat_I1_HRES.DELTA_J2000,X_WORLD,Y_WORLD,max_DIST1_D) 
        
        if INTERACT=='I': 
            # 1) Flux Low-resolution IRAC SExtractor Vs T-phot
            plt.plot(((23.9-cat_I1.MAG_AUTO[MATCH1C.IDX1])/2.5),np.log10(TPHOT_I1_FLUX[MATCH1C.IDX2]), color='red',marker='.',linestyle='none')
            plt.plot(np.array([-10,10]),np.array([-10,10]),color='blue')
            axes = plt.gca()
            axes.set_xlim([-1,5])
            axes.set_ylim([-1,5])
            plt.ylabel('I1 log(Flux[uJy]) HRES - Tphot')
            plt.xlabel('I1 log(Flux[uJy]) LRES')
            plt.show()

            # 1) Mag Low-resolution IRAC SExtractor Vs T-phot
            plt.plot(cat_I1.MAG_AUTO[MATCH1C.IDX1],TPHOT_I1_MAG[MATCH1C.IDX2], color='red',marker='.',linestyle='none')
            plt.plot(np.array([10,100]),np.array([10,100]),color='blue')
            axes = plt.gca()
            axes.set_xlim([10,30])
            axes.set_ylim([10,30])
            plt.ylabel('I1 log(Flux[uJy]) HRES - Tphot')
            plt.xlabel('I1 log(Flux[uJy]) LRES')
            plt.show()
    
            # 1B) AB magnitude uncertainty Low-resolution IRAC SExtractor Vs T-phot
            plt.plot(np.array([-10,10]),np.array([-10,10]),color='blue')
            plt.plot(cat_I1.MAGERR_AUTO[MATCH1C.IDX1],TPHOT_I1_MAGERR[MATCH1C.IDX2], color='red',marker='.',linestyle='none')
            axes = plt.gca()
            axes.set_xlim([-0.1,0.4])
            axes.set_ylim([-0.1,0.4])
            plt.ylabel('I1 mag uncertainty HRES - Tphot')
            plt.xlabel('I1 mag uncertainty LRES')
            plt.show()

            # 1B) Flux uncertainty Low-resolution IRAC SExtractor Vs T-phot
            FLUXERR_LRES=-(10**((23.9-cat_I1.MAG_AUTO)/2.5))  *((10**(-cat_I1.MAGERR_AUTO/2.5))-1)
            FLUXERR_LRES[FLUXERR_LRES<=0]=-99.
            plt.plot(np.log10(FLUXERR_LRES[MATCH1C.IDX1]),np.log10(TPHOT_I1_FLUXERR[MATCH1C.IDX2]), color='red',marker='.',linestyle='none')
            plt.plot(np.array([-10,10]),np.array([-10,10]),color='blue')
            axes = plt.gca()
            axes.set_xlim([-1,1])
            axes.set_ylim([-1,1])
            plt.ylabel('I1 log(Fluxerr[uJy]) HRES - Tphot')
            plt.xlabel('I1 log(Fluxerr[uJy]) LRES')
            plt.show()
        
        
        
            plt.plot(np.log10(TPHOT_I1_FLUX[MATCH1C.IDX2]), TPHOT_I1_FLUXERR[MATCH1C.IDX2]/TPHOT_I1_FLUX[MATCH1C.IDX2], color='red',marker='.',linestyle='none')
            plt.plot(np.log10(TPHOT_I1_FLUX[MATCH1C.IDX2]), FLUXERR_LRES[MATCH1C.IDX1]/TPHOT_I1_FLUX[MATCH1C.IDX2], color='blue',marker='.',linestyle='none')

            axes = plt.gca()
            # axes.set_xlim([-1,5])
            axes.set_ylim([0,0.4])
            plt.xlabel('I1 log(Flux[uJy]) HRES - Tphot')
            plt.ylabel('I1 Fluxerr/Flux  ')
            plt.show()


            # 2) High-resolution IRAC SExtractor Vs T-phot
            plt.plot(np.array([-10,10]),np.array([-10,10]),color='blue')
            plt.plot(((23.9-cat_I1_HRES.MAG_AUTO[MATCH1D.IDX1])/2.5),np.log10(TPHOT_I1_FLUX[MATCH1D.IDX2]), color='green',marker='.',linestyle='none')
            axes = plt.gca()
            axes.set_xlim([-1,5])
            axes.set_ylim([-1,5])
            plt.ylabel('I1 log(Flux[uJy]) HRES - Tphot')
            plt.xlabel('I1 log(Flux[uJy]) HRES - SExtractor')
            plt.show()

            
    ##############################################################
    # Creating folder.tar.gz containing main outputs 
    ##############################################################
    # folder name specified in the "files" function.
    print 'creating output folder: '+paths.pwd+files.output_folder
    make_output_dir(DS,paths,files,PARN)

    ##############################################################
    # Cleaning working folder from files saved elsewere 
    ##############################################################
    print 'Removing files already saved elsewere: '
    clean_after2(DS,paths,files,PARN)


    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    print "IRAC_Photometry terminated without interruptions"
    print "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"


if __name__ == '__main__':
    main()

