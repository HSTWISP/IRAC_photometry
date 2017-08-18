# List of functions:
#
# gaussfit       # Simple mono-dimensional gaussian fit
#------------------------------------------------------------
# image_depth    # measure depth (1sigma) of an image.fits
                 # using the random apertures method. 
#------------------------------------------------------------
# match          # matches a series of (integer) IDs
#------------------------------------------------------------
# mk_regionfile  # creates simple ds9 region files. 
#------------------------------------------------------------
# image_shift    # shift the RA, DEC center of an image
#------------------------------------------------------------

# used in gaussfit
from scipy.optimize import curve_fit
# used in img_depth:
from photutils import CircularAperture
from photutils import aperture_photometry
import matplotlib.pyplot as plt
# used in mk_regionfile and image_shift
import numpy as np
# used in image_shift:
from astropy.io import fits # to read fits files
import math
import os
# For test purposes:
from pdb import set_trace as stop


def gaussfit(xx,yy,xx0=None,yy0=None,sigma0=None):
    # ------------------------------------------------------------
    # Simple mono-dimensional gaussian fit.
    # ------------------------------------------------------------
    # xx= array of x values
    # yy= array of y values
    # xx0=initial estimate of gaussian mean x (optional)
    # yy0=initial estimate of gaussian height  (optional)
    # sigma0=initial estimate of gaussian sigma (optional)
    #
    # gaussfit returns a numpy.array A containing (Same order as IDL gaussfit)
    # A[0]=height of the Gaussian
    # A[1]=center of the Gaussian
    # A[2]=width (the standard deviation) of the Gaussian
    #################################################################
    xx=np.array(xx)
    yy=np.array(yy)
    # definition of the gaussian function used (3 terms only)
    def gauss(x, x0, y0, sigma):
        p = [x0, y0, sigma]
        return p[1]* np.exp(-((x-p[0])/p[2])**2)
    
    p0A=np.median(xx) # median of the x values to get x0 initial estimate
    p0B=np.median(yy[yy>=np.percentile(yy,80)]) # median of y values higer than the 80% percentile
    p0C_1=np.percentile(xx,16) # 16% percentile on x values
    p0C_2=np.percentile(xx,84) # 84% percentile on x values
    p0C=(p0C_2-p0C_1)/2.
    
    # If the parameters are set by the user:
    if xx0!=None:
        p0A=xx0
    if yy0!=None:
        p0B=yy0
    if sigma0!=None:
        p0C=sigma0

    print "gaussfit parameters initial estimates:", p0B, p0A, p0C
    p0 = [p0A, p0B, p0C] #initialization parameter
    fit, tmp = curve_fit(gauss, xx, yy, p0=p0)
    print "gaussfit output parameters initial estimate:",fit[1],fit[0],fit[2]
    
    return np.array([fit[1],fit[0],fit[2]]) # NOTE: SAME ORDER AS IDL GAUSSFIT




def image_depth(imagename,radius,xmin=None,xmax=None,ymin=None,ymax=None,low_percentile=None,high_percentile=None,INTERACT=None,title=None):
    # ------------------------------------------------------------
    # measure depth (1sigma) of a fits image using the random apertures method. 
    # ------------------------------------------------------------
    # imagename=input image name
    # radius= aperture radius to be used [pixels]
    # xmin= minimum x position for an aperture (optional)
    # xmax= maximum x position for an aperture (optional)
    # ymin= minimum y position for an aperture (optional)
    # ymax= maximum y position for an aperture (optional)
    # low_percentile= minimum percentile of aperture flux (default=1%)
    # high_percentile maximum percentile of aperture flux (default=80%)
    # INTERACT= if set to yes, plots will be shown
    # -----
    # The function creates a grid of apertures evenly distributed 
    # inside the image, or in the range specified using the 
    # xmin,xmax,ymin and ymax input parameters.
    # Using the low_percentile and high_percentile parameters
    # it is possible to exclude extreme values, making the fit of the
    # aperture fluxes histogram more precise. The default values are 
    # 1 and 70 (1% and 70%). The asimmetry allows us to exclude the high 
    # fluxes tail of the distribution, that is not used in the fit.
    # NOTES:
    #   - The image must be background subtracted; local background 
    #     can bring to wrong results;
    #   - Te depth is computed in the same units of the input image.
    # Example
    # from classes.Ivano_functions_1_1 import image_depth
    # from astropy.io import fits
    # imagename='../DATA_TEST_DIR/working_dir/WISPS322_ch1_mosaic.resamp_bck_sub.fits'
    # img=fits.open(imagename)
    # depth=image_depth(imagename,radius=15,INTERACT='I')
    #####################################################
    # load image (Hres)
    img=fits.open(imagename)
    #img.info()
    image=img[0].data # Image data are stored in this np array
    img.close()

    min_Nap=25 # Minimum number of apertures to compute depth
    
    if xmin==None: xmin=0
    if xmax==None: xmax=len(image[0,:])
    if ymin==None: ymin=0
    if ymax==None: ymax=len(image[:,0])
    if title==None: title='Depth computation'
    if low_percentile==None: low_percentile=2
    if high_percentile==None: high_percentile=85
    xmin=round(xmin)
    xmax=round(xmax)
    ymin=round(ymin)
    ymax=round(ymax)
    radius=round(radius)
    
    # N_ap_x= long(round(( (xmax-xmin)/(2.*radius) )-1)) # Number of apertures along x axis
    # N_ap_y= long(round(( (ymax-ymin)/(2.*radius) )-1)) # Number of apertures along x axis
    
    #positions = [( xmin+radius , ymin+radius )]

    Nap=0L #Number of apertures

    YY=ymin+radius #Aperture number
    while YY < ymax:
        XX=xmin+radius #Aperture number
        while XX+radius < xmax:
            if Nap==0 : positions =[( XX,YY )]
            if Nap!=0 : positions.append( (XX,YY) )
            Nap=Nap+1
            XX=XX+(radius*2.)
        YY=YY+(radius*2.)

    apertures = CircularAperture(positions, r=radius)
    phot_table = aperture_photometry(image, apertures)
    # print phot_table
    # print phot_table.keys()
    Ap_Fluxes=np.array(phot_table['aperture_sum']) # Aperture fluxes in a numpy array

    if Nap < min_Nap:
        print 'img_depth: Too few valid apertures to compute depth.'
        print 'Increase area or reduce aperture radius'
        return -1

    if Nap >= min_Nap:
        # create histogram of aperture fluxes
        min_hist=np.percentile(Ap_Fluxes,low_percentile) # 2% flux percentile
        max_hist=np.percentile(Ap_Fluxes,high_percentile)# 85% flux percentile
        # create array of bins
        NNN=2. # increase this number to increase the number of bins
        Nbin=long(NNN*np.sqrt(float(Nap)))
        step_bin=(max_hist-min_hist)/Nbin
        BINS=np.zeros(Nbin)
        BINS[0]=min_hist
        # print 'JJJJJJJJJJJJ' , min_hist
        BB=1l
        while BB<Nbin:
            BINS[BB]=BINS[BB-1]+step_bin
            BB=BB+1
        # print BINS
        
        F_distrib=Ap_Fluxes[(Ap_Fluxes>min_hist) & (Ap_Fluxes<max_hist)]
        if INTERACT=='I':
            plt.hist(F_distrib, bins=BINS) #'auto')
        # Histogram of aperture fluxes:
        HIST=np.histogram(F_distrib, bins=BINS) #'auto')
        HIST[1][:]=HIST[1][:]+((HIST[1][len(HIST[1])-1]-HIST[1][0])/len(HIST[1]))/2.
        # bins with 85% of the highest histogram y value
        tip=0.75*max(HIST[0])
        # tip=np.percentile(HIST[0],80) # influenced by high fluxes (higher threshold needed)
        if INTERACT=='I':
            plt.plot([min_hist,max_hist],[tip,tip],color='red',linestyle='dashed')
        high_idx=np.where(HIST[0] > tip)
        high_idx=high_idx[0]
        # Average x value of the y-highest histograms bins
        val_sym=np.median(HIST[1][high_idx])
        if INTERACT=='I':
            plt.plot(np.array([val_sym,val_sym]),np.array([0,1.1*np.max(HIST[0])]),color='red',linestyle='dashed')
            
        # Create new histogram (with symmetrized distribution  
        IDX=np.where(HIST[1]<=val_sym)
        IDX=IDX[0]
        RR=0L
        while RR < (2*len(IDX))+1:
            if RR==0:
                X_HIST=np.array([HIST[1][0]])
                Y_HIST=np.array([HIST[0][0]])
            if (RR<len(IDX)) and (RR!=0): 
                X_HIST=np.append(X_HIST,HIST[1][RR])
                Y_HIST=np.append(Y_HIST,HIST[0][RR])
#            if (RR>=len(IDX)) and (len(IDX)+(len(IDX)-RR)>=0):
            if (RR>=len(IDX)) and (len(IDX)+(len(IDX)-RR)-1>=0):
#                if RR==92: stop()
                if RR<len(HIST[1][:]):
                    X_HIST=np.append(X_HIST,HIST[1][RR])
                    Val=HIST[1][RR] # maximum x value reached
                if RR>=len(HIST[1][:]):
                    Stp=X_HIST[RR-1]-X_HIST[RR-2]
                    Val=Val+Stp
#                    print Val,Stp
                    X_HIST=np.append(X_HIST,Val)
#                Y_HIST=np.append(Y_HIST,HIST[0][len(IDX)+(len(IDX)-RR)])
                Y_HIST=np.append(Y_HIST,HIST[0][len(IDX)+(len(IDX)-RR)-1])
            RR=RR+1
        NEW_HIST=(Y_HIST,X_HIST)
        plt.plot(X_HIST,Y_HIST)
        
        
        # FIT OF A GAUSSIAN FUNCTION
        gfit=gaussfit(X_HIST,Y_HIST)
        # SHOW fit:
        
        if INTERACT=='I':
            def gauss(x, x0, y0, sigma):
                p = [x0, y0, sigma]
                return p[1]* np.exp(-((x-p[0])/p[2])**2)
            plt.plot(X_HIST, gauss(X_HIST, gfit[1], gfit[0], gfit[2]), 'b-')
            plt.ylabel('N apertures')
            plt.xlabel('Aperture total Counts')
            # plt.title='Counts distribution and fit'
            plt.title(title)
            plt.show()
            
        return gfit[2] # Depth 1sigma


    # from classes.Ivano_functions_1_1 import image_depth
    # from astropy.io import fits
    # imagename='../DATA_TEST_DIR/working_dir/WISPS322_ch1_mosaic.resamp_bck_sub.fits'
    # img=fits.open(imagename)
    # depth=image_depth(imagename,radius=32,INTERACT='I')
    # depth1=depth*3.0*8.461595
    # import numpy as np
    # depth_mag=-((2.5*np.log10(depth1))-23.9)
    # print depth_mag


def mk_regionfile(RA,DEC,RADIUS,file=None,width=None,color=None,font=None,label=None):
    #---------------------------------------------------
    # Function to create simple ds9 region files. 
    #---------------------------------------------------
    # RA =numpy array of RA coordinates in deg.decimals
    # DEC=numpy array of DEC coordinates in degdecimals
    # RADIUS=radius, in pixels of the output region files
    # color=color of the output region files
    # width= integer specifying the width of the output circles. 
    # font=string specifying the output style. 
    #      Example: font='helvetica 14 bold'
    # label=input -string- or -string array- or -float array-  
    #       converted to a string in the output file.
    #####################################################

    #----------------------------------------------------
    # Set defaults
    #----------------------------------------------------
    if file==None:
        file='ds9.reg'
    if color==None:
        color='red'
    #----------------------------------------------------
    RADIUS_str=str(RADIUS)

    WIDTH_STR=''
    if width !=None:
        if type(width)==str:
            WIDTH_STR=" width="+width
        if type(width)!=str:
            WIDTH_STR=" width=%i" % width

    FONT_STR=''
    if font !=None:
        FONT_STR=" font='"+font+"'"
        
    LABEL_STR=''
    if label ==None:
        ALLDATA = np.zeros(RA.size, dtype=[('RA', 'float'), ('DEC', float)])
        ALLDATA['RA'] = RA
        ALLDATA['DEC'] = DEC

    if label !=None:
        LABEL_STR=" text={%s}"
        ALLDATA = np.zeros(RA.size, dtype=[('RA', 'float'), ('DEC', float), ('label', "S"+str(len(label))) ])
        ALLDATA['RA'] = RA
        ALLDATA['DEC'] = DEC
        ALLDATA['label'] =label
    
    ALLDATA=ALLDATA.T # TRANSPOSE!!!
    STRING="circle %15.10fd %15.10fd "+RADIUS_str+" # "+WIDTH_STR+" color="+color+FONT_STR+LABEL_STR
    #----------------------------------------------------
    # Write region file
    #----------------------------------------------------
    np.savetxt(file,ALLDATA, fmt=STRING)


def image_shift(image_in,image_out=None,header_key_RA=None,header_key_DEC=None,delta_RA=None,delta_DEC=None,dist_RA=None,dist_DEC=None,new_RA=None,new_DEC=None):
    #------------------------------------------------------------
    #  image_in=image to be shifted
    #  image_out=[OPTIONAL] where to save the output (default: image_in
    #                is overwritten)
    #  header_key_RA=[OPTIONAL] 'header keyword specifying
    #                the reference pixel RA coordinate (default='CRVAL1') 
    #  header_key_DEC=[OPTIONAL] 'header keyword specifying 
    #                the reference pixel DEC coordinate (default='CRVAL2')
    #  ONE OF THE FOLLOWING OPTIONS MUST BE SPECIFIED: 
    #  FIRST OPTION ----------------------------------------------
    #  delta_RA= [OPTIONAL] to add to RA of the reference pixel  [arcseconds]
    #  delta_DEC=[OPTIONAL] to add to DEC of the reference pixel [arcseconds]
    #  SECOND OPTION ---------------------------------------------
    #  dist_RA= [OPTIONAL] to add ALONG the RA direction   [arcseconds]
    #  dist_DEC=[OPTIONAL] to add ALONG the DEC directoion [arcseconds]
    #  THIRD OPTION ----------------------------------------------
    #  new_RA= [OPTIONAL] # New RA center  [deg.decimals]
    #  new_DEC=[OPTIONAL] # New DEC center [deg.decimals]
    #  -----------------------------------------------------------
    # 
    # NOTES on defaults
    # 1) if the output image name is not set, the input
    #    image is overwritten.
    # 2) If the header keywords specifying the image refernce
    #    pixels are not specified in header_key_RA and header_key_DEC,
    #    header_key_RA is set to "CRVAL1" and header_key_DEC to "CRVAL2"
    # 4) Note that delta_RA is not the "angular distance along RA",
    #    but the actual delta_RA=(angular distance along RA)/cos(dec)
    #----------------------------------------------------
    if header_key_RA==None: header_key_RA='CRVAL1'
    if header_key_DEC==None: header_key_DEC='CRVAL2'
    if image_out==None: image_out=image_in
    #----------------------------------------------------
    # Read image
    img = fits.open(image_in)
    # Read image reference pixel from header
    RA_ref=img[0].header[header_key_RA]
    DEC_ref=img[0].header[header_key_DEC]

    ERR_MSG=' -------------------------- image_shift() --------------------------\n'+ \
    'ERROR: the operation to perform is unclear\n'+ \
    'Please specify only one of the following pair of inputs:\n'+ \
    '(delta_RA,delta_DEC) or (dist_RA,dist_DEC) or (new_RA,new_DEC)\n'+ \
    "This program is stopped here\n"+ \
    ' -------------------------------------------------------------------'
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    #  FIRST OPTION ----------------------------------------------
    if delta_RA!=None and delta_DEC!=None:
        if dist_RA!=None or dist_DEC!=None or new_RA!=None or new_DEC!=None:
            print ERR_MSG
            exit()
        # Compute new reference coordinates
        new_DEC_ref=DEC_ref+delta_RA/3600.
        new_RA_ref=RA_ref+delta_DEC/3600.

    #  SECOND OPTION ----------------------------------------------
    if dist_RA!=None and dist_DEC!=None:
        if delta_RA!=None or delta_DEC!=None or new_RA!=None or new_DEC!=None:
            print ERR_MSG
            exit()
        # Compute new reference coordinates
        new_DEC_ref=DEC_ref+dist_DEC/3600.
        new_RA_ref=RA_ref+(dist_RA/np.cos(math.pi*new_DEC_ref/180.))/3600.

    #  THIRD OPTION ----------------------------------------------
    if new_RA!=None and new_DEC!=None:
        if delta_RA!=None or delta_DEC!=None or dist_RA!=None or dist_DEC!=None:
            print ERR_MSG
            exit()
        # Compute new reference coordinates
        new_DEC_ref=new_DEC
        new_RA_ref=new_RA
    #xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
    # WRITE NEW KEYWORDS IN THE HEADER
    img[0].header[header_key_RA]=new_RA_ref
    img[0].header[header_key_DEC]=new_DEC_ref

    # WRITE OUTPUT IMAGE
    if image_in==image_out:
        print image_in
        os.system('rm '+image_in) # Remove input image if output not specified.
    img.writeto(image_out)  # Save changes into output file (default: overwrite input) 
    img.close()   # Close file

