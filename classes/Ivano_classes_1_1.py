# List of classes:
#----------------------------------------------------
# Deg2sex   # Given RA and DEC in input, they will be returnet in
#           # output converted from 
#           # decimal to sexagesimal or vice-versa (depending on
#           # intype and outype variables set).
#----------------------------------------------------
# Match     # matches a series of (integer) input IDs (like IDL match)
#----------------------------------------------------
# Match_cat # Matches coordinates [deg], in two lists of coordinates, 
#           # using an user specyfied searching radius dt[arcsec]
#           # For each source in the first catalog, the closest
#           # counterpart in the second catalog is associated.
#----------------------------------------------------
# Readcol # Read specified columns in an ascii file
#         # and returns them in the specified format
#----------------------------------------------------
####################################################################
# MODULES
####################################################################
### Modules for Readcol:
from pdb import set_trace as stop # for test purposes
import os
import numpy as np
### Modules for Match_cat:
# import numpy as np (already imported
from astropy.coordinates import ICRS
from astropy import units as u
import astropy.coordinates.representation
import math
import matplotlib.pyplot as plt # to plot

class Deg2sex(object):
    ##################################################################################
    # Given RA and DEC in input, they will be returnet in output converted from 
    # decimal to sexagesimal or vice-versa (depending on intype and outype variables.
    # The number of decimals after the seconds (ss.decimals) is set using ndec (default=2)
    #
    # LIST OF FUNCTIONS IMPLEMENTED SO FAR:
    # intype=1, outype=2
    # DEC) deg.decimals (float) --> deg:mm:ss.decimals (string format)
    # RA)  deg.decimals (float) --> hh:mm:ss.decimals  (string format)
    #
    # intype=outype --> input RA and DEC returned in the identical input form 
    #
    # Examples:
    # coord=Deg2sex(130.58137, 14.153013)
    # print coord.DEC(1,2,ndec=1)
    # > 14:09:10.8
    # print coord.RA(1,2,ndec=3)
    # > 08:42:19.529
    # print coord.RA(1,1,ndec=3)
    # > 130.58137
    ##################################################################################
    def __init__(self,RA_IN,DEC_IN):
        self.RA_IN=RA_IN
        self.DEC_IN=DEC_IN
    def DEC(self,intype,outype,ndec=2):
        self.intype=intype
        self.outype=outype
        self.ndec=ndec
        strdec="%."+str(ndec)+"f"
        if self.intype == outype:
            return self.DEC_IN
        # decimal -- > SEXAGESIMAL
        if self.intype == 1 and self.outype == 2:
            # returns the modulo division (does not work properly. Correction below)
            # dec1_int,dec1_dec=divmod(self.DEC_IN,1)
            # correction:---------------------
            dec1_int=abs(float(long(self.DEC_IN))) # Absolute - Integer part
            dec1_dec=abs(abs(self.DEC_IN)-dec1_int) # Absolute - decimal values
            #---------------------------------
            dec1_int_str=str(abs(dec1_int)) # transform into a string type - NO sign
            dec1_int_str=dec1_int_str.split(".")[0] # get what stays before "."
            if self.DEC_IN >= 0:
                if dec1_int >= 10: ref_DEC_center_sex=dec1_int_str
                if dec1_int < 10 : ref_DEC_center_sex='0'+dec1_int_str
            if self.DEC_IN < 0:
                if dec1_int >= 10: ref_DEC_center_sex='-'+dec1_int_str
                if dec1_int < 10 : ref_DEC_center_sex='-0'+dec1_int_str
            # dec1_int_A,dec1_dec_A=divmod(dec1_dec*60.,1)
            # correction:---------------------
            dec1_int_A=float(long(dec1_dec*60.))
            dec1_dec_A=dec1_dec*60. - dec1_int_A
            #---------------------------------
            dec1_int_A_str=str(dec1_int_A)
            dec1_int_A_str=dec1_int_A_str.split(".")[0]
            if dec1_int_A >= 10: ref_DEC_center_sex=ref_DEC_center_sex+':'+dec1_int_A_str
            if dec1_int_A < 10: ref_DEC_center_sex=ref_DEC_center_sex+':0'+dec1_int_A_str
            dec1_dec_A60=dec1_dec_A*60.
            # dec1_dec_A60=(strdec % dec1_dec_A60)
            # dec1_B_str=str(dec1_dec_A60)
            dec1_B_str=str(strdec % dec1_dec_A60)
            if dec1_dec_A60 >= 10: ref_DEC_center_sex=ref_DEC_center_sex+':'+dec1_B_str
            if dec1_dec_A60 < 10 : ref_DEC_center_sex=ref_DEC_center_sex+':0'+dec1_B_str
            return ref_DEC_center_sex
        return self.DEC
    
    def RA(self,intype,outype,ndec=2):
        self.intype=intype
        self.outype=outype
        self.ndec=ndec
        strdec="%."+str(ndec)+"f"
        if self.intype == outype:
            return self.RA_IN
        # decimal -- > SEXAGESIMAL
        if self.intype == 1 and self.outype == 2:
            RA_hr=24.*self.RA_IN/360.
            # returns the modulo division (does not work properly. Correction below)
            # ra1_int,ra1_dec=divmod(RA_hr,1)
            # correction:---------------------
            ra1_int=abs(float(long(RA_hr))) # Absolute - Integer part
            ra1_dec=abs(abs(RA_hr)-ra1_int) # Absolute - decimal values
            #---------------------------------
            ra1_int_str=str(abs(ra1_int)) # transform into a string type - NO sign
            ra1_int_str=ra1_int_str.split(".")[0] # get what stays before "."
            if self.RA_IN>=0:
                if ra1_int >= 10: ref_RA_center_sex=ra1_int_str
                if ra1_int < 10 : ref_RA_center_sex='0'+ra1_int_str
            if self.RA_IN<0:
                if ra1_int >= 10: ref_RA_center_sex='-'+ra1_int_str
                if ra1_int < 10 : ref_RA_center_sex='-0'+ra1_int_str
                
            # ra1_int_A,ra1_dec_A=divmod(ra1_dec*60.,1)
            # correction:---------------------
            ra1_int_A=float(long(ra1_dec*60.))
            ra1_dec_A=ra1_dec*60. - ra1_int_A
            #---------------------------------
            ra1_int_A_str=str(ra1_int_A)
            ra1_int_A_str=ra1_int_A_str.split(".")[0]
            if ra1_int_A >= 10: ref_RA_center_sex=ref_RA_center_sex+':'+ra1_int_A_str
            if ra1_int_A < 10 : ref_RA_center_sex=ref_RA_center_sex+':0'+ra1_int_A_str
            ra1_dec_A60=ra1_dec_A*60.
            # ra1_dec_A60=(strdec % ra1_dec_A60)
            # ra1_B_str=str(ra1_dec_A60)
            ra1_B_str=str(strdec % ra1_dec_A60)
            if ra1_dec_A60 >= 10: ref_RA_center_sex=ref_RA_center_sex+':'+ra1_B_str
            if ra1_dec_A60 < 10 : ref_RA_center_sex=ref_RA_center_sex+':0'+ra1_B_str
            return ref_RA_center_sex

# from classes.Ivano_classes_1_1 import Deg2sex
# coord=Deg2sex(1.753692550119E+01, -2.423583178468E+00)
# print coord.RA(1,2,ndec=3),coord.DEC(1,2,ndec=3)

class Match(object):
################################################
# Similar to the IDL function "match.pro".
# It matches the IDs (integer numbers) stored inside two input 
# numpy arrays or lists (ID1 and ID2).
# The indexes of the matching IDs are stored into two output 
# numpy arrays (IDX_1_2 and IDX_2_1)
# Example:
# > a=[1,3,5,7,9]
# > b=[5,6,7,8,9,10]
# > ID_MATCH=Match(A,B)
# > print ID_MATCH.IDX_1_2
# [2,3,4]
# > print ID_MATCH.IDX_2_1
# [0,2,4]
# 
# NOTE:
# - if the IDs in the input arrays are not unique, the output
# indexes will consider all the associations found.
# Example: AA=[10,11,12,10], BB=[70,50,30,10,20,10]
#          will give: IDX_1_2=[0,0,3,3] and IDX_2_1=[3,5,3,5]
# - non integer inputs will be rounded to the closest integer value
#   using the python approximation standard (that for 2.5 is not 3,
#   as they thought me at the elementary school, but 2. Fuck off python!).
################################################
    
    def __init__(self,ID1_in,ID2_in):
        ID1_in=np.around(np.array(ID1_in)) # round input arrays using python standard
        ID2_in=np.around(np.array(ID2_in)) # round input arrays using python standard
        self.ID1_in=np.array([long(x) for x in ID1_in]) # "type" is now "long" for all the elements
        self.ID2_in=np.array([long(x) for x in ID2_in]) # "type" is now "long" for all the elements
        
        if len(ID1_in)<=len(ID2_in):
            ID1=ID1_in
            ID2=ID2_in
        if len(ID2_in)<len(ID1_in): # Invert order to save number of cycles
            ID1=ID2_in
            ID2=ID1_in
            
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@        
        # output indexes:
        IDX_1_2=np.array([-1])
        IDX_2_1=np.array([-1])
        FF=0L
        while FF<len(ID1):
            IDX_A=np.where(ID2==ID1[FF]) # search for an association
            IDX_A=IDX_A[0]
            if len(IDX_A)>=1: # association found
        #---------------------------------------------------
                if IDX_1_2[0]==-1: # no previous associations found
                    IDX_2_1=np.array([IDX_A])
                    IDX_1_2=np.array([FF])
                    ii=0l
                    while ii <len(IDX_A)-1:
                        IDX_1_2=np.append(IDX_1_2,FF)
                        ii=ii+1
        #---------------------------------------------------
                #if IDX_1_2[0]!=-1: # previous associations found
                else:
                    IDX_2_1=np.append(IDX_2_1,[IDX_A])
                    ii=0l
                    while ii <len(IDX_A):
                        IDX_1_2=np.append(IDX_1_2,[FF])
                        ii=ii+1
            FF=FF+1
#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
        if len(ID1_in)<=len(ID2_in):
            self.IDX_1_2=IDX_1_2
            self.IDX_2_1=IDX_2_1
        if len(ID2_in)<len(ID1_in): # Invert order to save number of cycles
            self.IDX_1_2=IDX_2_1
            self.IDX_2_1=IDX_1_2

    def IDX_1_2(self):
        return self.IDX_1_2
    def IDX_2_1(self):
        return self.IDX_2_1


class Match_cat(object):
################################################
# Match coordinates [deg], in two lists of coordinates, 
# using an user specyfied searching radius dt[arcsec]
# For each source in the first catalog, the closest
# counterpart in the second catalog is associated.
# Some functions allow you to plot the results of the match
# Input parameters are numpy single dimension arrays.
# Example:
# > MATCH1=Match_cat(ra1,dec1,ra2,dec2,dist)
# > idx1=MATCH1.IDX1 # indexes of first catalog having a counterpart in the second 
# > idx2=MATCH1.IDX2 # indexes of second catalog that have a conterpart in the first
# > DIST=MATCH1.DIST # [arcsec] distance between counterparts (always<dist)
# > DIST_RA=MATCH1.DIST_RA # [arcsec] distance between counterparts (always<dist) along RA axis
# > DIST_DEC=MATCH1.DIST_DEC # [arcsec] distance between counterparts (always<dist) along DEC axis
# > MATCH1.plot_RA_DEC()  # plot Ra Vs Dec
# > MATCH1.plot_delta_RA_DEC()  # plot Delta Ra Vs Delta Dec
# > MATCH1.plot_delta2_RA_DEC() # plot Delta Ra and Delta Dec as a function of Ra and dec (4 plots in a window)
################################################
    # Class definition
    def __init__(self, RA1,DEC1,RA2,DEC2,dist):
        # Input definition
        self.RA1=RA1   # Ra first set of coordinates   (numpy array single dimension)
        self.DEC1=DEC1 # Dec first set of coordinates  (numpy array single dimension)
        self.RA2=RA2   # Ra second set of coordinates  (numpy array single dimension)
        self.DEC2=DEC2 # Dec second set of coordinates (numpy array single dimension)
        self.dist=dist # Maximum correlation distance in arcseconds (float)
        # 1) Transform np aray to "Angle" type
        RA1_angle=astropy.coordinates.representation.Longitude(self.RA1, unit=u.deg)
        DEC1_angle=astropy.coordinates.representation.Latitude(self.DEC1, unit=u.deg)
        RA2_angle=astropy.coordinates.representation.Longitude(self.RA2, unit=u.deg)
        DEC2_angle=astropy.coordinates.representation.Latitude(self.DEC2, unit=u.deg)
        # 2) create catalog classes
        cat1=ICRS(RA1_angle,DEC1_angle)#, unit=(u.degree, u.degree))
        cat2=ICRS(RA2_angle,DEC2_angle)#, unit=(u.degree, u.degree))
        # 3) Match catalogs using astropy.coordinates.match_coordinates_sky"
        d1d, d2d, d3d = astropy.coordinates.match_coordinates_sky(cat1,cat2)
        # Second catalog indexes of sources closer than "dist" ONLY! 
        self.IDX2=d1d[(d2d.degree *3600. < self.dist)]
        # First catalog indexes of sources closer than "dist" ONLY! 
        IDX1ALL = list(range(len(RA1)))
        IDX1ALL = np.array(IDX1ALL)
        self.IDX1=IDX1ALL[(d2d.degree *3600. < self.dist)]
        # Distances of sources closer than "dist" ONLY! 
        self.DIST=3600*d2d.degree[(d2d.degree*3600. < self.dist)]
        # Distance [arcsec] along RA (corrected for cos(dec).... it is not a DELTA_RA!)
        self.DIST_RA=3600.*(self.RA2[self.IDX2]-self.RA1[self.IDX1])*np.cos(math.pi*np.mean(np.array([self.DEC1[self.IDX1],self.DEC2[self.IDX2]])/180.))
        # Distance [arcsec] along DEC
        self.DIST_DEC=3600.*(self.DEC2[self.IDX2]-self.DEC1[self.IDX1])
    def DIST(self):
        return self.DIST
    def IDX1(self):
        return self.IDX1
    def IDX2(self):
        return self.IDX2
    def DIST_RA(self):
        return self.DIST_RA
    def DIST_DEC(self):
        return self.DIST_DEC
    
    def plot_RA_DEC(self,color1=None,marker1=None,markersize1=None,color2=None,marker2=None,markersize2=None,xlabel=None,ylabel=None,title=None):
        # Simple plotting function Ra Vs DEC
        #---------------------------------
        # Set defaults
        #---------------------------------
        if title ==None: title=' '
        if xlabel==None: xlabel='RA'
        if ylabel==None: ylabel='DEC'
        #---------------------------------
        if color1==None: color1='r'
        if marker1==None: marker1='.'
        if markersize1==None:markersize1=15
        #---------------------------------
        if color2==None: color2='b'
        if marker2==None: marker2='.'
        if markersize2==None: markersize2=7
        #---------------------------------
        plt.plot(self.RA1[self.IDX1],self.DEC1[self.IDX1],linestyle='none', color=color1,marker=marker1,markersize=markersize1)
        plt.plot(self.RA2[self.IDX2],self.DEC2[self.IDX2],linestyle='none', color=color2,marker=marker2,markersize=markersize2)
        plt.ylabel(ylabel)
        plt.xlabel(xlabel)
        plt.title(title)
        plt.show()

    def plot_delta_RA_DEC(self,color1=None,marker1=None,markersize1=None,xlabel=None,ylabel=None,crosscolor=None,color=None,marker=None,markersize=None,title=None):
        # Simple plotting function Delta_Ra Vs Delta_DEC
        #---------------------------------
        # Set defaults
        #---------------------------------
        if title ==None: title=' '
        if xlabel==None: xlabel='Distance along RA [arcsec]'
        if ylabel==None: ylabel='Distance along DEC [arcsec]'
        if crosscolor==None:crosscolor='b'
        if color==None: color='g'
        if marker==None: marker='.'
        if markersize==None:markersize=10
        #---------------------------------
        self.dist
        plt.plot(np.array([0,0]),np.array([-self.dist*1.1,self.dist*1.1]), color=crosscolor)
        plt.plot(np.array([-self.dist*1.1,self.dist*1.1]),np.array([0,0]), color=crosscolor)
        plt.plot(self.DIST_RA,self.DIST_DEC,linestyle='none', color=color,marker=marker,markersize=markersize)
        axes = plt.gca()
        axes.set_xlim([-self.dist*1.1,self.dist*1.1])
        axes.set_ylim([-self.dist*1.1,self.dist*1.1])
        plt.ylabel(xlabel)
        plt.xlabel(ylabel)
        plt.title(title)
        plt.show()

        
    def plot_delta2_RA_DEC(self,linecolor1=None,linecolor2=None,linewidth=None,color=None,marker=None,markersize=None,xlabel1=None,ylabel1=None,xlabel2=None,ylabel2=None,xlabel3=None,ylabel3=None,xlabel4=None,ylabel4=None,title=None):
        # plot Delta Ra and Delta Dec as a function of Ra and dec.
        # Plot divided in 4 parts showing Delta ra as a function 
        # of Ra and dec and delta dec as a function of  Ra and dec 
        #---------------------------------
        # Set defaults
        #---------------------------------
        if title ==None: title=' '
        if linecolor1==None: linecolor1='b'
        if linecolor2==None: linecolor2='r'
        if linewidth==None: linewidth=3
        if color==None: color='g'
        if marker==None: marker='.'
        if markersize==None:markersize=10
        if xlabel1==None: xlabel1='RA [arcsec]'
        if ylabel1==None: ylabel1='Distance along DEC [arcsec]'
        if xlabel2==None: xlabel2='DEC [arcsec]'
        if ylabel2==None: ylabel2='Distance along DEC [arcsec]'
        if xlabel3==None: xlabel3='RA [arcsec]'
        if ylabel3==None: ylabel3='Distance along RA [arcsec]'
        if xlabel4==None: xlabel4='DEC [arcsec]'
        if ylabel4==None: ylabel4='Distance along RA [arcsec]'
        #---------------------------------
        min_ra=np.min(self.RA1[self.IDX1])
        max_ra=np.max(self.RA1[self.IDX1])
        min_dec=np.min(self.DEC1[self.IDX1])
        max_dec=np.max(self.DEC1[self.IDX1])
        #---------------------------------
        # first plot 
        plt.subplot(221)
        plt.plot(np.array([min_ra-0.005,max_ra+0.005]),np.array([0,0]), color=linecolor1,linewidth=linewidth)
        plt.plot(np.array([min_ra-0.005,max_ra+0.005]),[np.median(self.DIST_DEC),np.median(self.DIST_DEC)],color=linecolor2,linestyle='dashed',linewidth=linewidth)
        plt.plot(self.RA1[self.IDX1],self.DIST_DEC,linestyle='none',color=color,marker=marker,markersize=markersize)
        plt.ylabel(ylabel1)
        plt.xlabel(xlabel1)
        plt.axis([min_ra-0.005,max_ra+0.005,-self.dist,self.dist])
        ax = plt.gca()
        ax.set_autoscale_on(False)
        #---------------------------------
        # second plot 
        plt.subplot(222)
        plt.plot(np.array([min_dec-0.005,max_dec+0.005]),np.array([0,0]), color=linecolor1,linewidth=linewidth)
        plt.plot(np.array([min_dec-0.005,max_dec+0.005]),[np.median(self.DIST_DEC),np.median(self.DIST_DEC)] ,color=linecolor2,linestyle='dashed',linewidth=linewidth)
        plt.plot(self.DEC1[self.IDX1],self.DIST_DEC,linestyle='none',color=color,marker=marker,markersize=markersize)
        plt.ylabel(ylabel2)
        plt.xlabel(xlabel2)
        plt.axis([min_dec-0.005,max_dec+0.005,-self.dist,self.dist])
        ax = plt.gca()
        ax.set_autoscale_on(False)
        #---------------------------------
        # third plot 
        plt.subplot(223)
        plt.plot(np.array([min_ra-0.005,max_ra+0.005]),np.array([0,0]), color=linecolor1,linewidth=linewidth)
        plt.plot(np.array([min_ra-0.005,max_ra+0.005]),[np.median(self.DIST_RA),np.median(self.DIST_RA)] ,color=linecolor2,linestyle='dashed',linewidth=linewidth)
        plt.plot(self.RA1[self.IDX1],self.DIST_RA,linestyle='none',color=color,marker=marker,markersize=markersize)
        plt.ylabel(ylabel3)
        plt.xlabel(xlabel3)
        plt.axis([min_ra-0.005,max_ra+0.005,-self.dist,self.dist])
        ax = plt.gca()
        ax.set_autoscale_on(False) 
        #---------------------------------
        # fourth plot 
        plt.subplot(224)
        plt.plot(np.array([min_dec-0.005,max_dec+0.005]),np.array([0,0]), color=linecolor1,linewidth=linewidth)
        plt.plot(np.array([min_dec-0.005,max_dec+0.005]),[np.median(self.DIST_RA),np.median(self.DIST_RA)] ,color=linecolor2,linestyle='dashed',linewidth=linewidth)
        plt.plot(self.DEC1[self.IDX1],self.DIST_RA,linestyle='none',color=color,marker=marker,markersize=markersize)
        plt.ylabel(ylabel4)
        plt.xlabel(xlabel4)
        plt.axis([min_dec-0.005,max_dec+0.005,-self.dist,self.dist])
        ax = plt.gca()
        ax.set_autoscale_on(False)
        #---------------------------------
        # Finally... plot
        plt.suptitle(title)
        plt.show()



class Readcol(object):
    ##################################################################################
    # Read specified columns in an ascii file and returns them in the specified format
    # Example of use:
    # cat=Readcol(path0+'fin_F160.cat','f,f,x,x,x,a',skipline=17,skipcol=7)
    # X=cat.col(0)# will be a float numpy array
    # Y=cat.col(1) # will be a float numpy array
    # Z=cat.col(2) # will be a string numpy array
    # # Note that the first 17 lines and the first 7 columns are not considered.
    # # For the columns skipped using "skipcol", the format should not be specified.
    # # The indexes "n" of the output columns, used in col(n), start from 0 and do not 
    # # consider skipped columns (skipcol or 'x').
    ##################################################################################
    def __init__(self, filename, form,skipline=0, skipcol=0, sep='default'):
        self.filename=filename # string containing path+filename
        self.form=form         # format of the elements in the columns. Example:
                               # 'f,f,f,x,a,i' --> first three columns (after the 
                               # skipped ones, specified using skipcol) will be
                               # returned as float, third column is jumped, fourth 
                               # column is returned as a character, fifth as an integer.   
        self.skipline=skipline # number of lines to skip (optional, default=0)
        self.skipcol=skipcol   # number of clumns to skip (optional, default=0)
        self.sep=sep           # Separator. Default are spaces. Other options 
                               # are ',' '\t' (tab). It accepts all the options
                               # allowed in string.split()
        if os.path.isfile(filename)==0:
            print 'File '+ filename +' not found'
        if os.path.isfile(filename)==True:
            FILE=open(filename,'r')
            ALL_COL_STR=FILE.readlines()
            FILE.close()
            
            FORMAT=np.array(form.split(',')) # Must be converted otherwise it doesn't work
            ncol=len(FORMAT[FORMAT != 'x'])  # Number of output columns 
            out_format=['' for x in xrange(ncol)]   # Format of output columns
            out_format=np.array(out_format, np.str) # numpy array of strings
            nlines=len(ALL_COL_STR)-skipline # Number of output lines 

            all_col=[['                              ' for x in xrange(ncol)] for x in xrange(nlines)]
            all_col=np.array(all_col, np.str) # numpy array of strings
            CN=skipcol # Input Column Number (also 'x' considered here)
            RCN=0      # Real (output) Column Number (no 'x')
            while CN < len(FORMAT)+skipcol:
                if FORMAT[CN-skipcol]!='x':
                    LN=skipline # Input line Number
                    RLN=0       # Real (output) line Number (no 'x')
                    while LN < len(ALL_COL_STR):
                        line=ALL_COL_STR[LN]
                        if sep=='default':
                            linesplit = line.split()
                        if sep!='default':
                            linesplit = line.split(self.sep)
                        #--------------------------------
                        all_col[RLN,RCN]=linesplit[CN]
                        #--------------------------------
                        LN=LN+1
                        RLN=RLN+1
                    out_format[RCN]=FORMAT[CN-skipcol]
                    RCN=RCN+1
                CN=CN+1
        self.out_format=out_format
        self.all_col=all_col        
    def col(self,coln):
        ###################################################
        # "coln" corresponds to the column number on the ascii
        # file (start from 0) minus "skipcol", minus the number
        # of "x" indicated in the input format string "form"  
        ###################################################
        if self.out_format[coln].lower()=='a':
            OUTCOL=np.array(self.all_col[:,coln], np.str)
        if self.out_format[coln].lower()=='i':
            OUTCOL=np.array(self.all_col[:,coln], np.int)
        if self.out_format[coln].lower()=='f':
            OUTCOL=np.array(self.all_col[:,coln], np.float)
        return OUTCOL

