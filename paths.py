class Get_paths(object):
    #---------------------------------------------------
    # Paths definitions
    #---------------------------------------------------
    def __init__(self, TS):
        self.TS=TS  # 1= Test run, 0=normal run

        #################################################################################################
        # REAL RUN (set the first 4 paths as required) 
        #################################################################################################
        if self.TS ==0:

            # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            # THESE PATHS SHOULD BE MODIFIED //////////////////////////////////////
            self.IRAC_DATA_PATH_I1='/Volumes/PROMISE_PEGASUS/IRAC_DATA/Data/CH1/' # IRAC 1 ORIGINAL DATA PATH
            self.IRAC_DATA_PATH_I2='/Volumes/PROMISE_PEGASUS/IRAC_DATA/Data/CH2/' # IRAC 2 ORIGINAL DATA PATH
            self.IRAC_PHOT_FOLDER='/Users/ivano/IRAC_photometry/'                 # location of IRAC_photometry.py
            self.TPHOT_LOCATION='/Users/ivano/Software/TPHOT/tphot/bin/'           # PATH to the TPHOT executable
            # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
            # @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

            # THE FOLLOWING PATHS SHOULD - NOT - BE MODIFIED  //////////////////////

            self.pwd='DATA/IRAC/'                        # Working directory for SWarp and tphot. Created by this program. Do not modify.
            self.p0='DATA/DIRECT_GRISM/'                 # J and H image & cat location (DATA/DIRECT_GRISM/)
            self.p1=self.IRAC_DATA_PATH_I1               # IRAC1 original image location. Set absolute path!   
            self.p2=self.IRAC_DATA_PATH_I2               # IRAC2 original image location. Set absolute path!               
            self.p3=self.pwd+'IRAC_SEx_cat/'             # output IRAC SExtractor catalogs location     
            self.p4=self.IRAC_PHOT_FOLDER+'SEx_config/'  # folder with SExtractor config files for this program 
                                                          # (do not confuse with "SEX" folder created by the
                                                          # wisp pipeline)
            self.p5='SEX/'                                # wisps pipeline "SEX" folder, containing segmentation maps.
            self.p6=self.IRAC_PHOT_FOLDER+'SWarp_config/' # SWARP configuration files location
            self.p7=self.pwd                              # SWARP orientation header location
            self.p8=self.IRAC_PHOT_FOLDER+'TPHOT_config/' # Tphot configuration files folder
            self.p9=self.TPHOT_LOCATION                        # location of "tphot" executable. Set absolute path!

        #################################################################################################
        # TEST RUN (set these pats to test folders if needed. Used only for debugging purposes.
        #################################################################################################
        if self.TS ==1:
            self.path_test='../DATA_TEST_DIR/'         # Defined ONLY for test purposes
            self.pwd='../DATA_TEST_DIR/working_dir/'   # Working directory for SWarp and tphot. Created by this program.
            self.p0='../DATA_TEST_DIR/DIRECT_GRISM/'   # J and H image & cat location (DATA/DIRECT_GRISM/)
            self.p1='../DATA_TEST_DIR/IRAC_images/CH1/'# IRAC1 original image location. Set absolute path!
            self.p2='../DATA_TEST_DIR/IRAC_images/CH2/'# IRAC2 original image location. Set absolute path!
            self.p3=self.pwd+'/IRAC_SEx_cat/'              # output IRAC SExtractor catalogs location
            self.p4='SEx_config/'                      # folder with SExtractor config files for this program 
                                                       # (do not confuse with "SEX" folder created by the
                                                       # wisp pipeline)
            self.p5='../DATA_TEST_DIR/SEX/'            # wisps pipeline "SEX" folder, containing segmentation. maps.
            self.p6='SWarp_config/'                    # SWARP configuration files location
            self.p7=self.pwd                           # SWARP orientation header location
            self.p8='TPHOT_config/'                    # Tphot configuration files folder
            self.p9='/home/baronchelli/software1/TPHOT/T-PHOT_2.0/tphot/bin/' # location of "tphot" executable.
                                                                              # Set absolute path!
        #################################################################################################

