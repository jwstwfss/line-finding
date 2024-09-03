import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import glob
import os

def gaussian(x, mu, sigma):
    return np.exp(-np.power(x - mu, 2.) / (2. * np.power(sigma, 2.))) # / (sigma * np.sqrt(2. * np.pi))



def is_number(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def read_config(config, availgrism='both'): 
    configfile = open(config, 'r')  
    config_pars = {} 
    for line in configfile:
        if ( (line[0] != '#') & (len(line)>1)): 
            name = line.split()[0] 
            
            if name == 'node_wave': 
                tmpnode = [float(s) for s in line.split()[1::]]
                if availgrism.lower() == 'g102':
                    nodelam = [x for x in tmpnode if x <= config_pars['transition_wave']]
                elif availgrism.lower() == 'g141':
                    nodelam = [x for x in tmpnode if x > config_pars['transition_wave']]
                else:
                    nodelam = tmpnode
                config_pars.setdefault('node_wave', []) 
                for l in nodelam:  
                    config_pars['node_wave'].append(l)
            elif name == 'mask_region1':
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region1', []) 
                config_pars['mask_region1'].append(masklam[0]) 
                config_pars['mask_region1'].append(masklam[1]) 
            elif name == 'mask_region2': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region2', []) 
                config_pars['mask_region2'].append(masklam[0]) 
                config_pars['mask_region2'].append(masklam[1])
            elif name == 'mask_region3': 
                masklam = [float(s) for s in line.split()[1::]]
                config_pars.setdefault('mask_region3', []) 
                config_pars['mask_region3'].append(masklam[0]) 
                config_pars['mask_region3'].append(masklam[1])  
            else:  
                val  = line.split()[1] 
                if is_number(val) == True:  val = float(val)
                if val == 'True': val = True
                if val == 'False': val = False 
                config_pars[line.split()[0]] = val 
    configfile.close()
    return config_pars    



# Creates .reg files for direct and dispersed images
# Written by Mason Huberty and adapted by KVN to work as part of the line finding software

def write_obj_region(parno, region_file_path, catalog, regfile_name, xoffset, yoffset, w, b_width, b_length):
    file = open(region_file_path + f'Par{parno}{regfile_name}', 'a')
    for i in range(len(catalog)):
        ra, dec = catalog['ra'][i], catalog['dec'][i]
        x, y = w.all_world2pix(ra, dec, 1)
        file.write(f'box({x-xoffset},{y-yoffset},{b_width},{b_length},0.0) # color=red text={{catalog["id"][i]}} font=times 10 bold italic textangle=30\n')
    file.close



# NOTE - The x & y offset values are specific for NIRISS. 
#        These will need to be updated for NIRCam data!
def create_regions(parno, args):

    spec_cat = glob.glob(args.speccat_file_path + 'Par*spec*.fits')
    hdul = fits.open(spec_cat[0])
    cat=hdul[1].data

    # create region file directory if needed
    os.makedirs(args.region_file_path, exist_ok=True)

    # Direct image region files
    f = open(args.region_file_path + f'Par{parno}regions_phot.reg','a')
    for i in range(len(cat)):
        f.write(f'WCS;circle({cat["ra"][i]},{cat["dec"][i]},0.5') # color=green text={'+str(cat['id'][i])+' z='+str(round(cat['redshift'][i],3))+'} font='times 10 bold italic' textangle=30\n')
    f.close()

    #This and subsequent are for the first order beams. Offsets taken from the config files
    # KVN :  Adding code to first check if the paths exist. Only create region file if filter/orientation is available 
    f115grism_R = glob.glob(args.drizzled_images_path + '*f115w*gr150r_drz_sci.fits')
    if len(f115grism_R) != 0:
        write_obj_region(parno, args.region_file_path, cat, 'F115r_grism.reg', 0.6548681566074263, 33.73739138173772, w = WCS(f115grism_R[0]),
                        b_width = 10.0, b_length = 93.54)
    
    f115grism_C = glob.glob(args.drizzled_images_path + '*f115w*gr150c_drz_sci.fits') # Added KVN 19-Aug-2024
    if len(f115grism_C) != 0:
        write_obj_region(parno, args.region_file_path, cat, 'F115c_grism.reg', 31.91107156101387, 1.3922939626209256, w = WCS(f115grism_C[0]),
                        b_width = 97.28751330105166, b_length = 10.0)

    f150grism_R = glob.glob(args.drizzled_images_path + '*f150w*gr150r_drz_sci.fits') # Added KVN 19-Aug-2024
    if len(f150grism_R) != 0:
        write_obj_region(parno, args.region_file_path, cat, 'F150r_grism.reg', 0.6548681566074263, 106.79254657227568, w = WCS(f150grism_R[0]),
                        b_width = 10.0, b_length = 93.54)
    
    f150grism_C = glob.glob(args.drizzled_images_path + '*f150w*gr150c_drz_sci.fits') # Added KVN 19-Aug-2024
    if len(f150grism_C) != 0:
        write_obj_region(parno, args.region_file_path, cat, 'F150c_grism.reg', 96.44444, 0.6548681566074263, w = WCS(f150grism_C[0]),
                        b_width = 93.54, b_length = 10.0)
    
    f200grism_R = glob.glob(args.drizzled_images_path + '*f200w*gr150r_drz_sci.fits') # Added KVN 19-Aug-2024
    if len(f200grism_R) != 0:    
        write_obj_region(parno, args.region_file_path, cat, 'F200r_grism.reg', 0.6548681566074263, 204.8370874255101, w = WCS(f200grism_R[0]),
                        b_width = 10.0, b_length = 131.78)        
    
    f200grism_C = glob.glob(args.drizzled_images_path + '*f200w*gr150c_drz_sci.fits') # Added KVN 19-Aug-2024
    if len(f200grism_C) != 0:
        write_obj_region(parno, args.region_file_path, cat, 'F200c_grism.reg', 200.9228, 0.6548681566074263, w = WCS(f200grism_C[0]),
                        b_width = 127.806, b_length = 10.0)        
    
