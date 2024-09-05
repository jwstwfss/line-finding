import astropy.io.fits as fits
from astropy.wcs import WCS
import numpy as np
import os
from passage_analysis import *
from distutils.sysconfig import *
from array import array
from glob import glob
import pandas as pd
from ast import literal_eval

# PASSAGE gui helpers
from guis_helpers import extract_image_extensions_key
from guis_helpers import display_images_in_DS9

# xpa resources:
# https://www.astro.louisville.edu/software/xmccd/archive/xmccd-4.1/xmccd-4.1e/docs/xpa/xpa.pdf
# http://ds9.si.edu/doc/ref/xpa.html#lock    


def showSpec2D_PASSAGE(parno, obid, path_to_spec2D=''):
    """Display spec2D cutouts in DS9.

    This can handle missed data in a given filter since there is a key
    to check against given the nature of how the spec2D files are setup.
    Frames will be made for missing data but will show nothing to indicate
    that the data is missing and that display layout is preserved.

    Modified by AA on Sep 2024 to allow flexible directory structures
    """

    spec2D_file = path_to_spec2D + f'Par{parno}_{obid:05d}.2D.fits'

    # for a given data file extract the data information that we want to be displayed in ds9
    spec2D_key_DS9 = extract_image_extensions_key(spec2D_file)

    # display in ds9 instance of a given "title"
    SPEC2D_TITLE_DS9 = "PASSAGE_spec2D"

    # now display the results.
    for key, _ in spec2D_key_DS9.items():

        # get key-values for frame id and extension
        frame_id = spec2D_key_DS9[key]["frame_id"]
        ext = spec2D_key_DS9[key]["ext"]

        # frame id
        command = f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id}"
        os.system(command)

        # if ext is None:
        #     command = f"xpaset -p {SPEC2D_TITLE_DS9} frame {frame_id} clear"
        #     os.system(command)

        # # display frame
        if ext is not None:
            command = f"xpaset -p {SPEC2D_TITLE_DS9} fits {spec2D_file}[{ext}]"
            os.system(command)
        else:
            command = f"xpaset -p {SPEC2D_TITLE_DS9} frame clear"
            os.system(command)

    # formating display into grid
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid mode manual")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} tile grid layout 5 3")

    # # config display properties
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame 1")
    # os.system("xpaset -p {SPEC2D_TITLE_DS9} lock frame physical")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock colorbar")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scale")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} lock scalelimits")
    os.system(f"xpaset -p {SPEC2D_TITLE_DS9} cmap grey")
    
    for fno in [1,2,3,4,5]:
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " +str(fno))
        # 39.5 and 32.5 are just 1/2 cutout length & 1/2 cutout width. Will need to be changed for NIRCam
        # The length is set to the length of the cutout -8 (for 4 pix on each end)
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 41 32.5 62 10# color=green} ")
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom to fit")

    for fno in [6,7,8,9,10]:
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " +str(fno))
        # 39.5 and 32.5 are just 1/2 cutout length & 1/2 cutout width + a small offset derived by eye. 
        # Will need to be changed for NIRCam
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 48 32.5 76 10# color=green} ")
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom to fit")

    for fno in [11,12,13,14,15]:
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} frame " +str(fno))
        # 70 and 33 are just 1/2 cutout length & 1/2 cutout width + a small offset derived by eye.
        # Will need to be changed for NIRCam
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9}"+" region command {box 71.5 33 108 10# color=green} ")
        os.system(f"xpaset -p {SPEC2D_TITLE_DS9} zoom to fit")


def find_file(directory, filename):
    """
    Scans a directory to find a file with the given name.

    Parameters:
    directory (str): The path of the directory to scan.
    filename (str): The name of the file to find.

    Returns:
    str: The path to the file if found, otherwise None.
    """
    for root, dirs, files in os.walk(directory):
        if filename in files:
            return os.path.join(root, filename)
    return None


def showDirect_PASSAGE(parno, path_to_drizzled_images='', path_to_region_files=''):
    '''
    Displays direct images for each of the filters
    Modified by AA on Sep 2024 to allow flexible directory structures
    '''

    ### KVN's quick fix to images having different names in different fields
    # specify the images to be displayed in DS9
    grism_file = glob(path_to_drizzled_images + '*gr150*_drz_sci.fits')
    grism_file_ext = os.path.basename(grism_file[0]).split('_')[1] # AA 2024/09/05: modified to include os.path.basename()


    # specify the images to be displayed in DS9
    images = {
        "f115w": [
            f"Par{parno}_f115w_drz_sci.fits",
            f"Par{parno}_{grism_file_ext}_f115w-gr150c_drz_sci.fits",
            f"Par{parno}_{grism_file_ext}_f115w-gr150r_drz_sci.fits",
        ],
        "f150w": [
            f"Par{parno}_f150w_drz_sci.fits",
            f"Par{parno}_{grism_file_ext}_f150w-gr150c_drz_sci.fits",
            f"Par{parno}_{grism_file_ext}_f150w-gr150r_drz_sci.fits",
        ],
        "f200w": [
            f"Par{parno}_f200w_drz_sci.fits",
            f"Par{parno}_{grism_file_ext}_f200w-gr150c_drz_sci.fits",
            f"Par{parno}_{grism_file_ext}_f200w-gr150r_drz_sci.fits",
        ],
    }

    print(f'\nDeb guis.py 146: parno={parno}, grism_file_ext={grism_file_ext}, images dict={images}') ##

    image_paths = {}
    for filter_name, filenames in images.items():
        paths = []
        for filename in filenames:
            full_path = find_file(path_to_drizzled_images, filename)
            if full_path:
                paths.append(full_path)
            else:
                print(f"Warning: File {filename} not found.")
                paths.append(None)
                
        image_paths[filter_name] = paths

    # specify the region files to show along with images, if any.
    region_files = {
        "f115w": [f"Par{parno}regions_phot.reg", f"Par{parno}F115c_grism.reg", f"Par{parno}F115r_grism.reg"],
        "f150w": [f"Par{parno}regions_phot.reg", f"Par{parno}F150c_grism.reg", f"Par{parno}F150r_grism.reg"],
        "f200w": [f"Par{parno}regions_phot.reg", f"Par{parno}F200c_grism.reg", f"Par{parno}F200r_grism.reg"],
    }

    # find the full paths of the region files
    region_paths = {}
    for filter_name, filenames in region_files.items():
        paths = []
        for filename in filenames:
            full_path = find_file(path_to_region_files, filename)
            if full_path:
                paths.append(full_path)
            else:
                print(f'Warning: File {filename} not found.')
                paths.append(None)
                
        region_paths[filter_name] = paths

    display_images_in_DS9(image_paths, region_files=region_paths)

    return


def panDirect_PASSAGE(ra, dec):
    # print(ra, dec)

    ds9_title = "PASSAGE_DIRECT"
    for fno in [1,4,7]:
        cmd = f"xpaset -p {ds9_title} frame " + str(fno)
        os.system(cmd)
        cmd = f"xpaset -p {ds9_title} pan to {ra[0]} {dec[0]} wcs degrees"
        os.system(cmd)

        # zoom to
        cmd = f"xpaset -p {ds9_title} zoom 2"
        os.system(cmd)


# written by KVN June 2024
# First pull info from region file
def getRegionFileInfo(parno, filt, path_to_region_files):
    '''
    Modified by AA on Sep 2024 to allow flexible directory structures
    '''

    region_file_name = path_to_region_files + f'Par{parno}{filt}_grism.reg'
    if os.path.exists(region_file_name):
            df = pd.read_csv(region_file_name, sep=' ', engine='python', header=None)

            # pull the id and x & y coord from the region file
            obj_ids = [i for i in df[3].apply(lambda z: int(z.split('=')[-1].strip('}{')))]
            obj_xs = df[0].apply(lambda z: literal_eval(z.lstrip('box'))[0])
            obj_ys = df[0].apply(lambda z: literal_eval(z.lstrip('box'))[1])
    else:
            print('No region file found for panning to object')
            obj_ids = None; obj_xs = None; obj_ys = None

    return obj_ids, obj_xs, obj_ys

# written by KVN June 2024
# Second, use info from region file to pan to each grism
def panDispersed_PASSAGE(objid, parno, path_to_drizzled_images, path_to_region_files):
    '''
    Modified by AA on Sep 2024 to allow flexible directory structures
    '''

    grism_file = glob(path_to_drizzled_images + '*gr150*_drz_sci.fits')
    grism_file_ext = grism_file[0].split('_')[1]

    ds9_title = "PASSAGE_DIRECT"
    # since tiles are already assigned, setup by tile rather than grism

    # AA added looping over filters and orientations to reduce blocks of code
    filters = ['115', '150', '200']
    orients = ['r', 'c']
    fno_dict = {'115c':2, '115r':3, '150c':5, '150r':6, '200c':8, '200r':9}

    for thisfilter in filters:
        for thisorient in orients:
            thisfile = path_to_drizzled_images + f'Par{parno}_{grism_file_ext}_f{thisfilter}w-gr150{thisorient}_drz_sci.fits'
            if os.path.exists(thisfile):
                fno = fno_dict[f'{thisfilter}{thisorient}']
                obj_ids, obj_xs, obj_ys = getRegionFileInfo(parno=parno, path_to_region_files=path_to_region_files, filt=f'{thisfilter}{thisorient}')

                # specify the frame
                cmd = f"xpaset -p {ds9_title} frame " + str(fno)
                os.system(cmd)

                # pan to the coords from the region file for obj with objid
                get_ind = [i for i in range(len(obj_ids)) if obj_ids[i] == objid]
                cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
                os.system(cmd)
                cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f " % (obj_xs[get_ind[0]], obj_ys[get_ind[0]])
                os.system(cmd)

# written by KVN June 2024
# DO NOT use this one. Most of below was testing. Faster/Better to get the coords from the region file
# Not the best way to do this... but the positions are in RA/DEC
# so I convert from RA/DEC to pixels then back to RA/DEC
def panDispersed_PASSAGE_doNOTuse(ra, dec, parno, path_to_wisp=""):

    grism_file = glob(path_to_data + '/Par'+str(parno)+'/DATA/*gr150*_drz_sci.fits')
    grism_file_ext = grism_file[0].split('_')[1]

    pix_per_um = 1 / (1e-4 * 46.934)  # GR150R 47.015 for GR150C
    tweak = 7

    ds9_title = "PASSAGE_DIRECT"    
    # since tiles are already assigned, setup by tile rather than grism
    for fno in [2,3]:
        if fno == 2 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150c_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150c_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position
            print('f115: ', xx,yy)

            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0# -6 + 0.135 * pix_per_um + tweak
            #updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)
        elif fno == 3 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150r_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f115w-gr150r_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position

            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset,  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)

    for fno in [5,6]:
        if fno == 5 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f150w-gr150c_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f150w-gr150c_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position
            print('f150: ', xx,yy)

            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0 #55 + 0.135 * pix_per_um + tweak
            #updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)
        else: # R orientation
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset,  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)


    for fno in [8,9]:
        if fno == 8 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150c_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150c_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position
            print('f200: ',xx,yy)
            
            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0#147 + 0.135 * pix_per_um + tweak
            #updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy, 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)

        elif fno == 9 and os.path.exists(path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150r_drz_sci.fits'):
            img_path = path_to_data + '/Par'+str(parno)+'/DATA/Par'+str(parno)+'_'+str(grism_file_ext)+ '_f200w-gr150r_drz_sci.fits'
            w = WCS(img_path)
            xx, yy = w.all_world2pix(ra, dec, 0) #this gives me the pixel coordinates of the object at (ra, dec) position

            grinfits = fits.open(img_path)
            hdr = grinfits[0].header
            if "PA_APER" in hdr:
                pa_aper = hdr["PA_APER"]
            else: pa_aper = 100

            pa_rad = -pa_aper / (180 / math.pi)
            cmd='xpaset -p PASSAGE_DIRECT frame ' + str(fno)
            os.system(cmd)

            obj_offset = 0#147 + 0.135 * pix_per_um + tweak
            updated_ra, updated_dec = w.all_pix2world(xx-obj_offset* math.sin(pa_rad),  yy-obj_offset* math.cos(pa_rad), 0)
            cmd = "xpaset -p PASSAGE_DIRECT pan to %f %f wcs degrees" % (updated_ra, updated_dec)
            os.system(cmd)
        

def show2dNEW(grism, parno, obid, zeroarr, user, trans, args, zran1=None, zran2=None):
    '''
    AA modified function header and body to include path variables as attributes of args
    '''
    # In version 1.0, will first look for wavelength-calibrated stamps in the G1??_DRIZZLE directories; failing this, will default to old stamps
    # zero and first order positions
    #    firstx = firstarr['x']
    #    firsty = firstarr['y']
    #    firstlen = firstarr['len']
    #    firstwid = firstarr['width']
    #    firstid = firstarr['objid']
    zerox = zeroarr["x"]
    zeroy = zeroarr["y"]
    zeroid = zeroarr["objid"]
    zerora = zeroarr["ra"]
    zerodec = zeroarr["dec"]
    zeromag = zeroarr["mag_auto"]

    dims = ()
    zrad = 10.0
    workingdir = os.getcwd()

    path2d = args.spec2D_path + f'Par{parno}_{obid:05d}.2D.fits'
    print(f'\nDeb guis.py 397: obid={obid:05d}, path2d={path2d}') ##

    pix_per_um = 1 / (1e-4 * 46.934)  # GR150R 47.015 for GR150C

    tweak = 7
    tweak_zp = 7
    if grism == "F200W":
        frameno = "3"
        maglimit = 26.0
        input_grism = glob(args.drizzled_images_path + f'Par{parno}_*f200w-gr150r_drz_sci.fits')[0]
        zp_offset = -210
        obj_offset = (147 + 0.238 * pix_per_um + tweak)  # First number is pixels to blue edge (poorly defined). Second num is wave distance to center in um.
    elif grism == "F150W":
        frameno = "2"
        maglimit = 24.5
        input_grism = glob(args.drizzled_images_path + f'Par{parno}_*f150w-gr150r_drz_sci.fits')[0]
        zp_offset = -219 + tweak_zp
        obj_offset = 55 + 0.171 * pix_per_um + tweak
    elif grism == "F115W":
        frameno = "1"
        maglimit = 23.5
        input_grism = glob(args.drizzled_images_path + f'Par{parno}_*f115w-gr150r_drz_sci.fits')[0]
        zp_offset = -225.5 + tweak_zp
        obj_offset = -6 + 0.135 * pix_per_um + tweak

    direct1 = (args.drizzled_images_path + f'Par{parno}_f115w_drz_sci.fits')
    direct2 = (args.drizzled_images_path + f'Par{parno}_f150w_drz_sci.fits')
    direct3 = (args.drizzled_images_path + f'Par{parno}_f200w_drz_sci.fits')

    if os.path.exists(direct1) == 1:
        input_direct = direct1
    elif os.path.exists(direct2) == 1:
        input_direct = direct2
    elif os.path.exists(direct3) == 1:
        input_direct = direct3

    tiltedgrism = 0
    extension_check = 0
    if os.path.exists(path2d) == 1:
        infits = fits.open(path2d)
        dal = len(infits)
        doug = list(range(1, dal))
        for tracker in doug:
            if grism == infits[tracker].header["EXTVER"]:
                extension_check = 1

    print(f'\nDeb guis.py 445: parno={parno}, extension_check={extension_check}, input_direct={input_direct}, input_grism={input_grism}') ##

    if extension_check:
        cutout_index = infits.index_of(("SCI", grism))
        grinfits = fits.open(input_grism)
        hdr = grinfits[0].header
        if "PA_APER" in hdr:
            pa_aper = hdr["PA_APER"]
        else:
            newfits = fits.open(input_direct)
            nhdr = newfits[0].header
            pa_aper = nhdr["PA_APER"]
            tiltedgrism = 1
        # print("PA_APER: ", pa_aper)
        wcs_in = WCS(header=grinfits[0].header)
        ra = array("d", zerora)
        dec = array("d", zerodec)
        x, y = wcs_in.all_world2pix(ra, dec, 0)
        x_obj = x[obid - 1]
        y_obj = y[obid - 1]

        # This section probably needs a grism orientation check G150R vs G150C
        if tiltedgrism:
            print(
                "Attention: Reference grism has been tilted from original orientation."
            )
            pa_rad = -pa_aper / (180 / math.pi)
            #  print(pa_rad,zp_offset*math.cos(pa_rad),zp_offset*math.sin(pa_rad))
            y_zp = y - zp_offset * math.cos(pa_rad)
            x_zp = x - zp_offset * math.sin(pa_rad)
            x_spec = x_obj - obj_offset * math.sin(pa_rad)
            y_spec = y_obj - obj_offset * math.cos(pa_rad)
            x_test = x - obj_offset * math.sin(pa_rad)
            y_test = y - obj_offset * math.cos(pa_rad)
        else:
            y_zp = y - zp_offset
            x_zp = x
            x_spec = x_obj
            y_spec = y_obj - obj_offset
            x_test = x
            y_test = y - obj_offset

        ### changing to read in 1st data extension ###
        # darr=infits[-1].data
        hdr = infits[cutout_index].header
        darr = infits[cutout_index].data
        rms = np.std(darr)
        dims = darr.shape

        y_c = 0.5 * (dims[0] + 1)
        x_c = 0.5 * (dims[1] + 1)
        x_delta = y_c - x_spec
        y_delta = x_c - y_spec
        x_zp_new = dims[1] - y_zp - y_delta
        y_zp_new = x_delta + x_zp
        x_test_new = dims[1] - y_test - y_delta
        y_test_new = x_delta + x_test
        #     print(x_test[obid-1],y_test[obid-1],x_test_new[obid-1],y_test_new[obid-1],dims[1],dims[0])
        xmax = dims[1]
        ymax = dims[0]
        infits.close()
        grinfits.close()

        if zran1 == None:
            zran1 = -1 * rms
        if zran2 == None:
            zran2 = 5 * rms

    elif (os.path.exists(path2d) == 0) or (extension_check == 0):
        print("%s stamp not found." % (grism))
        return False

    # COMMENTING OUT ZERO ORDER DISPLAY -START

    ### USING THE DRIZZLE TRANSFORMATIONS TO GET ZEROTH ORDERS ###
    #    _cx = np.array([xcoo for xcoo in zerox]) - hdr['BB0X'] - 1
    #    _cy = np.array([ycoo for ycoo in zeroy]) - hdr['BB0Y'] - 1
    #    cx = hdr['D001OUXC'] + (hdr['DRZ00'] + hdr['DRZ01']*(_cx-hdr['D001INXC']) + hdr['DRZ02']*(_cy-hdr['D001INYC']))
    #    cy = hdr['D001OUYC'] + (hdr['DRZ10'] + hdr['DRZ11']*(_cx-hdr['D001INXC']) + hdr['DRZ12']*(_cy-hdr['D001INYC']))
    # convert to (Angs,arcsec) coords
    #    cx = (cx - hdr['CRPIX1'])*hdr['CDELT1'] + hdr['CRVAL1']
    #    cy = (cy - hdr['CRPIX2'])*hdr['CDELT2'] + hdr['CRVAL2']
    #    rad = 5 * hdr['CD1_1']

    outcoo = (args.spec1D_path + "temp_zero_coords_%s_" % user + str(frameno) + ".reg")

    if os.path.exists(outcoo):
        os.unlink(outcoo)
    f = open(outcoo, "w")
    # f.write('wcs;\n')
    xmax1 = 3100
    ymax1 = 3100

    test = np.where(
        (x_zp_new > 0)
        & (x_zp_new < xmax)
        & (y_zp_new > 0)
        & (y_zp_new < ymax)
        & (zeromag < maglimit)
    )
    #  test=np.where(((x_zp > 0) & (x_zp < xmax1) & (y_zp > 0) & (y_zp < ymax1) & (zeromag < maglimit)))
    #    test=np.where(((x > 0) & (x < xmax1) & (y > 0) & (y < ymax1) & (zeromag < maglimit)))
    cross_size = 50
    rad = 5
    f.write(
        "point(%.2f,%.2f) # point=cross %i color=yellow text={%.1f}\n"
        % (x_test_new[obid - 1], y_test_new[obid - 1], cross_size, zeromag[obid - 1])
    )
    if test[0] != []:
        for j in range(0, len(test[0])):
            inputx = x_zp_new[test[0][j]]
            #    inputx=x_zp[test[0][j]]
            #    inputx=x[test[0][j]]
            inputy = y_zp_new[test[0][j]]
            #     inputy=y_zp[test[0][j]]
            #     inputy=y[test[0][j]]
            inputID = zeroid[test[0][j]]
            inputmag = zeromag[test[0][j]]
            #      testx=x_test[test[0][j]]
            #      testy=y_test[test[0][j]]
            print(inputx, inputy, rad, inputID)
            f.write(
                "circle(%.2f,%.2f,%.1f) # color=red text={%.1f;%s}\n"
                % (inputx, inputy, rad, inputmag, inputID)
            )
    #      f.write('circle(%.2f,%.2f,%.1f) # color=blue text={%.1f;%s}\n' % (testx,testy,rad,inputmag,inputID))
    f.close()

    # COMMENTING OUT ZERO ORDER DISPLAY -END

    return


def showDirectNEW(obid, parno, g102zeroarr, path_to_drizzled_images, load_images=False):
    """
    Removed lineno, which was only used to check whether the images
    should be reloaded.
    AA modified function header and body to include path variables as attributes of args
    """
    obid = int(obid)
    workingdir = os.getcwd()
    dirpts = workingdir.split("/")[1:-1]

    path110 = path_to_drizzled_images + f"Par{parno}_f115w_drz_sci.fits"
    path140 = path_to_drizzled_images + f"Par{parno}_f150w_drz_sci.fits"
    path160 = path_to_drizzled_images + f"Par{parno}_f200w_drz_sci.fits"

    inputx = g102zeroarr["x"]
    inputy = g102zeroarr["y"]

    if (os.path.exists(path110) == 0
        and os.path.exists(path140) == 0
        and os.path.exists(path160) == 0
    ):
        print("No Direct Images Found.")
        return 0

    hexcoo = [str(inputx[obid - 1]), str(inputy[obid - 1])]

    # Added getting PA angle of each image.

    '''
    # load the direct images
    if load_image:
        print("IMAGE IS BEING RELOADED AND ROTATION APPLIED!")
        if os.path.exists(path110) == 1:
            hdu110 = fits.open(path110)
            hdr = hdu110[0].header
            pa_aper = hdr["PA_APER"]
            cmd = "xpaset -p ds9 frame 4"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path110
            os.system(cmd)
            ### using F110_drz.reg with F110W_drz.fits
            #    cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F110_drz.reg'
            #    os.system(cmd)
            cmd = "xpaset -p ds9 rotate to " + str(-pa_aper + 90)
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to " + hexcoo[0] + " " + hexcoo[1]  # +' fk5'
            os.system(cmd)
        if os.path.exists(path140) == 1:
            hdu140 = fits.open(path140)
            hdr = hdu140[0].header
            pa_aper = hdr["PA_APER"]
            cmd = "xpaset -p ds9 frame 5"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path140
            os.system(cmd)
            #      cmd='xpaset -p ds9 regions file '+par_root_dir+ 'DATA/DIRECT_GRISM/F140_drz.reg'
            #     os.system(cmd)
            cmd = "xpaset -p ds9 rotate to " + str(-pa_aper + 90)
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to " + hexcoo[0] + " " + hexcoo[1]  # +' fk5'
            os.system(cmd)
        if os.path.exists(path160) == 1:  # was elif
            hdu160 = fits.open(path160)
            hdr = hdu160[0].header
            pa_aper = hdr["PA_APER"]
            cmd = "xpaset -p ds9 frame 6"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path160
            print(cmd)
            os.system(cmd)
            #      cmd='xpaset -p ds9 regions file '+par_root_dir+'DATA/DIRECT_GRISM/F160_drz.reg'
            #     os.system(cmd)
            cmd = "xpaset -p ds9 rotate to " + str(-pa_aper + 90)
            print(cmd)
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to " + hexcoo[0] + " " + hexcoo[1]  # +' fk5'
            print(cmd)
            os.system(cmd)
    '''
    # pan to the coordinates of this object
    if os.path.exists(path110):
        panDirect(hexcoo[0], hexcoo[1], grism="F115W")
    panDirect(hexcoo[0], hexcoo[1])
    if os.path.exists(path140):
        panDirect(hexcoo[0], hexcoo[1], grism="F150W")
    panDirect(hexcoo[0], hexcoo[1])
    if os.path.exists(path160):
        panDirect(hexcoo[0], hexcoo[1], grism="F200W")
    panDirect(hexcoo[0], hexcoo[1])


def showDispersed(obid, parno, path_to_dispersed, load_image=False):  # MB
    """
    Removed lineno, which was only used to check whether the images
    should be reloaded.
    AA modified function header and body to include path variables as attributes of args
    """
    obid = int(obid)
    workingdir = os.getcwd()
    dirpts = workingdir.split("/")[1:-1]

    ### Using G102.fits instead of G102_drz.fits ###
    path102 = (path_to_dispersed + "Par" + str(parno) + "_" + str(obid).zfill(5) + ".2D.fits")
    path141 = (path_to_dispersed + "Par" + str(parno) + "_" + str(obid).zfill(5) + ".2D.fits")
    path102_0reg = os.path.join(path_to_dispersed, "G102_0th.reg")
    path102_1reg = os.path.join(path_to_dispersed, "G102_1st.reg")
    path141_0reg = os.path.join(path_to_dispersed, "G141_0th.reg")
    path141_1reg = os.path.join(path_to_dispersed, "G141_1st.reg")
    if os.path.exists(path102) == 0 and os.path.exists(path141) == 0:
        print("Looking for objects here: ", path102)
        print("No Grism Images Found.")
        return 0
    # get center of 1st order
    ### Using same syntax as getzeroorders for consistency ###
    print(path102_1reg)
    if os.path.exists(path102_1reg) == 1:
        reg102 = open(path102_1reg, "r")
        x102, y102 = -1, -1
        for line in reg102:
            # using same syntax as getzeroorders
            linesplit = line.split()
            textid = int(re.search("\d+", linesplit[-1]).group(0))
            if textid == obid:
                x102 = float(linesplit[1].split(",")[0])
                y102 = float(linesplit[2].split(",")[0])
        reg102.close()
    if os.path.exists(path141_1reg) == 1:
        reg141 = open(path141_1reg, "r")
        x141, y141 = -1, -1
        for line in reg141:
            linesplit = line.split()
            textid = int(re.search("\d+", linesplit[-1]).group(0))
            if textid == obid:
                x141 = float(linesplit[1].split(",")[0])
                y141 = float(linesplit[2].split(",")[0])
        reg141.close()

    '''
    if load_image:
        if os.path.exists(path102) == 1:
            cmd = "xpaset -p ds9 frame 5"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path102
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path102_0reg
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path102_1reg
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to %f %f image" % (x102, y102)
            os.system(cmd)
        if os.path.exists(path141) == 1:
            cmd = "xpaset -p ds9 frame 6"
            os.system(cmd)
            cmd = "xpaset -p ds9 file " + path141
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path141_0reg
            os.system(cmd)
            cmd = "xpaset -p ds9 regions file " + path141_1reg
            os.system(cmd)
            cmd = "xpaset -p ds9 pan to %f %f image" % (x141, y141)
            os.system(cmd)
    '''
    # pan to the coordinates of this object
    if os.path.exists(path102):
        panDispersed(x102, y102, grism="G102")
    if os.path.exists(path141):
        panDispersed(x141, y141)


### AH: 1/26/28: I do not think we are using this?
def createAltGrismRegion(grism):
    workingdir = os.getcwd()
    par = os.path.dirname(workingdir)
    cat110 = os.path.join(par, "DATA/DIRECT_GRISM/fin_F110.cat")
    cat140 = os.path.join(par, "DATA/DIRECT_GRISM/fin_F140.cat")
    cat160 = os.path.join(par, "DATA/DIRECT_GRISM/fin_F160.cat")

    if grism == "G102":
        if os.path.exists(cat110) == 1:
            cat = np.genfromtxt(
                cat110,
                dtype=[("num", int), ("a_img", float), ("mag", float)],
                usecols=(1, 4, 12),
            )
        else:
            print(cat110)
            return 0
    if grism == "G141":
        if os.path.exists(cat140) == 1:
            cat = np.genfromtxt(
                cat140,
                dtype=[("num", int), ("a_img", float), ("mag", float)],
                usecols=(1, 4, 12),
            )
        elif os.path.exists(cat160) == 1:
            cat = np.genfromtxt(
                cat160,
                dtype=[("num", int), ("a_img", float), ("mag", float)],
                usecols=(1, 4, 12),
            )
        else:
            print("nope2")
            return 0

    f = open(os.path.join(workingdir, grism + "_temp.reg"), "w")
    f.write(
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
    )
    f.write("image\n")
    cenx = np.zeros(cat["num"].shape[0], dtype=float)
    ceny = np.zeros(cat["num"].shape[0], dtype=float)
    for i in range(cat["num"].shape[0]):
        objid = cat["num"][i]
        gfile = os.path.join(
            par, "%s_DRIZZLE/aXeWFC3_%s_mef_ID%i.fits" % (grism, grism, objid)
        )
        if os.path.isfile(gfile):
            hdr = fits.getheader(gfile, 1)
            # coords of bounding box
            boxx = np.array([hdr["bb0x"], hdr["bb1x"]])
            boxy = np.array([hdr["bb0y"], hdr["bb1y"]])
            cenx[i] = boxx[0] + (boxx[1] - boxx[0]) / 2.0
            ceny[i] = boxy[0] + (boxy[1] - boxy[0]) / 2.0
            slitwidth = hdr["slitwidt"]
            sx = 184
            sy = 8
            #            sy = slitwidth * cat['a_img'][i]
            f.write(
                "box(%f,%f,%i,%i,0) # text={%i}\n" % (cenx[i], ceny[i], sx, sy, objid)
            )
    f.close()
    return cenx, ceny


def panDirect(ra, dec, grism="F200W"):
    # Pan to coords in frame
    if grism == "F115W":
        fno = "4"
    elif grism == "F150W":
        fno = "5"
    elif grism == "F200W":
        fno = "6"
    cmd = "xpaset -p ds9 frame " + fno
    os.system(cmd)
    cmd = "xpaset -p ds9 pan to " + ra + " " + dec  # +' fk5'
    os.system(cmd)


def panDispersed(xx, yy, grism="G141"):  # MB
    # Pan to coords in frame
    if grism == "G141":
        fno = "6"
    else:
        fno = "5"
    cmd = "xpaset -p ds9 frame " + fno
    os.system(cmd)
    cmd = "xpaset -p ds9 pan to %f %f image" % (xx, yy)
    os.system(cmd)


def reloadReg(region_file_path):
    '''
    Modified by AA on Sep 2024 to allow flexible directory structures
    '''
    # reload direct image region files
    # AA modified below to reduce blocks of code

    region_file1 = region_file_path + f'F110_drz.reg'
    region_file2 = region_file_path + f'F160_drz.reg'
    region_file3 = region_file_path + f'F140_drz.reg'

    if os.path.exists(region_file1):
        cmd = "xpaset -p ds9 frame 3"
        os.system(cmd)
        cmd = (f"xpaset -p ds9 regions file {thisregion_file}")
        os.system(cmd)

    if os.path.exists(region_file2): thisregion_file = region_file2
    elif os.path.exists(region_file3): thisregion_file = region_file3
    else: thisregion_file = None

    if thisregion_file is not None:
        cmd = "xpaset -p ds9 frame 4"
        os.system(cmd)
        cmd = (f"xpaset -p ds9 regions file {thisregion_file}")
        os.system(cmd)
