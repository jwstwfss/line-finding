'''
    Filename: mainPASSAGE.py
    Notes: Interactive emission line finding routine for PASSAGE fields
    Author : Kalina Nedkova
    Modified: Ayan Acharyya
    Last modified: 3 September 2024
    Examples: run mainPASSAGE.py --user ayan_gdrive
             run mainPASSAGE.py --verbose
'''

try:
    import passage_analysis as passage
    from passage_analysis.utilities import *
except:
    print('passage_analysis could not be imported due to path issues, will remedy this soon.')

import os
import re
import glob
import sys
import argparse
import subprocess

import warnings
warnings.filterwarnings('ignore')

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    By AA; Sep 2024
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')
    parser.add_argument('--user', metavar='user', type=str, action='store', default='knedkova', help='Which user template to follow for directory structures?')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=False, help='Maximise prints to screen? Default is no.')
    parser.add_argument('--clobber_region', dest='clobber_region', action='store_true', default=False, help='Re-make *.reg files? Default is no.')
    parser.add_argument('--clobber_1D', dest='clobber_1D', action='store_true', default=False, help='Make *1D.dat outputs? Default is no.')
    parser.add_argument('--clobber_RC', dest='clobber_RC', action='store_true', default=False, help='Make *_R.fits and *_C.fits outputs? Default is no.')
    parser.add_argument('--clobber_linelist', dest='clobber_linelist', action='store_true', default=False, help='Re-make linelist file? Default is no.')
    args = parser.parse_args()

    return args

# -------------------------------------------------------------------------------------------------------
def get_user_directory_structure(args):
    '''
    Determines paths and directory structures paths based on the given user template
    Simply add your paths and directory structure as an additional if block if it does not exist already
    Returns the object args, with paths and path placeholders as attributes
    By AA; Sep 2024
    '''
    if args.user == 'knedkova':
        args.args.code_dir = '/Users/knedkova/Work/2024PASSAGE/passage_analysis/'
        args.args.data_dir = '/Users/knedkova/Work/2024PASSAGE/output/'
        args.output_dir = '/Users/knedkova/Work/2024PASSAGE/data/'

        # in the following paths, "FIELD" is a placeholder, later to be replaced by the relevant field name
        args.speccat_file_path = args.data_dir + 'ParFIELD/DATA/DIRECT_GRISM/' # this is where it will look for *speccat.fits
        args.photcat_file_path = args.data_dir + 'ParFIELD/DATA/DIRECT_GRISM/' # this is where it will look for *photcat.fits
        args.region_file_path = args.data_dir + 'ParFIELD/DATA/' # this is where it will look for any *.reg file
        args.drizzled_images_path = args.data_dir + 'ParFIELD/DATA/' # this is where it will look for any *drz*.fits file
        args.spectra_path = args.data_dir + 'ParFIELD/Spectra/' # this is where it will look for Par*G*.dat spectra files
        args.spec1D_path = args.data_dir + 'ParFIELD/spec1D/' # this is where it will look for *.1D.fits files
        args.spec2D_path = args.data_dir + 'ParFIELD/spec2D/' # this is where it will look for *.2D.fits files
        args.beam_files_path = args.data_dir + 'ParFIELD/beams/' # this is where it will look for *beam.fits files

        args.linelist_path = args.output_dir + 'linelist/' # this is where it will store the Par*linelist.dat files
        args.stored_fits_path = args.output_dir + 'ParFIELD/' # this is where it will make Par*_output*/ directory, and fitdata/ and figs/ within it, to store the outputs

        args.is_fieldname_padded = False # set to True if files/folders are named like Par005 instead of Par5 and so on

    elif 'ayan' in args.user:
        if args.user == 'ayan_local':
            args.root_dir = '/Users/acharyya/Work/astro/passage/'
        elif args.user == 'ayan_hd':
            args.root_dir = '/Volumes/Elements/acharyya_backup/Work/astro/passage/'
        elif args.user == 'ayan_gdrive':
            args.root_dir = '/Users/acharyya/Library/CloudStorage/GoogleDrive-ayan.acharyya@inaf.it/My Drive/passage/'

        args.code_dir = '/Users/acharyya/Work/astro/passage/line-finding/passage_analysis/'
        args.data_dir = args.root_dir + 'passage_data/'
        args.output_dir = args.root_dir + 'passage_output/'

        # in the following paths, "FIELD" is a placeholder, later to be replaced by the relevant field name
        args.speccat_file_path = args.data_dir + 'ParFIELD/Products/' # this is where it will look for *speccat.fits
        args.photcat_file_path = args.data_dir + 'ParFIELD/Products/' # this is where it will look for *photcat.fits
        args.region_file_path = args.data_dir + 'ParFIELD/Regions/' # this is where it will look for any *.reg file
        args.drizzled_images_path = args.data_dir + 'ParFIELD/Products/' # this is where it will look for any *drz*.fits file
        args.spectra_path = args.data_dir + 'ParFIELD/Products/Spectra/' # this is where it will look for Par*G*.dat spectra files
        args.spec1D_path = args.data_dir + 'ParFIELD/Products/spec1D/' # this is where it will look for *.1D.fits files
        args.spec2D_path = args.data_dir + 'ParFIELD/Products/spec2D/' # this is where it will look for *.2D.fits files
        args.beam_files_path = args.data_dir + 'ParFIELD/Extractions/' # this is where it will look for *beam.fits files

        args.linelist_path = args.output_dir + 'linelist/' # this is where it will store the Par*linelist.dat files
        args.stored_fits_path = args.output_dir + 'ParFIELD/' # this is where it will make Par*_output*/ directory, and fitdata/ and figs/ within it, to store the outputs

        args.is_fieldname_padded = True # set to True if files/folders are named like Par005 instead of Par5 and so on

    return args

# -------------------------------------------------------------------------------------------------------
def substitute_fieldname_in_paths(args, field, placeholder='FIELD'):
    '''
    Converts all pseudo-paths available in the attributes of args to proper paths by substituting field names appropriately
    Returns args
    By AA; Sep 2024
    '''
    all_attributes = list(args.__dict__.keys())
    pseudo_paths = [item for item in all_attributes if placeholder in str(args.__dict__[item])]

    for this_path in pseudo_paths:
        args.__dict__[this_path] = args.__dict__[this_path].replace(placeholder, field)

    return args

# ------------------------------------------------------------------------------------------------------
def close(procs):
    '''
    Tiny function to close any process opened by subprocess (in this case, the ds9 windows) given a list of process
    By AA; Sep 2024
    '''
    for proc in procs: proc.terminate()

# -------------------------------------------------------------------------------------------------------
if __name__ == "__main__":
    # -----------------------------Take field name from user input--------------------------------------
    parno = input('\033[94m' + "Enter the number of the parallel field you want to analyze.\n> " + '\033[0m')
    # pull out only the number, in case the user entered e.g. "Par1"
    while True:
        try:
            parno = int(re.findall(r'\d+', str(parno))[0])
        except:
            parno = input(
                '\033[94m' + "A parallel field number is required. Enter the number of the parallel field you want to analyze.\n> " + '\033[0m')
            continue
        else:
            break

    # -----get user name and other user args------
    args = parse_args()

    # -----get associated directory structure and appropriate parno------
    args = get_user_directory_structure(args)
    if args.is_fieldname_padded: parno = f'{parno:03}' # pad field name if needed
    args.parno = parno
    args = substitute_fieldname_in_paths(args, parno) # convert pseudo paths to proper paths

    # -----import passage_analysis if not successful earlier------
    if 'passage' not in locals():
        sys.path.insert(1, args.code_dir)
        import passage_analysis as passage
        from passage_analysis.utilities import *
        print('passage_analysis now successfully imported.')

    # ---------check if region files exist. If not, run code to create necessary region files---------
    regionfiles = glob.glob(args.region_file_path + '*.reg')
    if len(regionfiles) == 0:
        print('\033[94m' + "No region files found, creating those for you now."  + '\033[0m')
        create_regions(parno, args)
    else:
        print(f'\nFound some region files in {args.region_file_path}, so proceeding to the next step. If you want to re-make the region files please rerun mainPASSAGE.py with --clobber_region option.')

    # ---------check if headers for 2D spectra need to be updated, and update if necessary-----------
    if fits_headers_need_updating(args.spec2D_path):
        update_fits_headers(args.spec2D_path)

    # ---------check if .dat spectra files exist. If not, make them------------------
    spec_files = sorted(glob.glob(args.spectra_path + '*1D.dat'))
    if len(spec_files) == 0 or args.clobber_1D:
        convert_1Dspectra_fits_to_dat(args)
    else:
        print(f'\nFound 1D.dat files in {args.spectra_path}, so proceeding to the next step. If you want make *1D.dat files please rerun mainPASSAGE.py with --clobber_1D option.')

    # ---------check if R and C spectra files exist. If not, make them------------------
    C_files = sorted(glob.glob(args.spectra_path + '*C.dat'))
    R_files = sorted(glob.glob(args.spectra_path + '*R.dat'))
    if (len(C_files) == 0 and len(R_files) == 0) or args.clobber_RC:
        make_1D_spectra_per_orientation(args)
    else:
        print(f'\nFound R & C files in {args.spectra_path}, so proceeding to the next step. If you want to make *_R/C.fits files please rerun mainPASSAGE.py with --clobber_RC option.')

    # -------check if line list exists. If not, run code to create the linelist------------
    linelist_file = args.linelist_path + 'Par'+str(parno)+'lines.dat'
    if not os.path.exists(linelist_file):
        print('\033[94m' + "No line list file found, creating the line list for you now." + '\033[0m')
        passage.loop_field_cwt(parno, args)
    else:
        print(f'\nFound linelist file in {args.linelist_path}, so proceeding to the next step. If you want to re-make the linelist file please rerun mainPASSAGE.py with --clobber_linelist option.')

    # -----------------------------Open two ds9 windows--------------------------------------
    #os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_DIRECT &')
    #os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_spec2D &')
    ds9_direct = subprocess.Popen(['/Applications/SAOImageDS9.app/Contents/MacOS/ds9', '-title', 'PASSAGE_DIRECT'])
    ds9_2d = subprocess.Popen(['/Applications/SAOImageDS9.app/Contents/MacOS/ds9', '-title', 'PASSAGE_spec2D'])
    ds9 = [ds9_direct, ds9_2d] # accumulating process IDs in a list, so that we can close ALL if needed (just do 'close(ds9)' on terminal)

    # -----------run the measure_z_interactive codes-------------------
    args.show_dispersed = True
    args.print_colors = True
    passage.measure_z_interactive(parno, args)


