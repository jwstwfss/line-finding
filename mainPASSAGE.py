'''
    Filename: mainPASSAGE.py
    Notes: Interactive emission line finding routine for PASSAGE fields
    Author : Kalina Nedkova
    Modified: Ayan Acharyya
    Last modified: 3 September 2024
    Examples: run mainPASSAGE.py --user ayan_hd
             run mainPASSAGE.py --verbose
'''

try:
    import passage_analysis as passage
    from passage_analysis import utilities
except:
    print('passage_analysis could not be imported due to path issues, will remedy this soon.')

import os
import re
import glob
import sys
import argparse

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')
    parser.add_argument('--user', metavar='user', type=str, action='store', default='knedkova', help='Which user template to follow for directory structures?')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=False, help='Maximise prints to screen? Default is no.')
    args = parser.parse_args()

    return args

# -------------------------------------------------------------------------------------------------------
def get_paths_by_user(args):
    '''
    Determines paths and directory structures paths based on the given user template
    Simply add your paths and directory structure as an additional if block if it does not exist already
    Returns the object args, with paths and path placeholders as attributes
    By Ayan Acharyya
    '''
    if args.user == 'knedkova':
        args.args.code_dir = '/Users/knedkova/Work/2024PASSAGE/passage_analysis/'
        args.args.data_dir = '/Users/knedkova/Work/2024PASSAGE/output/'
        args.output_dir = '/Users/knedkova/Work/2024PASSAGE/data/'

        # in the following paths, "FIELD" is a placeholder, later to be replaced by the relevant field name
        args.speccat_file_path = args.data_dir + 'ParFIELD/DATA/DIRECT_GRISM/'
        args.photcat_file_path = args.data_dir + 'ParFIELD/DATA/DIRECT_GRISM/'
        args.region_file_path = args.data_dir + 'ParFIELD/DATA/'
        args.drizzled_images_path = args.data_dir + 'ParFIELD/DATA/'
        args.spec1D_path = args.data_dir + 'ParFIELD/Spectra/'
        args.spec2D_path = args.data_dir + 'ParFIELD/spec2D/'
        args.linelist_path = args.output_dir + 'linelist/'
        args.stored_fits_path = args.output_dir + 'ParFIELD/'

        args.is_fieldname_padded = False # set to True if files/folders are named like Par005 instead of Par5 and so on
        args.spectra_available_for_individual_filters = True # set to True if spectra files are named as *G115_1D.dat

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

        args.speccat_file_path = args.data_dir + 'ParFIELD/Products/'
        args.photcat_file_path = args.data_dir + 'ParFIELD/Products/'
        args.region_file_path = args.data_dir + 'ParFIELD/Regions/'
        args.drizzled_images_path = args.data_dir + 'ParFIELD/Products/'
        args.spec1D_path = args.data_dir + 'ParFIELD/Products/spec1D/'
        args.spec2D_path = args.data_dir + 'ParFIELD/Products/spec2D/'
        args.linelist_path = args.output_dir + 'linelist/'
        args.stored_fits_path = args.output_dir + 'ParFIELD/'

        args.is_fieldname_padded = True
        args.spectra_available_for_individual_filters = False

    return args

# -------------------------------------------------------------------------------------------------------
def substitute_fieldname_in_paths(args, field, placeholder='FIELD'):
    '''
    Converts all pseudo-paths available in the attributes of args to proper paths by substituting field names appropriately
    Returns args
    By Ayan Acharyya
    '''
    all_attributes = list(args.__dict__.keys())
    pseudo_paths = [item for item in all_attributes if placeholder in str(args.__dict__[item])]

    for this_path in pseudo_paths:
        args.__dict__[this_path] = args.__dict__[this_path].replace(placeholder, field)

    return args

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

    # -----get user name and associated directory structure------
    args = parse_args()
    args = get_paths_by_user(args)

    # -----pad field name if needed------
    if args.is_fieldname_padded: parno = f'{parno:03}'

    # -----convert pseudo paths to proper paths------
    args = substitute_fieldname_in_paths(args, str(parno))

    # -----import passage_analysis if not successful earlier------
    if 'passage' not in locals():
        sys.path.insert(1, args.code_dir)
        import passage_analysis as passage
        from passage_analysis import utilities
        print('passage_analysis now successfully imported.')

    # --------move to the directory where you want your outputs-----------
    #os.chdir(args.output_dir)

    # ---------check if region files exist. If not, run code to create necessary region files---------
    regionfiles = glob.glob(args.region_file_path + '*.reg')
    if len(regionfiles) == 0:
        print('\033[94m' + "No region files found, creating those for you now."  + '\033[0m')
        utilities.create_regions(parno, args)

    # -------check if line list exists. If not, run code to create the linelist------------
    linelist_file = args.linelist_path + 'Par'+str(parno)+'lines.dat'
    if not os.path.exists(linelist_file):
        print('\033[94m' + "No line list file found, creating the line list for you now." + '\033[0m')
        passage.loop_field_cwt(parno, args)
    
    # -----------------------------Open two ds9 windows--------------------------------------
    os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_DIRECT &')
    os.system('/Applications/SAOImageDS9.app/Contents/MacOS/ds9 -title PASSAGE_spec2D &')

    # -----------run the measure_z_interactive codes-------------------
    args.show_dispersed = True
    args.print_colors = True
    passage.measure_z_interactive(parno, args)


