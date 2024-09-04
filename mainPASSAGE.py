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
from datetime import datetime, timedelta
from astropy.io import fits
from astropy.table import Table
import numpy as np
from grizli import multifit

import warnings
warnings.filterwarnings('ignore')

# --------------------------------------------------------------------------------------------------------------------
def parse_args():
    '''
    Parse command line arguments. Returns args object
    '''

    parser = argparse.ArgumentParser(description='Produces emission line maps for JWST-PASSAGE data.')
    parser.add_argument('--user', metavar='user', type=str, action='store', default='knedkova', help='Which user template to follow for directory structures?')
    parser.add_argument('--verbose', dest='verbose', action='store_true', default=False, help='Maximise prints to screen? Default is no.')
    parser.add_argument('--clobber_1D', dest='clobber_1D', action='store_true', default=False, help='Over-write *1D.dat outputs? Default is no.')
    parser.add_argument('--clobber_RC', dest='clobber_RC', action='store_true', default=False, help='Over-write *_R.fits and *_C.fits outputs? Default is no.')
    args = parser.parse_args()

    return args

# -------------------------------------------------------------------------------------------------------
def get_user_directory_structure(args):
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
        args.spectra_available_for_individual_filters = True # set to True if spectra files are named as *G115_1D.dat and False if they are *.1D.fits

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
        args.spectra_available_for_individual_filters = True # set to True if spectra files are named as *G115_1D.dat and False if they are *.1D.fits

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
def make_table(tab):
    '''
    Modifies and returns input table, by changing normalisation of flux, etc
    Adapted by Ayan Acharyya from passage_analysis/passage_convert_data.ipynb
    '''
    tab.rename_columns(['err'], ['error'])

    for column in ['flux', 'error', 'contam']: tab[column] /= tab['flat']
    tab['zeroth'] = 0
    tab = tab.filled(0.0)  # Replace nans with zeros

    # Spectra dispersed beyond the chip have zero fluxes that must be replaced to prevent crashes in fitting.
    tab['flux'][np.where(tab['flux'] == 0.0)] = np.median(tab['flux'][np.where(tab['flux'] != 0.0)])
    tab['error'][np.where(tab['error'] == 0.0)] = np.median(tab['error'][np.where(tab['error'] != 0.0)])
    tab = tab['wave', 'flux', 'error', 'contam', 'zeroth']

    return tab

# -------------------------------------------------------------------------------------------------------
def convert_1Dspectra_fits_to_dat(args):
    '''
    Converts all existing *.1D.fits to .dat files, for each available filter, for a given field
    Adapted by Ayan Acharyya from passage_analysis/passage_convert_data.ipynb
    '''
    start_time = datetime.now()
    print(f'Starting convert_1Dspectra_fits_to_dat..')

    spec1d_files = sorted(glob.glob(args.spec1D_path + '*1D.fits'))
    os.makedirs(args.spectra_path, exist_ok=True)

    # ------looping over all 1D spectra files--------------
    for index, spec1d_filename in enumerate(spec1d_files):
        data = fits.open(spec1d_filename)
        print(f'\nDoing {index + 1} of {len(spec1d_files)} files..')

        for ext in range(1, len(data)):
            filter = data[ext].header['EXTNAME'].strip()

            if filter in ['F115W', 'F150W', 'F200W']:
                outfilename = args.spectra_path + os.path.basename(spec1d_filename).replace('1D.fits', f'G{filter[1:-1]}_1D.dat')

                if not os.path.exists(outfilename) or args.clobber_1D:
                    tab = Table(data[filter].data)
                    tab = make_table(tab)

                    if filter == 'F200W':
                        try:
                            wave_lim = np.max(Table(data[ext - 1].data)['wave'])
                            tab = tab[tab['wave'] > wave_lim]
                        except:
                            pass

                    # Write out the updated files.
                    tab.write(outfilename, format='ascii.fixed_width_two_line', overwrite=True)
                    print(f'Written {outfilename}')
                else:
                    print(f'It appears the .dat files were already created for this object/filter. Skipping this filter for this object.')

    print(f'convert_1Dspectra_fits_to_dat completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

# -------------------------------------------------------------------------------------------------------
def make_1D_spectra_per_orientation(args):
    '''
    Makes *1D_C.dat and *1D_R.dat spectra files, corresponding to each grism orientations, from *beams.fits files, for a given field
    Adapted by Ayan Acharyya from passage_analysis/passage_convert_data.ipynb
    '''
    print(f'Starting make_1D_spectra_per_orientation..')
    start_time = datetime.now()

    speccat_file = glob.glob(args.speccat_file_path + '*_speccat.fits')[0]
    speccat = Table.read(speccat_file)

    beam_files = sorted(glob.glob(args.beam_files_path + '*beams.fits'))

    # ------looping over all beam files--------------
    for index, beam_filename in enumerate(beam_files):
        print(f'\nDoing {index + 1} of {len(beam_files)} files..')
        objid = int(os.path.basename(beam_filename).split('_')[1].split('.')[0])
        outfile_thisobj = glob.glob(args.spectra_path + f'Par{args.parno}_{objid:05d}.*_1D_*.dat')

        if len(outfile_thisobj) == 0: # if no R/C file exists for this object id
            z = speccat[speccat['id'] == objid]['redshift'].value[0]

            mb = multifit.MultiBeam(beam_filename, fcontam=0.1, sys_err=0.02, min_sens=0.05, MW_EBV=-1, group_name='', verbose=False)

            Cgrism_beams = [mb.beams[k] for k in range(len(mb.beams)) if mb.beams[k].grism.filter == 'GR150C']
            Rgrism_beams = [mb.beams[k] for k in range(len(mb.beams)) if mb.beams[k].grism.filter == 'GR150R']

            if len(Cgrism_beams) > 0:
                mb_C = multifit.MultiBeam(beams=Cgrism_beams, fcontam=0.1, sys_err=0.02, min_sens=0.05, MW_EBV=-1, group_name='')

                # this catches cases where spectrum contains only zeros (very rare)
                # The fit will crash in such cases
                try:
                    tfitC = mb_C.template_at_z(z, fitter='bounded')
                    keys_C = mb_C.oned_spectrum(tfit=tfitC, bin=1).keys()
                except:
                    keys_C = mb_C.oned_spectrum(bin=1).keys()

                for c in keys_C:
                    t_out = make_table(mb_C.oned_spectrum(tfit=tfitC, bin=1)[c])
                    outfilename_C = args.spectra_path + os.path.basename(beam_filename).replace('beams.fits', c + '_1D_C.dat')
                    t_out.write(outfilename_C, format='ascii.fixed_width_two_line', overwrite=True)
                    print(f'Written {outfilename_C}')

            if len(Rgrism_beams) > 0:
                mb_R = multifit.MultiBeam(beams=Rgrism_beams, fcontam=0.1, sys_err=0.02, min_sens=0.05, MW_EBV=-1, group_name='')

                # this catches cases where spectrum contains only zeros (very rare)
                # The fit will crash in such cases
                try:
                    tfitR = mb_R.template_at_z(z, fitter='bounded')
                    keys_R = mb_R.oned_spectrum(tfit=tfitR, bin=1).keys()
                except:
                    keys_R = mb_R.oned_spectrum(bin=1).keys()

                for r in keys_R:
                    t_out = make_table(mb_R.oned_spectrum(tfit=tfitR, bin=1)[r])
                    outfilename_R = args.spectra_path + os.path.basename(beam_filename).replace('beams.fits', c + '_1D_R.dat')
                    t_out.write(outfilename_R, format='ascii.fixed_width_two_line', overwrite=True)
                    print(f'Written {outfilename_R}')
        else:
            print(f'It appears the R/C files were already created for this object. Skipping object {objid:05d}.')

    print(f'make_1D_spectra_per_orientation completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

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
    args = get_user_directory_structure(args)

    # -----pad field name if needed------
    if args.is_fieldname_padded: parno = f'{parno:03}'
    args.parno = parno

    # -----convert pseudo paths to proper paths------
    args = substitute_fieldname_in_paths(args, str(parno))

    # -----import passage_analysis if not successful earlier------
    if 'passage' not in locals():
        sys.path.insert(1, args.code_dir)
        import passage_analysis as passage
        from passage_analysis import utilities
        print('passage_analysis now successfully imported.')

    # ---------check if region files exist. If not, run code to create necessary region files---------
    regionfiles = glob.glob(args.region_file_path + '*.reg')
    if len(regionfiles) == 0:
        print('\033[94m' + "No region files found, creating those for you now."  + '\033[0m')
        utilities.create_regions(parno, args)
    else:
        print(f'\nFound region files in {args.region_file_path}, proceeding..')

    # ---------check if .dat spectra files exist. If not, make them------------------
    spec_files = sorted(glob.glob(args.spectra_path + '*1D.dat'))
    if len(spec_files) == 0 or args.clobber_1D:
        convert_1Dspectra_fits_to_dat(args)
    else:
        print(f'\nFound 1D.dat files in {args.spectra_path}, proceeding..')

    # ---------check if R and C spectra files exist. If not, make them------------------
    C_files = sorted(glob.glob(args.spectra_path + '*C.dat'))
    R_files = sorted(glob.glob(args.spectra_path + '*R.dat'))
    if (len(C_files) == 0 and len(R_files) == 0) or args.clobber_RC:
        make_1D_spectra_per_orientation(args)
    else:
        print(f'\nFound R & C files in {args.spectra_path}, proceeding..')

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


