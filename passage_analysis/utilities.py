import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import glob
import os

from datetime import datetime, timedelta
from astropy.io import fits
from astropy.table import Table
from grizli import multifit

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
    
# -------------------------------------------------------------------------------------------------------
def make_table(tab):
    '''
    Modifies and returns input table, by changing normalisation of flux, etc
    Adapted by AA from passage_analysis/passage_convert_data.ipynb; Sep 2024
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
    Adapted by AA from passage_analysis/passage_convert_data.ipynb; Sep 2024
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

                if not os.path.exists(outfilename):
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
    Adapted by AA from passage_analysis/passage_convert_data.ipynb; Sep 2024
    '''
    print(f'Starting make_1D_spectra_per_orientation..')
    start_time = datetime.now()

    speccat_file = glob.glob(args.speccat_file_path + '*_speccat.fits')[0]
    speccat = Table.read(speccat_file)

    beam_files = sorted(glob.glob(args.beam_files_path + '*beams.fits'))

    # ------looping over all beam files--------------
    for index, beam_filename in enumerate(beam_files):
        objid = int(os.path.basename(beam_filename).split('_')[1].split('.')[0])
        print(f'\nDoing object {objid} which is {index + 1} of {len(beam_files)} files..')
        outfile_thisobj = glob.glob(args.spectra_path + f'Par{args.parno}_{objid:05d}.*_1D_*.dat')

        if len(outfile_thisobj) == 0: # if no R/C file exists for this object id
            try:
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
            except Exception as e:
                print(f'Skipping object {objid} because it ran into error {e}'.)
                continue
        else:
            print(f'It appears the R/C files were already created for this object. Skipping object {objid:05d}.')

    print(f'make_1D_spectra_per_orientation completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')

