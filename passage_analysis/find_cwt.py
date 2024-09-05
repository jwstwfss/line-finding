import os
import pdb
from glob import glob
from passage_analysis import *
from astropy.table import Table
from astropy.io import ascii as asc
from astropy.io import fits
from datetime import datetime, timedelta
import pandas as pd

def find_cwt(lam, flux, err, fwhm_est_pix, beam, config_pars, plotflag = True):
    '''
    AA removed 'zeros' from the function definition as it seemed to be an unused variable
    '''

    cont_medfilt = int(config_pars['cont_medfilt']) # window for median filter
    max_width = config_pars['maxwidth'] * fwhm_est_pix
    min_width = config_pars['minwidth']
    dw  = (max_width - min_width) / config_pars['nwidths']
    widths = min_width + np.arange(config_pars['nwidths']) * dw
    max_distance_ridge = widths * config_pars['max_dist_ridge_scl'] + config_pars['max_dist_ridge_const'] # if peak in a row of the cwt matrix is off by more than this many pixels new/separate ridge
    gap_allowed_between_ridges = config_pars['gap_btw_ridges'] # gap between ridges can be no more than N pixels, otherwise new/separate ridge/peak
    snr_cwt = config_pars['snr_cwt'] # snr for the cwt line finder
    noise_cut_cwt = config_pars['noise_cut_cwt'] # noise cut for cwt to estimate noise
    min_length_cwt = config_pars['min_length_ridge'] # minimum length of a cwt ridge to be considered real
    edge_reject = config_pars['edge_reject'] # reject cwt detections within 5 pixcels of edge
    sn_thresh_cont_check = config_pars['n_sigma_above_cont'] # step2 requires cwt line candidates to have npix_thresh abvove sn_thresh_cont_check
    npix_thresh = config_pars['npix_thresh']
    min_line_contrast = config_pars['min_line_contrast'] # minimum allowed for rejecting low EW lines.
    
    # print('this is the max:', config_pars['lambda_max'])

    if plotflag == True:
        f, axarr = plt.subplots(2, 1, figsize=(8, 8))
        w=np.where((lam > config_pars['lambda_min']) & (lam < config_pars['lambda_max']))
        spec_max  = np.max(flux[w])
        axarr[1].plot(lam,flux, ls='steps-mid', color = 'k')
        axarr[1].axis([np.min(lam), np.max(lam), -0.5e-19, 1.3 * spec_max])

    # run the order filter
    cont_filter = si.medfilt(flux, cont_medfilt)

    if plotflag ==True :
        # continuum model
        axarr[1].plot(lam, cont_filter+err*sn_thresh_cont_check, color= 'orange')
        axarr[1].plot(lam, cont_filter)

    # calculate and show continuous wavelet transform array
    cwarray = si.cwt(flux, si.ricker, widths)

    if plotflag ==True:
        axarr[0].imshow(cwarray, vmax = .6 * spec_max, vmin  = -0.2 * spec_max, aspect = 'auto')

    # find the peaks and overplot them
    peaks= si.find_peaks_cwt(flux, widths, wavelet = si.ricker, max_distances=max_distance_ridge, gap_thresh=gap_allowed_between_ridges,
               min_snr=snr_cwt, noise_perc = noise_cut_cwt, min_length=min_length_cwt)

    if plotflag == True:
        axarr[1].plot(lam[peaks], flux[peaks], 'ro', ms=7)

    peaks = np.array(peaks)

    # reject peaks near the edge
    w = np.where((peaks > edge_reject) & (peaks < np.size(lam) - edge_reject))
    peaks = peaks[w[0]]

    # KVN testing:
    #print('peaks: ', peaks)
    #print('division: ', (flux[peaks] - cont_filter[peaks])/cont_filter[peaks])
    
    # reject lines with presumably low EWs
    try:
        peak_contrast = (flux[peaks] - cont_filter[peaks])/cont_filter[peaks]
    except IndexError:
        # added this step in case peaks fail the above conditional
        peaks = []
    else:
        w = np.where(peak_contrast > min_line_contrast)
        peaks = peaks[w[0]]

    if np.size(peaks) > 0:

        # count contiguous pixels above the noise threshold
        snr_thresh = sn_thresh_cont_check
        npix_peak = []
        line_snr_guess = []

        for i in peaks:
            # first of all, is peak above threshold
            if flux[i] > cont_filter[i] + err[i] * snr_thresh:
                pixel_count = 1
                cond = 0
                j = i + 1

                while ((cond == 0.) & (j < np.size(flux) - 1)):
                    if flux[j] > cont_filter[j] + err[j] * snr_thresh:
                        pixel_count = pixel_count + 1
                        j = j + 1
                    else:
                        cond = 1.

                cond = 0
                j = i-1

                while ((cond == 0) & (j > 0)):
                    if flux[j] > cont_filter[j] + err[j] * snr_thresh:
                         pixel_count = pixel_count + 1
                         #cond = 0
                         j = j-1
                    else:
                        cond = 1.
            else:
                pixel_count = 0

            npix_peak.append(pixel_count)

            # crudely estimate the snr of each line candidate
            if lam[i] > config_pars['transition_wave1']:
                disp_est = config_pars['dispersion_red']
            else:
                disp_est = config_pars['dispersion_blue']

            fwhm_est = fwhm_est_pix * disp_est

            w = np.where((lam > lam[i] - 0.5 * fwhm_est) & (lam < lam[i] + 0.5 * fwhm_est))
            line_signal_guess = np.sum((flux[w] - cont_filter[w]) * disp_est)
            line_noise_guess = np.sqrt(np.sum((err[w] * disp_est)**2))
            line_snr_guess.append(line_signal_guess/line_noise_guess)

        npix_peak = np.array(npix_peak)
        line_snr_guess = np.array(line_snr_guess)
        w = np.where(npix_peak >= npix_thresh)
        real_peaks = peaks[w[0]]
        npix_real = npix_peak[w[0]]
        snr_real = line_snr_guess[w[0]]
    else:
         real_peaks = []
         npix_real = []
         snr_real = []
         peaks = []

    if plotflag==True:
            axarr[1].plot(lam[real_peaks], flux[real_peaks], 'rs', ms=9, markeredgecolor= 'r', markerfacecolor = 'none', markeredgewidth=2)
            plt.title(f'{beam}')
            plt.show(block=True)

    return [lam[real_peaks], flux[real_peaks], npix_real, snr_real, cwarray, cont_filter, lam[peaks], flux[peaks]]

# -------------------------------------------------------------------------------------------------------------------
def loop_field_cwt(parno, args):
    '''
    Function to create linelist for a given PASSAGE field; incorporates each available filter for each available object
    Saves linelist as ASCII file

    Original code by M.D.R and K.V.N
    KVN updating this to write the linelist to the 'output' directory... take path to data as input
    AA updated this, for flexibility, to take args as input, which contains the directory structure
    AA further modified this whole function to make it simpler and more concise
    '''

    start_time = datetime.now()
    filters = ['115', '150', '200'] # AA added on 3 Sep 2024

    # -------displaying/creating relevant paths----------------
    print(f'\nSearching for default.config at: {args.code_dir}')
    print(f'\nSearching for grism spectra files in {args.spectra_path}')

    print(f'\nSearching for catalogs at: {args.photcat_file_path}')
    catalogs = glob(args.photcat_file_path + f'Par{parno}_phot*.fits') # get list of available catalogs
    catalogs.sort()
    print(f'\nI found the following catalogs: {catalogs}..')
    cat = Table.read(catalogs[0])
    print(f'..and successfully opened the first one!')
    # M.D.R. - 10/08/2020

    if not os.path.exists(args.linelist_path): os.mkdir(args.linelist_path)
    outfilename = args.linelist_path + f'Par{parno}lines.dat'
    output_df = pd.DataFrame()

    # -------looping over all three PASSAGE filters----------------
    for index, thisfilter in enumerate(filters):
        print(f'Doing filter G{thisfilter} which is {index + 1} out of {len(filters)} filters..')

        thisfilter_files = glob(args.spectra_path + f'*{thisfilter}_1D.dat')
        thisfilter_files.sort()

        # ------------looping over all spectra files of the current filter------------
        for index2, filename in enumerate(thisfilter_files):
            print(f'Starting obj id = {os.path.basename(filename)} which is {index2 + 1} of {len(thisfilter_files)}..')

            # -------read in config_pars freshly for each iteration (in case it gets changed within this loop)----------
            config_pars = read_config(str(args.code_dir)+'/default.config')
            config_pars['transition_wave1'] = 13000. # MDR 2022/08/16

            # ------------get spectral data-----------------
            spdata = asc.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])

            if '115' in thisfilter: trimmed_spec = trim_spec(spdata, None, None, config_pars)
            elif '150' in thisfilter: trimmed_spec = trim_spec(None, spdata, None, config_pars)
            elif '200' in thisfilter: trimmed_spec = trim_spec(None, None, spdata, config_pars)

            trimmed_df = pd.DataFrame({'lambda': trimmed_spec[0], 'flux':trimmed_spec[1], 'error':trimmed_spec[2], 'contam':trimmed_spec[3], 'zero':trimmed_spec[4]}) # converting into a dataframe for easier handling and readability
            if len(trimmed_df) < config_pars['min_spec_length']: continue # check if spectrum is long enough to proceed
            trimmed_df['flux_corrected'] = trimmed_df['flux'] - trimmed_df['contam']

            # -----------------look up the object in the se catalog and grab the a_image-----------------
            beam = int(os.path.basename(filename).split('_')[1].split('.')[0])
            fwhm_est_pix = cat[cat['id'] == beam]['a_image'][0] * 2.0

            # -----------------cwt it and unpack and append the results to dataframe-----------------
            thisfilter_cwt = find_cwt(trimmed_df['lambda'].values, trimmed_df['flux_corrected'].values, trimmed_df['error'].values, fwhm_est_pix, beam, config_pars, plotflag=False)
            cwt_df = pd.DataFrame({'lambda': thisfilter_cwt[0], 'flux':thisfilter_cwt[1], 'npix':thisfilter_cwt[2], 'snr':thisfilter_cwt[3], 'beam':beam, 'filter':thisfilter}) # converting into a dataframe for easier handling and readability
            output_df = pd.concat([output_df, cwt_df])

            # -----------------repeat above cwt step if needed-----------------
            if config_pars['n_sigma_for_2pix_lines']:
                thisfilter_cwt = find_cwt(trimmed_df['lambda'].values, trimmed_df['flux_corrected'].values, trimmed_df['error'].values, fwhm_est_pix, beam, config_pars, plotflag=False)
                cwt_df = pd.DataFrame({'lambda': thisfilter_cwt[0], 'flux': thisfilter_cwt[1], 'npix': thisfilter_cwt[2], 'snr': thisfilter_cwt[3], 'beam':beam, 'filter':thisfilter})  # converting into a dataframe for easier handling and readability
                output_df = pd.concat([output_df, cwt_df])

    # -----------------save the final linelist file----------------------
    output_df['par'] = parno
    output_df = output_df.sort_values(['beam', 'filter', 'lambda'])
    output_df['filter'] = output_df['filter'].map(lambda x: f'G{x}')
    output_df = output_df[['par', 'filter', 'beam', 'lambda', 'npix', 'snr']] # to re-order columsn before saving to file
    output_df.to_csv(outfilename, sep='\t', index=None, header=False)

    print(f'loop_find_cwt() completed in {timedelta(seconds=(datetime.now() - start_time).seconds)}')
