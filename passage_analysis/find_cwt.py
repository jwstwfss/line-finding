import os
import pdb
from glob import glob
from passage_analysis import *
from astropy.table import Table
from astropy.io import ascii as asc

def find_cwt(lam, flux, err, fwhm_est_pix, beam_name, config_pars, plotflag = True):
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
    print('peaks: ', peaks)
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
            plt.title(beam_name)
            plt.show(block=True)

    return [lam[real_peaks], flux[real_peaks], npix_real, snr_real, cwarray, cont_filter, lam[peaks], flux[peaks]]

def loop_field_cwt(parno, args):
    # no inputs and run from inside the data directory
    # KVN updating this to write the linelist to the 'output' directory... take path to data as input
    # AA updated this, for flexibility, to take args as input, which contains the directory structure
    if not os.path.exists(args.linelist_path):
        os.mkdir(args.linelist_path)

    print('Looking for spectra here: ', args.spec1D_path)
    g115files = glob(args.spec1D_path + '*G115_1D.dat') # looking for 3 spectra for PASSAGE
    g115files.sort()
    g150files = glob(args.spec1D_path + '*G150_1D.dat')
    g150files.sort()

    # M.D.R. - 10/08/2020
    print(f'\nSearching for default.config at: {args.code_dir}')
    config_pars = read_config(str(args.code_dir) + '/default.config')

    photcat_file_path = args.photcat_file_path + f'Par{parno}_phot*.fits'
    print(f'Searching for catalogs at: {photcat_file_path}')
    catalogs = glob(photcat_file_path) # get list of available catalogs
    catalogs.sort()

    print(f'\nI found the following catalogs: {catalogs}')
    cat = Table.read(catalogs[0])

    print(f'\nCatalog opened successfully: {catalogs[0]}')
    # M.D.R. - 10/08/2020

    a_images = cat['a_image']
    beam_se = cat['id']

    tempfilename = args.linelist_path + 'temp'
    config_pars['transition_wave'] = 13000. # MDR 2022/08/16

    print('\nSearching for grism files...')

    with open(tempfilename, 'w') as outfile:
        # AA added looping over all three PASSAGE filters, as opposed separate code block for each filter
        filters = ['115', '150', '200']
        for index, thisfilter in enumerate(filters):
            print(f'Doing filter G{thisfilter} which is {index + 1} out of {len(filters)} filters..')
            thisfilter_files = glob(args.spec1D_path + f'*{thisfilter}_1D.dat')
            thisfilter_files.sort()

            # looping over all spectra files of the current filter
            for filename in thisfilter_files:
                print('starting obj id = ', filename)
                # get spectral data
                spdata = asc.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])
                trimmed_spec = trim_spec(spdata, None, None, config_pars)

                # look up the object in the se catalog and grab the a_image
                beam = float(filename.split('_')[1].split('.')[0])
                parno = parno #os.getcwd().split('/')[-2].split('Par')[-1] # fixed parallel field number to zero for the mudf program
                print('Par Number: ', parno)

                w = np.where(beam_se == beam)
                w = w[0] # because of tuples
                a_image = a_images[w][0]
                fwhm_est_pix = a_image * 2.0

                # unpack spectrum and check that it is long enough to proceed
                lam = trimmed_spec[0]
                flux_corr = trimmed_spec[1] - trimmed_spec[3]
                err = trimmed_spec[2]

                if len(lam) < config_pars['min_spec_length']:
                    continue

                # cwt it and unpack and write results
                thisfilter_cwt = find_cwt(lam, flux_corr, err, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
                lam_cwt = thisfilter_cwt[0]
                flam_cwt = thisfilter_cwt[1]
                npix_cwt = thisfilter_cwt[2]
                snr_cwt = thisfilter_cwt[3]

                for i in np.arange(len(lam_cwt)):
                    print(beam, f'G{thisfilter}', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
                    outfile.write(f'{parno}  {thisfilter}  {int(beam)}  {lam_cwt[i]}  {npix_cwt[i]}  {snr_cwt[i]}\n')

                if config_pars['n_sigma_for_2pix_lines'] != False:
                    config_pars['npix_thresh'] = 2
                    config_pars['n_sigma_above_cont'] = config_pars['n_sigma_for_2pix_lines']
                    thisfilter_cwt = find_cwt(lam, flux_corr, err, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
                    lam_cwt = thisfilter_cwt[0]
                    flam_cwt = thisfilter_cwt[1]
                    npix_cwt = thisfilter_cwt[2]
                    snr_cwt = thisfilter_cwt[3]
                    for i in np.arange(len(lam_cwt)):
                        print(beam, f'G{thisfilter}', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
                        outfile.write(f'{parno}  {thisfilter}  {int(beam)}  {lam_cwt[i]}  {npix_cwt[i]}  {snr_cwt[i]}\n')

                # go back to the beginning with the old config pars
                config_pars = read_config(str(args.code_dir)+'/default.config')
                config_pars['transition_wave1'] = 13000. # MDR 2022/08/16

    tab = asciitable.read(tempfilename, format = 'no_header')
    par = tab['col1']
    grism = tab['col2']
    beam = tab['col3']
    wave = tab['col4']
    npix = tab['col5']
    snr = tab['col6']
    s = np.argsort(beam)
    beam = beam[s]
    grism = grism[s]
    wave =wave[s]
    npix = npix[s]
    snr = snr[s]
    par = par[0]
    beams_unique = np.unique(beam)
    outfilename = args.linelist_path + f'Par{par}lines.dat'

    with open(outfilename, 'w') as outfile:
        for b in beams_unique:
            # AA adding a loop over filters, instead of separate blocks of code per filter
            for thisfilter in filters:
                w = (beam == b) & (grism == f'G{thisfilter}')
                waves_obj = wave[w]
                npix_obj = npix[w]
                snr_obj = snr[w]
                waves_uniq, ind = np.unique(waves_obj, return_index = True)
                npix_uniq = npix_obj[ind]
                snr_uniq = snr_obj[ind]
                s = np.argsort(waves_uniq)
                waves_final_thisfilter = waves_uniq[s]
                npix_final_thisfilter = npix_uniq[s]
                snr_final_thisfilter = snr_uniq[s]

                for lam, npx, sn in zip(waves_final_thisfilter, npix_final_thisfilter, snr_final_thisfilter):
                    outfile.write(f'{par}  G{thisfilter}  {b}  {lam} {npx} {sn}\n')
