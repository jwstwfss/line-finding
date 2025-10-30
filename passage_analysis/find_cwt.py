import os
import pdb
from glob import glob
from passage_analysis import *
from astropy.table import Table
from astropy.io import ascii as asc
from copy import copy

def find_cwt(lam, flux, err, zeros, fwhm_est_pix, beam_name, config_pars, plotflag = True):

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
    edge_reject = config_pars['edge_reject'] # reject cwt detections within 5 pixels of edge
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
            plt.title(beam_name)
            plt.show(block=True)

    return [lam[real_peaks], flux[real_peaks], npix_real, snr_real, cwarray, cont_filter, lam[peaks], flux[peaks]]

def loop_field_cwt(path_to_data, path_to_code, parno):
    # no inputs and run from inside the data directory
    # KVN updating this to write the linelist to the 'output' directory... take path to data as input
    if not os.path.exists(os.path.join(path_to_data, 'linelist')):
        os.system(f"mkdir -p {os.path.join(path_to_data, 'linelist'):s}")

    print('Looking for spectra here: ', str(path_to_data)+'Par'+str(parno)+'/Spectra/*G115_1D.dat')
    g115files = glob(str(path_to_data)+'Par'+str(parno)+'/Spectra/*G115_1D.dat') # looking for 3 spectra for PASSAGE
    g115files.sort()
    g150files = glob(str(path_to_data)+'Par'+str(parno)+'/Spectra/*G150_1D.dat')
    g150files.sort()
    g200files = glob(str(path_to_data)+'Par'+str(parno)+'/Spectra/*G200_1D.dat')
    g200files.sort()

    # M.D.R. - 10/08/2020
    print('')
    print('Searching for default.config at: ' + str(path_to_code))
    config_pars = read_config(str(path_to_code)+'/default.config')
    
    print('Searching for catalogs at: ' + str(path_to_data) +'Par'+str(parno)+'/DATA/DIRECT_GRISM/Par'+str(parno)+'_phot*.fits')
    catalogs = glob(str(path_to_data) +'Par'+str(parno)+'/DATA/DIRECT_GRISM/Par'+str(parno)+'_phot*.fits') # get list of available catalogs
    catalogs.sort()
    print('')
    print('I found the following catalogs: ' + str(catalogs))
    cat = Table.read(catalogs[0])
    print('')
    print('Catalog opened successfully: ' + str(catalogs[0]))

    if len(catalogs) > 1:
        cat2 = Table.read(catalogs[1])

    a_images = cat['a_image']
    beam_se = cat['id']
    
    outfile = open(os.path.join(path_to_data, 'linelist/temp'), 'w')
    #config_pars['transition_wave'] = 11700.
    config_pars['transition_wave'] = 13000. # MDR 2022/08/16

    print('')
    print('Searching for grism files...')

    for filename in g115files:
        print(filename)
        print('starting obj id = ', filename)
        # get spectral data
        # spdata = asc.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])
        spdata = Table.read(filename, format="ascii")
        col_names = spdata.colnames
        if 'err' in col_names and 'ferr' not in col_names:
            spdata.rename_columns(["wave", "error", "zeroth"], ["lambda", "ferror", "zero"])
        elif 'ferr' in col_names:
            spdata.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        trimmed_spec = trim_spec(spdata, None, None, config_pars)

        # look up the object in the se catalog and grab the a_image
        # beam = float(filename.split('_')[1].split('.')[0])
        beam = float(filename.split('Spectra/Par')[1].split('_')[1].split('.')[0])
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
        zeros = trimmed_spec[4]

        if len(lam) < config_pars['min_spec_length']:
            continue

        # cwt it and unpack and write results
        g115_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
        lam_cwt = g115_cwt[0]
        flam_cwt = g115_cwt[1]
        npix_cwt = g115_cwt[2]
        snr_cwt = g115_cwt[3]

        for i in np.arange(len(lam_cwt)):
            #print(beam, 'G102', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
            print(beam, 'G115', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
            outfile.write(str(parno)+'  G115  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) + '  ' + str(npix_cwt[i]) + '  ' + str(snr_cwt[i]) + '\n')

        if config_pars['n_sigma_for_2pix_lines'] != False:
            config_pars['npix_thresh'] = 2
            config_pars['n_sigma_above_cont'] = config_pars['n_sigma_for_2pix_lines']
            g115_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
            lam_cwt = g115_cwt[0]
            flam_cwt = g115_cwt[1]
            npix_cwt = g115_cwt[2]
            snr_cwt = g115_cwt[3]
            for i in np.arange(len(lam_cwt)):
                #print(beam, 'G102', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
                print(beam, 'G115', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
                outfile.write(str(parno)+'  G115  ' +str(int(beam)) + '  ' + str(lam_cwt[i]) + '  ' + str(npix_cwt[i]) + '  ' + str(snr_cwt[i]) + '\n')

        # go back to the beginning with the old config pars
        config_pars = read_config(str(path_to_code)+'/default.config')
        #config_pars['transition_wave'] = 11700.
        config_pars['transition_wave1'] = 13000. # MDR 2022/08/16

    #config_pars['transition_wave'] = 11100.
    config_pars['transition_wave1'] = 13000. # MDR 2022/08/16

    for filename in g150files:
        print('starting obj id = ', filename)
        # spdata = asc.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])
        spdata = Table.read(filename, format="ascii")
        col_names = spdata.colnames
        if 'err' in col_names and 'ferr' not in col_names:
            spdata.rename_columns(["wave", "error", "zeroth"], ["lambda", "ferror", "zero"])
        elif 'ferr' in col_names:
            spdata.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        trimmed_spec = trim_spec(None, spdata, None, config_pars)
        # beam = float(filename.split('_')[1].split('.')[0])
        beam = float(filename.split('Spectra/Par')[1].split('_')[1].split('.')[0])
        parno = parno
        print('Par Number: ', parno)
        w = np.where(beam_se == beam)
        w = w[0]    # because of tuples
        a_image = a_images[w][0]
        lam = trimmed_spec[0]
        flux_corr = trimmed_spec[1] - trimmed_spec[3]
        err = trimmed_spec[2]
        zeros = trimmed_spec[4]

        if len(lam) < config_pars['min_spec_length']:
            continue
        fwhm_est_pix = a_image * 2
        config_pars

        g150_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix,str(int(beam)), config_pars, plotflag=False)
        lam_cwt = g150_cwt[0]
        flam_cwt = g150_cwt[1]
        npix_cwt = g150_cwt[2]
        snr_cwt = g150_cwt[3]

        for i in np.arange(len(lam_cwt)):
            print(beam, 'G150', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
            outfile.write(str(parno)+'  G150  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) + '  ' + str(npix_cwt[i]) + '  ' + str(snr_cwt[i]) + '\n')

        if config_pars['n_sigma_for_2pix_lines'] != False:
            config_pars['npix_thresh'] = 2
            config_pars['n_sigma_above_cont'] = config_pars['n_sigma_for_2pix_lines']
            g150_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
            lam_cwt = g150_cwt[0]
            flam_cwt = g150_cwt[1]
            npix_cwt = g150_cwt[2]
            snr_cwt = g150_cwt[3]

            for i in np.arange(len(lam_cwt)):
                print(beam, 'G150', lam_cwt[i], npix_cwt[i], snr_cwt[i])
                outfile.write(str(parno)+'  G150  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) + '  ' + str(npix_cwt[i]) + '  ' + str(snr_cwt[i]) + '\n')

        # go back to the beginning with the old config pars
        config_pars = read_config(str(path_to_code)+'/default.config')
        #config_pars['transition_wave'] = 11200.
        config_pars['transition_wave'] = 13000. # MDR 2022/08/16

    for filename in g200files:
        print('starting obj id = ', filename)
        # spdata = asc.read(filename, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])
        spdata = Table.read(filename, format="ascii")
        col_names = spdata.colnames
        print(col_names);print('err' in col_names); print('ferr' in col_names)
        if 'err' in col_names and 'ferr' not in col_names:
            spdata.rename_columns(["wave", "error", "zeroth"], ["lambda", "ferror", "zero"])
        elif 'ferr' in col_names:
            spdata.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        trimmed_spec = trim_spec(None, None, spdata, config_pars)
        # beam = float(filename.split('_')[1].split('.')[0])
        beam = float(filename.split('Spectra/Par')[1].split('_')[1].split('.')[0])
        parno = parno # fixed parallel field number to zero for the mudf program
        print('Par Number: ', parno)
        w = np.where(beam_se == beam)
        w = w[0]    # because of tuples
        a_image = a_images[w][0]
        lam = trimmed_spec[0]
        flux_corr = trimmed_spec[1] - trimmed_spec[3]
        err = trimmed_spec[2]
        zeros = trimmed_spec[4]

        if len(lam) < config_pars['min_spec_length']:
            continue
        fwhm_est_pix = a_image * 2
        config_pars

        g200_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix,str(int(beam)), config_pars, plotflag=False)
        lam_cwt = g200_cwt[0]
        flam_cwt = g200_cwt[1]
        npix_cwt = g200_cwt[2]
        snr_cwt = g200_cwt[3]

        for i in np.arange(len(lam_cwt)):
            print(beam, 'G200', lam_cwt[i], npix_cwt[i], fwhm_est_pix, snr_cwt[i])
            outfile.write(str(parno)+'  G200  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) + '  ' + str(npix_cwt[i]) + '  ' + str(snr_cwt[i]) + '\n')

        if config_pars['n_sigma_for_2pix_lines'] != False:
            config_pars['npix_thresh'] = 2
            config_pars['n_sigma_above_cont'] = config_pars['n_sigma_for_2pix_lines']
            g200_cwt= find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(int(beam)), config_pars, plotflag=False)
            lam_cwt = g200_cwt[0]
            flam_cwt = g200_cwt[1]
            npix_cwt = g200_cwt[2]
            snr_cwt = g200_cwt[3]

            for i in np.arange(len(lam_cwt)):
                print(beam, 'G200', lam_cwt[i], npix_cwt[i], snr_cwt[i])
                outfile.write(str(parno)+'  G200  ' + str(int(beam)) + '  ' + str(lam_cwt[i]) + '  ' + str(npix_cwt[i]) + '  '+ str(snr_cwt[i]) + '\n')

        # go back to the beginning with the old config pars
        config_pars = read_config(str(path_to_code)+'/default.config')
        #config_pars['transition_wave'] = 11200.
        config_pars['transition_wave'] = 13000. # MDR 2022/08/16

    outfile.close()

    tab = asciitable.read(os.path.join(path_to_data, 'linelist/temp'), format = 'no_header')
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
    outfile = open(os.path.join(path_to_data, 'linelist/Par'+str(parno) + 'lines.dat'), 'w')

    for b in beams_unique:
        # do the g102 for b
        w = (beam == b) & (grism == 'G115')
        waves_obj = wave[w]
        npix_obj = npix[w]
        snr_obj = snr[w]
        waves_uniq, ind = np.unique(waves_obj, return_index = True)
        npix_uniq = npix_obj[ind]
        snr_uniq = snr_obj[ind]
        s = np.argsort(waves_uniq)
        waves_final_g115 = waves_uniq[s]
        npix_final_g115 = npix_uniq[s]
        snr_final_g115 = snr_uniq[s]

        for lam, npx, sn in zip(waves_final_g115, npix_final_g115, snr_final_g115):
            outfile.write(str(parno) + '  G115  ' + str(b) + '  ' + str(lam) + '  ' + str(npx) + '  ' + str(sn) + '\n')

        # do the g141 for b  
        w = (beam == b) & (grism == 'G150')
        waves_obj = wave[w]
        npix_obj = npix[w]
        snr_obj = snr[w]
        waves_uniq, ind = np.unique(waves_obj, return_index = True)
        npix_uniq = npix_obj[ind]
        snr_uniq = snr_obj[ind]
        s = np.argsort(waves_uniq)
        waves_final_g150 = waves_uniq[s]
        npix_final_g150 = npix_uniq[s]
        snr_final_g150 = snr_uniq[s]

        for lam, npx, sn in zip(waves_final_g150, npix_final_g150, snr_final_g150):
            outfile.write(str(parno) + '  G150  ' + str(b) + '  ' + str(lam) + '  ' + str(npx) + '  ' + str(sn) + '\n')

        # KVN: Adding third grism filter: do the g200 for b  
        w = (beam == b) & (grism == 'G200')
        waves_obj = wave[w]
        npix_obj = npix[w]
        snr_obj = snr[w]
        waves_uniq, ind = np.unique(waves_obj, return_index = True)
        npix_uniq = npix_obj[ind]
        snr_uniq = snr_obj[ind]
        s = np.argsort(waves_uniq)
        waves_final_g200 = waves_uniq[s]
        npix_final_g200 = npix_uniq[s]
        snr_final_g200 = snr_uniq[s]

        for lam, npx, sn in zip(waves_final_g200, npix_final_g200, snr_final_g200):
            outfile.write(str(parno) + '  G200  ' + str(b) + '  ' + str(lam) + '  ' + str(npx) + '  ' + str(sn) + '\n')

    outfile.close()

def test_obj_cwt(parno, beamno, configfile):
    blue_se = asciitable.read('../DATA/DIRECT_GRISM/fin_F110.cat')
    red_se  = asciitable.read('../DATA/DIRECT_GRISM/fin_F160.cat')
    a_image_blue = blue_se['col5']
    a_image_red = red_se['col5']
    beam_se = blue_se['col2']
    config_pars = read_config(configfile)
    bluefile = 'Par'+str(parno) + '_G102_BEAM_'+str(beamno)+'A.dat'
    redfile =  'Par'+str(parno) + '_G141_BEAM_'+str(beamno)+'A.dat'

    spdata_blue = asciitable.read(bluefile, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])
    trimmed_spec_blue= trim_spec(spdata_blue, None, config_pars)

    # do the blue side
    lam = trimmed_spec_blue[0]
    flux_corr = trimmed_spec_blue[1] - trimmed_spec_blue[3]
    err = trimmed_spec_blue[2]
    zeros = trimmed_spec_blue[4]
    #config_pars['transition_wave']  = 11700.
    config_pars['transition_wave'] = 13000. # MDR 2022/08/16

    if len(lam) < config_pars['min_spec_length']:
        print('Short spec. skip it!')
    else:
        w = np.where(beam_se == beamno)
        w = w[0] # because of tuples
        a_image = a_image_blue[w][0]
        fwhm_est_pix = a_image * 2
        g102_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(beamno), config_pars, plotflag=True)
        print(g102_cwt[0], g102_cwt[1], g102_cwt[2], fwhm_est_pix)

    # do the red side
    #config_pars['transition_wave'] = 11200.
    #config_pars['transition_wave'] = 9353. # MDR 2022/08/16
    spdata_red = asciitable.read(redfile, names = ['lambda', 'flux', 'ferror', 'contam', 'zero'])
    trimmed_spec_red = trim_spec(None, spdata_red, config_pars)
    lam = trimmed_spec_red[0]
    flux_corr = trimmed_spec_red[1] - trimmed_spec_red[3]
    err = trimmed_spec_red[2]
    zeros = trimmed_spec_red[4]

    if len(lam) < config_pars['min_spec_length']:
        print('Short spec. skip it!')
    else:
        w = np.where(beam_se == beamno)
        w = w[0] # because of tuples
        a_image = a_image_red[w][0]
        fwhm_est_pix = a_image * 2
        g141_cwt = find_cwt(lam, flux_corr, err, zeros, fwhm_est_pix, str(beamno), config_pars, plotflag=True)
        print(g141_cwt[0], g141_cwt[1], g141_cwt[2], fwhm_est_pix)
