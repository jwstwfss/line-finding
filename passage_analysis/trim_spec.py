from passage_analysis import *

def initialize_arrays(data, bluecut, redcut):
    """Define and trim initial arrays from a table of input data.

    Define the wavelength, flux, flux error, contamination and zeroth order 
    flags for an object's spectrum in one of the grisms. The first and 
    last 3 pixels of the spectrum are removed to avoid the mess caused when 
    spectra fall off the edge of the detector. This step is usually 
    redundant, but the trimming and masking steps may not catch these 
    problems at the edges.

        ### mask edges in wavelength where stuff gets crazy (i.e. low throughput)
    The arrays are then trimmed in wavelength to remove the areas where 
    stuff gets crazy (i.e. low throughput). The elements corresponding to
    wavelengths outside of range defined by (bluecut,redcut) are removed.    

    Args:
        data (astropy.io.ascii table): table of data from the 1-D dat file
        bluecut (float): the minimum wavelength 
        redcut (float): the maximum wavelength

    Returns:
        (tuple): tuple containing:
            lam (float): wavelength array
            flux (float): flux array
            error (float): flux error array
            contam (float): contamination array
            zeros (int): zeroth order flags
    """
    lam = data['lambda'][3:-3]
    flux = data['flux'][3:-3]
    error = data['ferror'][3:-3]
    contam = data['contam'][3:-3]
    zeros = data['zero'][3:-3]

    cut = (lam > bluecut) & (lam < redcut)
    lam = lam[cut]
    flux = flux[cut]
    error = error[cut]
    contam = contam[cut]
    zeros = zeros[cut]

    return lam,flux,error,contam,zeros


def newmask(array, maskedarray):
    newarray = np.ma.masked_where(np.ma.getmask(maskedarray), array)
    return newarray

def check_grism(tbdata, min_wav, max_wav):
    if tbdata is not None:
        _bluecut = min_wav
        _redcut = max_wav
        _specdata = initialize_arrays(tbdata, _bluecut, _redcut)
        
        _lam_spec = _specdata[0] 
        _flux_spec = _specdata[1]
        _error_spec = _specdata[2]
        _contam_spec = _specdata[3]
        _zero_spec = _specdata[4]

        ### only fit finite data
        _flux_spec = np.ma.masked_where(np.logical_or(~(np.isfinite(_flux_spec)),~(np.isfinite(_error_spec))), _flux_spec)
        _error_spec = newmask(_error_spec, _flux_spec)
        _contam_spec = newmask(_contam_spec, _flux_spec)

    else: _lam_spec=[]; _flux_spec=[]; _error_spec=[]; _contam_spec=[]; _zero_spec=[]

    return _lam_spec, _flux_spec, _error_spec, _contam_spec, _zero_spec



def trim_spec(tbdata_blue, tbdata_green, tbdata_red, config_pars, mask_zeros=False, return_masks=True):
    """Create one array of spectra, etc. from multiple grism 1d files
    
    If masking: 
        don't need to mask wavelenght or zeroth orders, keep these for 
        writing out file

    create an array keeping track of the masked regions

    Args:
        tbdata_blue ():  blue filter
        tbdata_red ():   red filter
        tbdata_green (): green/middle filter
        config_pars ():
        return_masks (Optional[bool]):

    Returns:
    """
    ### check each grism separately
    ### only fit finite data
    lam_spec_blue, flux_spec_blue, error_spec_blue, contam_spec_blue, zero_spec_blue = check_grism(tbdata_blue, config_pars['lambda_min'], config_pars['transition_wave1'])
    lam_spec_green, flux_spec_green, error_spec_green, contam_spec_green, zero_spec_green = check_grism(tbdata_green, config_pars['transition_wave1'], config_pars['transition_wave2'])
    lam_spec_red, flux_spec_red, error_spec_red, contam_spec_red, zero_spec_red = check_grism(tbdata_red, config_pars['transition_wave2'], config_pars['lambda_max'])

    # If one filter doesn't exist, then arrays are already empty. Just add them. 
    lam_spec = np.array(list(lam_spec_blue) + list(lam_spec_green) + list(lam_spec_red))
    flux_spec = np.array(list(flux_spec_blue) + list(flux_spec_green) + list(flux_spec_red))
    error_spec = np.array(list(error_spec_blue) + list(error_spec_green) + list(error_spec_red))
    contam_spec = np.array(list(contam_spec_blue) + list(contam_spec_green) + list(contam_spec_red))
    zero_spec = np.array(list(zero_spec_blue) + list(zero_spec_green) + list(zero_spec_red))
    # The above does what the code below does. I didn't want to go through
    # all possible combinations of which filter might be none the way the
    # code originally did. 
    
    # ### concatenate.
    # if tbdata_blue is None: 
    #     lam_spec = lam_spec_red 
    #     flux_spec = flux_spec_red
    #     error_spec = error_spec_red
    #     contam_spec = contam_spec_red 
    #     zero_spec = zero_spec_red
    # if tbdata_red is None: 
    #     lam_spec = lam_spec_blue
    #     flux_spec = flux_spec_blue 
    #     error_spec = error_spec_blue 
    #     contam_spec = contam_spec_blue 
    #     zero_spec = zero_spec_blue
    
    # if (tbdata_red is not None) & (tbdata_blue is not None): 
    #     lam_spec = np.append(lam_spec_blue, lam_spec_red) 
    #     flux_spec = np.ma.append(flux_spec_blue, flux_spec_red) 
    #     error_spec = np.ma.append(error_spec_blue, error_spec_red) 
    #     contam_spec = np.ma.append(contam_spec_blue, contam_spec_red) 
    #     zero_spec = np.append(zero_spec_blue, zero_spec_red) 


    if mask_zeros:
        ### remove bad 0th orders 
        flux_spec = np.ma.masked_where(zero_spec == 3, flux_spec)
        error_spec = newmask(error_spec, flux_spec)
        contam_spec = newmask(contam_spec, flux_spec)

    ### create an array of flags for masked regions
    masked_regions = np.zeros(lam_spec.shape, dtype=int)
    ### removed masked regions
    for k,v in config_pars.items():
        if 'mask_region' in k:
            # print(k, v) ; this prints the masked regions for testing (KVN, Jan2025)
            bluecut = v[0]
            redcut = v[1]
            mask = np.logical_and(lam_spec > bluecut, lam_spec < redcut)
            flux_spec = np.ma.masked_where(mask, flux_spec)
            error_spec = newmask(error_spec, flux_spec)
            contam_spec = newmask(contam_spec, flux_spec)
            masked_regions[mask] = 1
    
    if return_masks:
        # wavelength and zeroth orders are not masked, 
        # this is fine for measure_z_interactive
        return [lam_spec, flux_spec, error_spec, contam_spec, zero_spec, masked_regions] 
    else:
        # need to apply mask to wavelength and zeroth orders so all the 
        # arrays have the same shape
        # this is necessary for the cwt code
        lam_spec = newmask(lam_spec, flux_spec)
        zero_spec = newmask(zero_spec, flux_spec)
        # compress the arrays 
        outlam = np.ma.compressed(lam_spec)
        outflux = np.ma.compressed(flux_spec)
        outerror = np.ma.compressed(error_spec)
        outcontam = np.ma.compressed(contam_spec)
        outzero = np.ma.compressed(zero_spec)
        return outlam,outflux,outerror,outcontam,outzero


