#!/usr/bin/env python
##########################################################################
##########################################################################
# get_z_interactive.py By Nathaniel R. Ross, UCLA, nross@astro.ucla.edu
# usage: python get_z_interactive.py line_list
# Reads in list of emission lines from the WISP survey HST/WFC3/IR Grisms and
# plots the spectra, iterates through the detected emission lines, allowing the
# user to make line identifications or reject spurious lines quickly and
# efficiently.
#
# Version 1.0 updates the program to look for wavelength-calibrated 2d grism
# stamps first. Also, will now look for default line list name and make
# linelistfile an optional parameter. Furthermore, I have added the option
# to save the incomplete line list to file as you progress, but default is NOT to
# do this.
# Optional arguments for go() are 1. linelistfile (String) - path to line list file
#                                 2. save_temp (Bool) - Save progress in
#                                    linelist/Par???lines_with_redshifts.incomplete.dat
#                                 3. recover_temp (Bool) - Recover progress from previous
#                                    session using the .incomplete.dat files
#
# ** Major change in version 1.0:
# We now use xpa instead of IRAF to display the 2d grism stamps and full 2d
# direct images (instead of cutout stamp). Reason for this change was the desire
# to show the grism stamps that have a wavelength solution applied. When using
# the IRAF display command, ds9 would fail to recognize this coordinate system.
# The XPA system can be downloaded here: http://hea-www.harvard.edu/RD/xpa/index.html
#
##########################################################################
import os
import shutil
from glob import glob
import tarfile
import distutils
import numpy as np  # MDR 2022/05/17
import fileinput
import scipy
#### Added by KVN: this changes the matplotlib backend
import matplotlib
matplotlib.use('TkAgg')
####
import matplotlib.pylab as plt
from scipy.interpolate import splrep as spline
from astropy.table import Table
from distutils.sysconfig import *  # question-- what is this for?
import sys
from matplotlib import gridspec
import matplotlib.transforms as mtransforms
import matplotlib.patheffects as PathEffects
from copy import copy

# Explicitly import readline to make the text entry process easier on OSX
import readline

# SQLLite database support for data persistence
from passage_analysis.WISPLFDatabaseManager import WISPLFDatabaseManager as WDBM
from passage_analysis import *
import passage_analysis as passage
from passage_analysis.guis import showSpec2D_PASSAGE, showDirect_PASSAGE, panDirect_PASSAGE, panDispersed_PASSAGE 

verbose = True  # MDR 2022/05/17

if verbose == True:
    print("\nCompiling measure_z_interactive...\n")  # MDR 2022/05/17
    print("The current working director is...\n")  # MDR 2022/05/17
    print(os.getcwd())  # MDR 2022/05/17

#######################################################################

# Define the emission line vacuum wavelengths.
la_1216_vac = 1215.670
n5_1238_vac = 1238.821
n5_1242_vac = 1242.804
c4_1548_vac = 1548.203
c4_1550_vac = 1550.777
h2_1640_vac = 1640.420
o3_1660_vac = 1660.8092
o3_1666_vac = 1666.1497
s3_1883_vac = 1882.707
s3_1892_vac = 1892.030
c3_1907_vac = 1906.680
c3_1909_vac = 1908.730
m2_2796_vac = 2796.352
m2_2803_vac = 2803.531
o2_3727_vac = 3727.092
o2_3730_vac = 3729.875
hg_4342_vac = 4341.684
o3_4363_vac = 4364.436
h2_4686_vac = 4687.020
hb_4863_vac = 4862.683
o3_4959_vac = 4960.295
o3_5007_vac = 5008.240
o1_6300_vac = 6302.046
o1_6363_vac = 6365.536
n2_6550_vac = 6549.850
ha_6565_vac = 6564.610
n2_6585_vac = 6585.280
s2_6716_vac = 6718.290
s2_6731_vac = 6732.670
s3_9069_vac = 9071.100
s3_9532_vac = 9533.200
he10830_vac = 10832.86
pg_10941_vac = 10941.1
pb_12822_vac = 12821.6
pa_18756_vac = 18756.1

# Make all catalog header and data write commands loops over the 'flux_strings' variable.

"""
Note that the code expects that the strings in 'flux_strings' to match those in
'fitresults' with +'_flux' appended. This way the results can be looped over.
All emission line labels and parameters are defined at the top of the code,
except for the contamination flags, which must be defined within a function so
they reset for each object.
"""

supported_lines = [
    la_1216_vac,
    n5_1238_vac,
    c4_1548_vac,
    h2_1640_vac,
    o3_1660_vac,
    s3_1883_vac,
    c3_1907_vac,
    m2_2796_vac,
    o2_3727_vac,
    hg_4342_vac,
    o3_4363_vac,
    h2_4686_vac,
    hb_4863_vac,
    o3_5007_vac,
    o1_6300_vac,
    ha_6565_vac,
    s2_6716_vac,
    s3_9069_vac,
    s3_9532_vac,
    he10830_vac,
    pg_10941_vac,
    pb_12822_vac,
    pa_18756_vac,
]

# These lines are close to their doublets so are not plotted in ax1.
# However, their fluxes and fitresults are still saved to the catalog.
supported_lines_extra = [
    n5_1242_vac,
    c4_1550_vac,
    o3_1666_vac,
    s3_1892_vac,
    c3_1909_vac,
    m2_2803_vac,
    o2_3730_vac,
    o3_4959_vac,
    o1_6363_vac,
    n2_6550_vac,
    n2_6585_vac,
    s2_6731_vac,
]

supported_lines_strings = [
    r"Ly$\alpha$",
    "N V",
    "C IV",
    "He II",
    "O III]",
    "Si III]",
    "C III]",
    "Mg II",
    "[O II]",
    r"H$\gamma$",
    "[O III]",
    "He II",
    r"H$\beta$",
    "[O III]",
    "[O I]",
    r"H$\alpha$",
    "[S II]",
    "[S III]",
    "[S III]",
    "He I",
    r"P$\gamma$",
    r"P$\beta$",
    r"P$\alpha$",
]

flux_strings_1gauss = [
    "la_1216_wing",
    "n5_1238_1242",
    "c4_1548_1550",
    "h2_1640",
    "o3_1660_1666",
    "s3_1883_1892",
    "c3_1907_1909",
    "m2_2796_2803",
    "o2_3727_3730",
    "hg_4342",
    "o3_4363",
    "h2_4686",
    "hb_4863",
    "o3_4959_5007",
    "o1_6300_6363",
    "ha_6550_6565_6585",
    "s2_6716_6731",
    "s3_9069",
    "s3_9532",
    "he10830",
    "pg_10941",
    "pb_12822",
    "pa_18756",
]

flux_strings_2gauss = [
    "la_1216_wing",
    "n5_1238_1242",
    "c4_1548_1550",
    "h2_1640tot", "h2_1640nar", "h2_1640bro", 
    "o3_1660_1666",
    "s3_1883_1892",
    "c3_1907_1909",
    "m2_2796_2803",
    "o2_3727_3730",
    "hg_4342tot", "hg_4342nar", "hg_4342bro",
    "o3_4363tot", "o3_4363nar", "o3_4363bro",
    "h2_4686tot", "h2_4686nar", "h2_4686bro",
    "hb_4863tot", "hb_4863nar", "hb_4863bro",
    "o3_4959_5007",
    "o1_6300_6363",
    "ha_6550_6565_6585",
    "s2_6716_6731",
    "s3_9069",
    "s3_9532",
    "he10830tot", "he10830nar", "he10830bro",
    "pg_10941tot", "pg_10941nar", "pg_10941bro",
    "pb_12822tot", "pb_12822nar", "pb_12822bro",
    "pa_18756tot", "pa_18756nar", "pa_18756bro",
]
# , 'pg_10941','pb_12822','pa_18756']

contam_flags_string = "la_1216, n5_1238, n5_1242, c4_1548, c4_1550, h2_1640, \
o3_1660, o3_1666, s3_1883, s3_1892, c3_1907, c3_1909, \
m2_2796, m2_2803, o2_3727, o2_3730, hg_4342, o3_4363, \
h2_4686, hb_4863, o3_4959, o3_5007, o1_6300, o1_6363, \
n2_6550, ha_6565, n2_6585, s2_6716, s2_6731, s3_9069, \
s3_9532, he10830, pg_10941, pb_12822, pa_18756, cont"

# colors to help split up terminal output
# helpmsg = light blue
# obj = purple/magenta
# heading = green
# accept = green
# interim = cyan
# endc = clear color
setcolors = {
    "helpmsg": "\033[94m",
    "obj": "\033[95m",
    "heading": "\033[92m",
    "accept": "\033[92m",
    "working": "\033[0m",
    "interim": "\033[36m",
    "random": "\033[31m",
    "endc": "\033[0m"}
#######################################################################


def isFloat(string):
    try:
        float(string)
        return True
    except ValueError:
        return False


def getzeroorders(zeroorderpath, g="G141", magcut=23.5):  # MB: changed from 23.0
    if verbose == True:
        print("Running getzeroorders...\n")  # MDR 2022/05/17
    """
    Changing to return a table
    """
    zop = open(zeroorderpath, "r")
    zox = []
    zoy = []
    zora = []
    zodec = []
    zoid = []
    zmag = []
    for line in zop:
        if len(line) > 60:
            linesplit = line.split()
            zox.append(float(linesplit[1][0:-1]))
            zoy.append(float(linesplit[2][0:-3]))
            zox.append(float(linesplit[16][0:-3]))
            zoy.append(float(linesplit[17][0:-3]))
            zoid.append(int(linesplit[-2].split("{")[-1]))
            zmag.append(float(linesplit[-1][1:-2]))  # get mag from reg file

    zop.close()
    zoid = np.array(zoid)
    zoy = np.array(zoy)
    zox = np.array(zox)
    zora = np.array(zora)
    zodec = np.array(zodec)
    zmag = np.array(zmag)
    cond = zmag <= magcut
    t = Table(
        [zox[cond], zoy[cond], zoid[cond], zora[cond], zodec[cond]],
        names=("x", "y", "objid"),
    )
    return t


def getzeroorders_from_cat(zeroorderpath, g="G141", magcut=23.5):  
    # KN: updated function to use phot catalog
    #if verbose == True:
        # print("Running getzeroorders_from_cat...\n")  # MDR 2022/05/17
    """
    Changing to return a table
    """
    zeroorder = Table.read(zeroorderpath)
    t = Table(
        [
            zeroorder["x"],
            zeroorder["y"],
            zeroorder["id"],
            zeroorder["ra"],
            zeroorder["dec"],
            zeroorder["mag_auto"],
        ],
        names=("x", "y", "objid", "ra", "dec", "mag_auto"),
    )
    return t


def getfirstorders(firstorderpath):
    # if verbose == True:
    #     print("Running getfirstorders...\n")  # MDR 2022/05/17
    """
    Changing to return a table
    """
    fop = open(firstorderpath, "r")
    fox = []
    foy = []
    folen = []
    fowid = []
    #    foid=[] nothing is done with 1st order IDs anymore
    for line in fop:
        # if line[0]!='#':
        linesplit = line.split()
        fox.append(float(linesplit[1][0:-1]))  # [0:-1] strips off the comma.
        foy.append(float(linesplit[2][0:-1]))
        folen.append(float(linesplit[3][0:-1]))
        # python is weird.
        fowid.append(float(linesplit[-1].split("{")[-1].split("}")[0]))
    print(foid)
    t = Table([fox, foy, folen, fowid], names=("x", "y", "len", "width", "objid"))
    return t


def get_remaining_objects(full_obj_list, objid_done):
    if verbose == True:
        print("Running get_remaining_objects...")  # MDR 2022/05/17
    # array of flags
    wdone = np.in1d(np.array(full_obj_list), objid_done)
    remaining = np.copy(full_obj_list)
    mask = np.ones(full_obj_list.shape, dtype=bool)
    mask[wdone] = False
    remaining = remaining[mask]
    return remaining



def print_help_message():
    """
    Just putting this here to keep it out of the way.
    -- modified by FH 11/18/24
    """
    msg = setcolors["helpmsg"] + "Available Options:\n"
    msg += setcolors["heading"] + "\tOBJECT SPECIFIC OPTIONS:\n"
    msg += (
        setcolors["helpmsg"] + "\ta = accept object fit\n"
        "\tac = accept object fit, noting contamination\n"
        "\tr = reject object\n"
        "\tc = add comment\n"
        "\tuser = toggle between previously saved fits\n"
        "\tcontam = specify contamination to line flux and/or continuum\n"
        "\treset = reset interactive options back to default for this object\n"
        "\ts = print the (in progress) object summary\n\n"
    )
    # Look at the full line list and select which 
    msg += setcolors["heading"] + "\tEMISSION LINE SPECIFIC OPTIONS:\n"
    msg += (
        setcolors["helpmsg"] + "\tz = enter a different z guess\n"
        "\tw = enter a different emission line wavelength guess\n"
        "\tdz = change the allowable redshift difference between lines\n"
        "\tn = skip to next brightest line found in this object\n"
        "\tp = move to previous brightest line found in this object\n"
        "\t2gauss = double gaussian profile for the line being fitted\n"
        "\t1gauss = option to go back to 1 gaussian fit after selecting 2 gaussian fit\n"
        "\tha, hb, hg, o31, o32, o2, s2, s31, s32, lya, c4, pb, pa, pg = change strongest emission line\n"
        "\tPlease see the README file in the top-level directory for the full list of lines and their corresponding commands\n\n"

    )
    msg += setcolors["heading"] + "\tSPECTRUM SPECIFIC OPTIONS:\n"
    msg += (
        setcolors["helpmsg"] + "\tfw = change the fwhm guess in pixels\n"
        "\tt1, t2 = change transition wavelength between F115W and F150W (t1) and F150W and F200W (t2)\n"
        "\tm1, m2, m3, to m8 = mask up to eight discontinuous wavelength regions\n"
        "\tnodes = change the wavelengths for the continuum spline nodes\n"
        "\taddnodes = add wavelengths for the continuum spline nodes\n"
        "\trmnodes = remove wavelengths from the continuum spline nodes\n"
        "\tshiftallnodes = SHIFT ALL nodes used for the continuum spline by some wavelength \n"
        "\tbluecut = change the blue cutoff of the F115W grism\n"
        "\tredcut  = change the red cutoff of the F200W grism\n"
        "\tlincont = fit continuum as a line\n"
        "\tpolycont = fit continuum as a higher-order polynomial\n"
        "\tsplinecont = fit continuum as a spline (piecewise) polynomial\n"
        "\tgrismr = use only Grism-R spectrum for line-fitting\n"
        "\tgrismrcontam = use only Grism-R spectrum (with contamination) for line-fitting\n"
        "\tgrismc = use only Grism-C spectrum for line-fitting\n"
        "\tgrismccontam = use only Grism-C spectrum (with contamination) for line-fitting\n"
        "\tcomb = Use combined spectrum (default)\n"
        "\tcombcontam = Use combined spectrum with contamination\n\n"
    )
    msg += setcolors["heading"] + "\tDS9 SPECIFIC OPTIONS:\n"
    msg += (
        setcolors["helpmsg"] + "\tlin = linear z-scale\n"
        "\tlog = logarithmic\n"
        "\tzs102 = z1,z2 comma-separated range for G102 zscale\n"
        "\tzs141 = z1,z2 comma-separated range for G141 zscale\n"
        "\tdc = recenter direct images\n"
        "\treload = reload direct images\n"
        "\tdr = reload direct image reg files\n\n"
    )
    msg += setcolors["heading"] + "\tSOFTWARE SPECIFIC OPTIONS: \n"
    msg += (
        setcolors["helpmsg"] + "\th = print this message\n"
        "\tq = quit\n" + setcolors["endc"]
    )

    print(msg)



def check_masked_lines(fitresults, snr_meas_array, spdata, flux_strings):
    if verbose == True:
        print("Running check_masked_lines...\n")  # MDR 2022/05/17
    spec_lam = spdata[0]

    z = fitresults["redshift"]
    fwhm = fitresults["fwhm_g141"]

    for i, (wave, line) in enumerate(zip(supported_lines, flux_strings)):
        waveobs = wave * (1.0 + z)
        w = (
            spdata[1].mask[
                (spec_lam >= (waveobs - fwhm / 2.0))
                & (spec_lam <= (waveobs + fwhm / 2.0))
            ]
        ).any()
        if w:
            for dtype in ["flux", "error", "ew_obs"]:
                fitresults["%s_%s" % (line, dtype)] = -1.0
            snr_meas_array[i] = 0.0

    return fitresults, snr_meas_array


def comment_out_obj(par, obj, catalogname):
    if verbose == True:
        print("Running comment_out_obj...\n")  # MDR 2022/05/17
    # if a row already exists for this object, comment it out
    # objstr = '{:<8d}'.format(par) + '{:<6d}'.format(obj)
    objstr = "{:<6d}".format(obj)
    if os.path.exists(catalogname):
        for line in fileinput.input(catalogname, inplace=True):
            if objstr in line:
                print(
                    "#%s" % line,
                )
            else:
                print(
                    "%s" % line,
                )


def print_prompt(prompt, prompt_type="obj"):
    print(setcolors[prompt_type] + prompt + setcolors["endc"])


def write_object_summary(par, obj, fitresults, snr_meas_array, contamflags, comp_fit, summary_type="accept"):
    if verbose == True:
        print("Running write_object_summary...\n")  # MDR 2022/05/17
    """ """
    # string names for output
    # ---------
    # Added by KVN 13-Aug-2024 because double gaussian fit has extra lines. 
    if comp_fit == True: flux_strings = flux_strings_2gauss
    elif comp_fit == False: flux_strings = flux_strings_1gauss
    # ---------

    linefluxes = np.array([fitresults["%s_flux" % fs] for fs in flux_strings])

    # initial message
    msg = setcolors[summary_type] + "#" * 72
    msg += "\n## Par{} Obj {}:\n##   Fit Redshift: z = {:.4f}\n".format(
        par, obj, fitresults["redshift"]
    )

    # lines with S/N > 3
    good_snr = np.where(snr_meas_array > 3)
    msg += "##   Lines fit with S/N > 3:\n"
    for gsnr in good_snr[0]:
        # msg += '##\t%s: Flux = %.3e    S/N = %.2f\n'%(supported_lines_strings[gsnr],
        #                         linefluxes[gsnr], snr_meas_array[gsnr])
        msg += "##\t%s: Flux = %.3e    S/N = %.2f\n" % (
            flux_strings[gsnr],
            linefluxes[gsnr],
            snr_meas_array[gsnr],
        )

    cfout = [
        "%s:%i" % (cf, contamflags[cf]) for cf in contamflags if contamflags[cf] > 0
    ]
    msg += "##   Contamination flags set:\n##\t" + ", ".join(cfout) + "\n"
    msg += "#" * 72 + setcolors["endc"]
    print(msg)


def make_tarfile(outdir, path_to_code):
    if verbose == True:
        print("Running make_tarfile...\n")  # MDR 2022/05/17
    """ """
    # 
    print(path_to_code+"/default.config")
    # copy default.config into output directory to keep a copy
    shutil.copy(path_to_code+"/default.config", outdir)
    with tarfile.open("%s.tar.gz" % outdir, "w:gz") as tar:
        tar.add(outdir, arcname=os.path.basename(outdir))


# DEFINE FUNCTION FOR ADDING EMISSION LINE LABELS TO ZOOM-IN PLOTS.
def add_line_labels(ax, axtrans, xmin, xmax):
    # if (verbose == True): print 'Calling add_line_labels()...\n'
    for li, lstring, sn_meas in zip(
        lamobs, supported_lines_strings, snr_meas_array):
        if (li > xmin + 1) & (li < xmax - 1):
            ax.axvline(x=li, color="b")
            stringplot = lstring + "  (" + str(round(sn_meas, 1)) + ") "
            # use data coordinates for x-axis and axes coords for y-axis
            ax.text(li, 0.92, stringplot, ha="right", fontsize="14", transform=axtrans)
                # ax.set_title(stringplot, fontsize='14')
    return "Added line label to", str(ax)

def add_extra_lines(ax, axtrans, xmin, xmax):
    # if (verbose == True): print 'Calling add_extra_lines()...\n'
    for li in lamobs_extra:
        if (li > xmin + 1) & (li < xmax - 1):
            ax.axvline(x=li, color="b")  # , linestyle='--')
            # stringplot = lstring + '  (' + str(round(sn_meas, 1)) + ') '
            # use data coordinates for x-axis and axes coords for y-axis
            # ax.text(li, 0.92, stringplot, ha='right', fontsize='14', transform=axtrans)
            # ax.set_title(stringplot, fontsize='14')
    return "Added line label to", str(ax)


def plot_chooseSpec(spdata1, spdata2, spdata3, config_pars, plottitle, outdir, zset=None, orientation=""):
    if verbose == True:
        print("\nRunning plot_chooseSpec...\n")  # adapted from plot_object by KVN 2024/07/22
    """
    # Written by KVN to allow the user to see spectra in all orientations
    # save the figure for everything
    # previous figures are overwritten
    """
    # define the filename that will be used for figures.
    plotfilename = os.path.join(outdir, "figs", "%s_fit_spectra.png" % (plottitle))

    # initialize the wavelength, flux, uncertainty, contamination, and zero order arrays.
    if spdata1 != None:
        spec_lam1 = spdata1[0]; spec_val1 = spdata1[1]; spec_unc1 = spdata1[2];
        spec_con1 = spdata1[3]; spec_zer1 = spdata1[4];
    else: spec_lam1=[]; spec_val1 =[]; spec_unc1=[]; spec_con1=[]; spec_zer1=[]
    
    if spdata2 != None:
        spec_lam2 = spdata2[0]; spec_val2 = spdata2[1]; spec_unc2 = spdata2[2]; 
        spec_con2 = spdata2[3]; spec_zer2 = spdata2[4]; 
    else: spec_lam2=[]; spec_val2 =[]; spec_unc2=[]; spec_con2=[]; spec_zer2=[]

    if spdata3 != None:
        spec_lam3 = spdata3[0]; spec_val3 = spdata3[1]; spec_unc3 = spdata3[2]; 
        spec_con3 = spdata3[3]; spec_zer3 = spdata3[4];
    else: spec_lam3=[]; spec_val3 =[]; spec_unc3=[]; spec_con3=[]; spec_zer3=[]

    # apply the mask to the wavelength array
    masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val1), spec_lam1)

    # generate the plot grid.
    plt.ion()
    fig = plt.figure(1, figsize=(11, 12), dpi=75)
    plt.clf()
    # gs = gridspec.GridSpec(3, 4)
    gs = gridspec.GridSpec(3, 1)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    xmin = np.ma.min(spec_lam1) - 10.0  # 200.0 - M.D.R. - 10/22/2020
    xmax = np.ma.max(spec_lam1) + 10.0  # 200.0 - M.D.R. - 10/22/2020
    ymin = -0.2 * np.ma.max(spec_val1)
    ymax = 1.2 * np.ma.max(spec_val1)

    # the line widths for the data and overlaid fit.
    lw_data = 2.0
    lw_fits = 1.75

    ax1.plot(spec_lam1, spec_val1, "k", spec_lam1, spec_con1, "hotpink", spec_lam1, spec_val1+spec_con1, "cornflowerblue", drawstyle="steps-mid", lw=lw_data)
    ax2.plot(spec_lam2, spec_val2, "k", spec_lam2, spec_con2, "hotpink", spec_lam2, spec_val2+spec_con2, "cornflowerblue", drawstyle="steps-mid", lw=lw_data)
    ax3.plot(spec_lam3, spec_val3, "k", spec_lam3, spec_con3, "hotpink", spec_lam3, spec_val3+spec_con3, "cornflowerblue", drawstyle="steps-mid", lw=lw_data)
    # ax1.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)


    ax2.plot(np.linspace(0,1), np.linspace(0,1), 'cornflowerblue', label = 'With Contamination')
    ax2.plot(np.linspace(0,1), np.linspace(0,1), 'k', label = 'Without Contamination')
    ax2.plot(np.linspace(0,1), np.linspace(0,1), 'hotpink', label = 'Contamination')
    ax2.legend(loc=1)
    # plot any masked regions
    for mr in ["mask_region1", "mask_region2", "mask_region3", "mask_region4", "mask_region5", "mask_region6", "mask_region7", "mask_region8"]:
        if (config_pars[mr][0] != 0.0) & (config_pars[mr][1] != 0.0):
            for ax in [ax1, ax2, ax3]:  # [ax1, ax2]:
                trans = mtransforms.blended_transform_factory(
                    ax.transData, ax.transAxes
                )
                handles, labels = ax.get_legend_handles_labels()
                if "masked regions" in labels:
                    maskedlabel = None
                else:
                    maskedlabel = "masked regions"
                ax.fill_between(config_pars[mr], 0, 1, color="grey", alpha=0.3, transform=trans, label=maskedlabel)
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 0:
        ax1.legend(bbox_to_anchor=[1.0, 1.08], loc="upper right")

    for ax in [ax1, ax2, ax3]:
        ax.set_xlabel("Observed Wavelength ($\AA$)", size="xx-large")
        # ax1.set_ylabel(r'F$_\lambda$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$', size='xx-large')
        ax.set_ylabel("Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", size="xx-large")
        ax.set_xlim([xmin, xmax])
        ax.set_ylim([ymin, ymax])
    ax1.set_title(plottitle+" Combined", size="xx-large")
    ax2.set_title(plottitle+" R", size="xx-large")
    ax3.set_title(plottitle+" C", size="xx-large")

    plt.tight_layout() 
    fig.savefig(plotfilename)
    plt.draw()
    

def plot_object(zguess, zfit, spdata, config_pars, snr_meas_array, snr_tot_others, full_fitmodel, full_contmodel, broad_fitmodel, current_lam, lamlines_found, index_of_strongest_line, contmodel, plottitle, outdir, zset=None, orientation=""):
    if verbose == True:
        print("\nRunning plot_object...\n")  # MDR 2022/05/17
    """
    # save the figure for everything, junk objects and all
    # previous figures are overwritten
    """
    # the expected wavelengths of emission lines given the zguess
    lamobs = (1 + zguess) * np.array(supported_lines)
    lamobs_extra = (1 + zguess) * np.array(supported_lines_extra)

    # define the filename that will be used for figures.
    plotfilename = os.path.join(outdir, "figs", "%s_fit%s.png" % (plottitle, orientation))

    # initialize the wavelength, flux, uncertainty, contamination, and zero order arrays.
    spec_lam = spdata[0]
    spec_val = spdata[1]
    spec_unc = spdata[2]
    spec_con = spdata[3]
    spec_zer = spdata[4]

    # apply the mask to the wavelength array
    masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)

    # limit the S/N array to MUSE so we can zoom-in on those lines.
    snr_muse = snr_meas_array[np.argwhere(lamobs < config_pars["mask_region1"][0])]

    # replace the occasional NaN value with zero to prevent gridspec crashes.
    snr_muse = np.nan_to_num(snr_muse)

    # determine the number of emission lines with S/N > 3.
    num_gridspec_plots = len(np.where(snr_muse >= 100.0)[0])

    # the maximum number of zoom-in plots is 7 for visibility.
    if num_gridspec_plots > 7:
        num_gridspec_plots = 7

    # if (verbose == True): print 'There are '+str(num_gridspec_plots)+' lines with S/N > 3 in MUSE.\n'

    # generate the plot grid.
    plt.ion()
    fig = plt.figure(2, figsize=(11, 8), dpi=75)
    plt.clf()
    # gs = gridspec.GridSpec(3, 4)
    gs = gridspec.GridSpec(3, num_gridspec_plots + 1)
    ax1 = fig.add_subplot(gs[0:2, :])
    ax2 = fig.add_subplot(gs[2:, 0:1])

    if num_gridspec_plots == 1:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        all_ax = [ax1, ax2, ax3]
    elif num_gridspec_plots == 2:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        all_ax = [ax1, ax2, ax3, ax4]
    elif num_gridspec_plots == 3:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        all_ax = [ax1, ax2, ax3, ax4, ax5]
    elif num_gridspec_plots == 4:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6]
    elif num_gridspec_plots == 5:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        ax7 = fig.add_subplot(gs[2:, 5:6])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7]
    elif num_gridspec_plots == 6:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        ax7 = fig.add_subplot(gs[2:, 5:6])
        ax8 = fig.add_subplot(gs[2:, 6:7])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8]
    elif num_gridspec_plots >= 7:
        ax3 = fig.add_subplot(gs[2:, 1:2])
        ax4 = fig.add_subplot(gs[2:, 2:3])
        ax5 = fig.add_subplot(gs[2:, 3:4])
        ax6 = fig.add_subplot(gs[2:, 4:5])
        ax7 = fig.add_subplot(gs[2:, 5:6])
        ax8 = fig.add_subplot(gs[2:, 6:7])
        ax9 = fig.add_subplot(gs[2:, 7:8])
        all_ax = [ax1, ax2, ax3, ax4, ax5, ax6, ax7, ax8, ax9]
    else:
        all_ax = [ax1, ax2]

    xmin = np.ma.min(spec_lam) - 10.0  # 200.0 - M.D.R. - 10/22/2020
    xmax = np.ma.max(spec_lam) + 10.0  # 200.0 - M.D.R. - 10/22/2020
    ymin = -0.2 * np.ma.max(spec_val)
    ymax = 1.2 * np.ma.max(spec_val)

    # the line widths for the data and overlaid fit.
    lw_data = 2.0
    lw_fits = 1.75

    ax1.plot(spec_lam, spec_val, "k", spec_lam, spec_con, "hotpink", drawstyle="steps-mid", lw=lw_data)
    # ax1.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)

    # transforms for plotting in data and axes coordinates
    ax1trans = mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes)
    ax2trans = mtransforms.blended_transform_factory(ax2.transData, ax2.transAxes)

    # contamination model
    ax1.fill_between(spec_lam, spec_con, -1, color="#ff69b4", alpha=0.1, step="pre")

    # plot observed wavelengths of all the possible lines.
    for li, lstring, sn_meas in zip(lamobs, supported_lines_strings, snr_meas_array):
        if (li > xmin + 100) & (li < xmax - 100):
            for ax in all_ax:  # [ax1, ax2, ax3, ax4]: # for ax in [ax1, ax2]: # - M.D.R. - 10/22/2020
                if ax != ax2:  # skip s/n plot.
                    ax.axvline(x=li, color="b")
                    # add blue lines for close doublets without text label.
                    # for li in lamobs_extra:
                    #     if (li > xmin + 1) & (li < xmax - 1):
                    #         ax.axvline(x=li, color='b', linestyle='--')
                    #
            stringplot = lstring  # + '  (' + str(round(sn_meas, 1)) + ') '

            # use data coordinates for x-axis and axes coords for y-axis
            ax1.text(li, 0.85, stringplot, rotation="vertical",
                ha="right", fontsize="16", transform=ax1trans)

    ax1.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
    ax1.plot(spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits)
    ax1.plot(spec_lam, broad_fitmodel, color="r", linestyle="-", lw=0.5)

    # plot 0th orders
    w = np.where(spec_zer == 3)
    spec_zero_bad = spec_zer * 0 - 1
    spec_zero_bad[w] = 1.0
    # mild zeroth orders
    w = np.where(spec_zer == 2)
    spec_zero_mild = spec_zer * 0 - 1
    spec_zero_mild[w] = 1.0
    for ax in [ax1, ax2]:
        # use data coordinates for x-axis and axes coords for y-axis
        trans = mtransforms.blended_transform_factory(ax.transData, ax.transAxes)
        if np.any(spec_zero_bad[spec_zero_bad != -1]):
            ax.fill_between(spec_lam, 0, 1, where=spec_zero_bad == 1,color="red", alpha=0.3, transform=trans, label="Major 0th order contam" )
        if np.any(spec_zero_mild[spec_zero_mild != -1]):
            ax.fill_between(spec_lam, 0, 1, where=spec_zero_mild == 1, color="orange", alpha=0.3, transform=trans, label="Minor 0th order contam")

    # plot any masked regions
    for mr in ["mask_region1", "mask_region2", "mask_region3", "mask_region4", "mask_region5", "mask_region6", "mask_region7", "mask_region8"]:
        if (config_pars[mr][0] != 0.0) & (config_pars[mr][1] != 0.0):
            for ax in all_ax:  # [ax1, ax2]:
                trans = mtransforms.blended_transform_factory(
                    ax.transData, ax.transAxes
                )
                handles, labels = ax.get_legend_handles_labels()
                if "masked regions" in labels:
                    maskedlabel = None
                else:
                    maskedlabel = "masked regions"
                ax.fill_between(config_pars[mr], 0, 1, color="grey", alpha=0.3, transform=trans, label=maskedlabel)
    handles, labels = ax.get_legend_handles_labels()
    if len(labels) > 0:
        ax1.legend(bbox_to_anchor=[1.0, 1.08], loc="upper right")

    # find values of spec_lam nearest to the nodes
    nodelam = config_pars["node_wave"]
    nl_arr = []
    cont_node = []
    for nl in nodelam:
        w = np.argmin(np.abs(spec_lam - nl))
        nl_arr.append(spec_lam[w])
        cont_node.append(full_contmodel[w])
    ax1.plot(nl_arr, cont_node, "ko", ms=9)

    # repeat for line_candidates
    lf_lam = []
    lf_cont = []
    for lf in lamlines_found:
        w = np.argmin(np.abs(spec_lam - lf))
        lf_lam.append(spec_lam[w])
        lf_cont.append(full_contmodel[w])
    ax1.plot(lf_lam, lf_cont, "bo", ms=9)

    #   indicate "current" line
    #   current_lam = lamlines_found[index_of_strongest_line]
    current_cont = contmodel[np.argmin(np.abs(np.ma.compressed(masked_spec_lam) - current_lam))]
    ax1.plot(current_lam, current_cont, "ro", ms=10)

    ax1.set_xlabel("Observed Wavelength ($\AA$)", size="xx-large")
    # ax1.set_ylabel(r'F$_\lambda$ erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$', size='xx-large')
    ax1.set_ylabel("Flux (erg s$^{-1}$ cm$^{-2}$ $\AA^{-1}$)", size="xx-large")
    ax1.set_xlim([xmin, xmax])
    ax1.set_ylim([ymin, ymax])
    ax1.set_title(plottitle, size="xx-large")

    # second panel for s/n
    s2n = (spec_val - full_contmodel) / spec_unc
    s2n_lam = spec_lam
    mask = np.logical_and(s2n > -10000.0, s2n < 10000.0)
    s2n = s2n[mask]
    s2n_lam = s2n_lam[mask]
    ax2.plot(s2n_lam, s2n, "k-", drawstyle="steps-mid", lw=lw_data)
    ymin = -0.2 * s2n.max()
    ymax = 1.2 * s2n.max()  # ymax = 1.5 * s2n.max() # M.D.R 01/07/2021
    ax2.axhline(y=config_pars["n_sigma_above_cont"], c="r")
    # for li in lamobs:
    #     ax2.axvline(x=li, color='b')
    # UNCOMMENT BELOW IF YOU WANT THE CYAN LINE INDICATING THE TRANSITION WAVE
    # ax2.axvline(x=config_pars['transition_wave'], c='c', linestyle=':', lw=3)
    # ax2.set_xlabel(r'$\lambda$ ($\AA$)', size='xx-large')
    # ax2.set_ylabel(r'S/N', size='xx-large')
    ax2.set_title(r"S/N", size="xx-large")
    ax2.set_xlim([xmin, xmax])
    ax2.set_ylim(ymin, ymax)
    # fig = plt.gcf() a

    # - M.D.R. - 10/22/2020
    # find max S/N in list of detected lines.
    # pull index location, determine lambda of that line.
    # set ymin/max based on centroid +/-i range based on dispersion.
    # print('snr_meas_array =', snr_meas_array)  # M.D.R 01/27/2021
    # print(np.max(snr_meas_array)) # M.D.R 01/27/2021
    # print(np.argmax(snr_meas_array)) # M.D.R 01/27/2021
    # lamobs_maxsnr = lamobs[np.argmax(snr_meas_array)]
    lamobs_blue = lamobs[np.argwhere(lamobs < config_pars["mask_region1"][0])]
    # lamobs_red = lamobs[np.argwhere(lamobs > config_pars['mask_region1'][1])]
    # array of the indices for lines with s/n > 3 from strongest to weakest.
    snr_muse_idx = np.flip(np.argsort(snr_muse, axis=0))

    # loop over lines with s/n > 3 and create zoom-in plots.
    for x in range(num_gridspec_plots):
        # print 'x =', x
        if snr_muse[snr_muse_idx[x]][0] >= 10.0:

            idx_i = snr_muse_idx[x][0]
            snr_i = snr_muse[snr_muse_idx[x]][0][0]
            lam_i = lamobs_blue[snr_muse_idx[x]][0][0]
            # print 'idx_i = ', idx_i
            # print 'snr_i = ', snr_i
            # print 'lam_i = ', lam_i

            # the wavelength array index corresponding to the detected line wavelength.
            snr_muse_idy = min(enumerate(spec_lam), key=lambda x: abs(lam_i - x[1]))[0]
            # half-width of zoom-in plots in angstroms
            half_xrange_ang = 1000
            # half-width of zoom-in plots in pixels
            half_xrange_pix = int((half_xrange_ang / config_pars["dispersion_blue"]))

            xmin = lam_i - half_xrange_ang
            xmax = lam_i + half_xrange_ang

            # update to a sigma rejected min?
            # print(spec_val)
            # print(np.shape(spec_val))
            # print(snr_muse_idy - half_xrange_pix, snr_muse_idy + half_xrange_pix)
            ymin = 0.9 * np.nanmin(spec_val[snr_muse_idy - half_xrange_pix : snr_muse_idy + half_xrange_pix])
            ymax = 1.2 * np.nanmax(spec_val[snr_muse_idy - half_xrange_pix : snr_muse_idy + half_xrange_pix])

            if x == 0:
                ax3trans = mtransforms.blended_transform_factory(ax3.transData, ax3.transAxes)
                ax3.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax3.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax3.plot(spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits)
                ax3.set_xlim([xmin, xmax])
                ax3.set_ylim([ymin, ymax])
            #                 add_line_labels(ax3, ax3trans, xmin, xmax)
            #                 add_extra_lines(ax3, ax3trans, xmin, xmax)
            elif x == 1:
                ax4trans = mtransforms.blended_transform_factory(ax4.transData, ax4.transAxes)
                ax4.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax4.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax4.plot(spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits)
                ax4.set_xlim([xmin, xmax])
                ax4.set_ylim([ymin, ymax])
            #                 add_line_labels(ax4, ax4trans, xmin, xmax)
            #                 add_extra_lines(ax4, ax4trans, xmin, xmax)
            elif x == 2:
                ax5trans = mtransforms.blended_transform_factory(
                    ax5.transData, ax5.transAxes
                )
                ax5.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax5.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax5.plot(
                    spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits
                )
                ax5.set_xlim([xmin, xmax])
                ax5.set_ylim([ymin, ymax])
            #                 add_line_labels(ax5, ax5trans, xmin, xmax)
            #                 add_extra_lines(ax5, ax5trans, xmin, xmax)
            elif x == 3:
                ax6trans = mtransforms.blended_transform_factory(
                    ax6.transData, ax6.transAxes
                )
                ax6.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax6.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax6.plot(
                    spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits
                )
                ax6.set_xlim([xmin, xmax])
                ax6.set_ylim([ymin, ymax])
            #                 add_line_labels(ax6, ax6trans, xmin, xmax)
            #                 add_extra_lines(ax6, ax6trans, xmin, xmax)
            elif x == 4:
                ax7trans = mtransforms.blended_transform_factory(
                    ax7.transData, ax7.transAxes
                )
                ax7.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax7.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax7.plot(
                    spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits
                )
                ax7.set_xlim([xmin, xmax])
                ax7.set_ylim([ymin, ymax])
            #                 add_line_labels(ax7, ax7trans, xmin, xmax)
            #                 add_extra_lines(ax7, ax7trans, xmin, xmax)
            elif x == 5:
                ax8trans = mtransforms.blended_transform_factory(ax8.transData, ax8.transAxes)
                ax8.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax8.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax8.plot(spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits)
                ax8.set_xlim([xmin, xmax])
                ax8.set_ylim([ymin, ymax])
            #                 add_line_labels(ax8, ax8trans, xmin, xmax)
            #                 add_extra_lines(ax8, ax8trans, xmin, xmax)
            elif x == 6:
                ax9trans = mtransforms.blended_transform_factory(ax9.transData, ax9.transAxes)
                ax9.plot(spec_lam, spec_val, "k", lw=lw_data)
                ax9.plot(spec_lam, full_fitmodel, color="r", lw=lw_fits)
                ax9.plot(spec_lam, full_contmodel, color="b", linestyle="--", lw=lw_fits)
                ax9.set_xlim([xmin, xmax])
                ax9.set_ylim([ymin, ymax])
            #                 add_line_labels(ax9, ax9trans, xmin, xmax)
            #                 add_extra_lines(ax9, ax9trans, xmin, xmax)
            else:
                print("Failed to set XY limits.")

    if zset is None:
        addtext = ("In progress (z = {:.4f}".format(zfit) +", "+"wSNR = {:.2f})".format(snr_tot_others))
        addtextcolor = "orange"
    elif zset == 0:
        addtext = "Rejected"
        addtextcolor = "red"
    elif zset == 1:
        addtext = ("Accepted (z = {:.4f}".format(zfit) + ", " + "wSNR = {:.2f})".format(snr_tot_others))
        addtextcolor = "green"

    fig.text(0.3, 0.96, addtext, ha="right", va="bottom", color=addtextcolor,
        fontsize=18, fontweight=500, path_effects=[PathEffects.withStroke(linewidth=0.5,foreground="k")])

    plt.tight_layout()  # MDR 2022/05/19 - Added to fill screen when fitting.
    fig.savefig(plotfilename)
    plt.draw()

def inspect_object(
    user,
    par,
    obj,
    objinfo,
    lamlines_found,
    ston_found,
    g102zeros,
    g141zeros,
    linelistoutfile,
    commentsfile,
    remaining,
    allobjects,
    show_dispersed=True,
    stored_fits=False,
    path_to_data=" ",
    path_to_output=" ",
    path_to_code=" ",
    orientation='combined',
    comp_fit = False,
    polycont_fit = False,
    lincont_fit = False):
    
    if verbose == True: print("Running inspect_object...\n")  # MDR 2022/05/17
    """
    An attempt to move all object-specific tasks
    """
    # set up and filenames
    outdir = "Par%s_output_%s" % (par, user)

    if path_to_data == " ":
        specnameg1 = ("Par%i_" + str(obj).zfill(5) + ".specG115_1D.dat" % (par))
        specnameg2 = ("Par%i_" + str(obj).zfill(5) + ".specG150_1D.dat" % (par))
        specnameg3 = ("Par%i_" + str(obj).zfill(5) + ".specG200_1D.dat" % (par))
        specnameg1_R = ("Par%i_" + str(obj).zfill(5) + ".specG115_1D_R.dat" % (par))
        specnameg2_R = ("Par%i_" + str(obj).zfill(5) + ".specG150_1D_R.dat" % (par))
        specnameg3_R = ("Par%i_" + str(obj).zfill(5) + ".specG200_1D_R.dat" % (par))
        specnameg1_C = ("Par%i_" + str(obj).zfill(5) + ".specG115_1D_C.dat" % (par))
        specnameg2_C = ("Par%i_" + str(obj).zfill(5) + ".specG150_1D_C.dat" % (par))
        specnameg3_C = ("Par%i_" + str(obj).zfill(5) + ".specG200_1D_C.dat" % (par))

    else:
        base_path = path_to_data+ "Par"+ str(par)+ "/Spectra/Par"+ str(par)+ "_" + str(obj).zfill(5)
        specnameg1 = (base_path + ".specG115_1D.dat")
        specnameg2 = (base_path + ".specG150_1D.dat")
        specnameg3 = (base_path + ".specG200_1D.dat")
        specnameg1_R = (base_path + ".specG115_1D_R.dat")
        specnameg2_R = (base_path + ".specG150_1D_R.dat")
        specnameg3_R = (base_path + ".specG200_1D_R.dat")
        specnameg1_C = (base_path + ".specG115_1D_C.dat")
        specnameg2_C = (base_path + ".specG150_1D_C.dat")
        specnameg3_C = (base_path + ".specG200_1D_C.dat")

    plottitle = "PASSAGE_%i" % (obj)
    fitdatafilename = os.path.join(outdir, "fitdata/%s_fitspec" % plottitle)
    availgrism = ""
    # read in 1D spectrum
    if os.path.exists(specnameg1):
        availgrism += "g115"
        tab_blue = Table.read(specnameg1, format="ascii")
        tab_blue.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_blue_cont = np.copy(tab_blue)
        tab_blue_cont['flux'] = tab_blue['flux'] + tab_blue['contam']
    else: tab_blue = None; tab_blue_cont = None 
    if os.path.exists(specnameg2):
        availgrism += "g150"
        tab_mid = Table.read(specnameg2, format="ascii")
        tab_mid.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_mid_cont = np.copy(tab_mid)
        tab_mid_cont['flux'] = tab_mid['flux'] + tab_mid['contam']
    else: tab_mid = None; tab_mid_cont = None
    if os.path.exists(specnameg3):
        availgrism += "g200"
        tab_red = Table.read(specnameg3, format="ascii")
        tab_red.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_red_cont = np.copy(tab_red)
        tab_red_cont['flux'] = tab_red['flux'] + tab_red['contam']
    else: tab_red = None; tab_red_cont = None

    if availgrism == "g115g150g200": availgrism = "both"  

    if os.path.exists(specnameg1_R):
        tab_blue_R = Table.read(specnameg1_R, format="ascii")
        tab_blue_R.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_blue_R_cont = np.copy(tab_blue_R)
        tab_blue_R_cont['flux'] = tab_blue_R['flux'] + tab_blue_R['contam']
    else: tab_blue_R = None; tab_blue_R_cont = None
    if os.path.exists(specnameg2_R):
        tab_mid_R = Table.read(specnameg2_R, format="ascii")
        tab_mid_R.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_mid_R_cont = np.copy(tab_mid_R)
        tab_mid_R_cont['flux'] = tab_mid_R['flux'] + tab_mid_R['contam']
    else: tab_mid_R = None; tab_mid_R_cont= None
    if os.path.exists(specnameg3_R):
        tab_red_R = Table.read(specnameg3_R, format="ascii")
        tab_red_R.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_red_R_cont = np.copy(tab_red_R)
        tab_red_R_cont['flux'] = tab_red_R['flux'] + tab_red_R['contam']
    else: tab_red_R = None; tab_red_R_cont = None

    if os.path.exists(specnameg1_C):
        tab_blue_C = Table.read(specnameg1_C, format="ascii")
        tab_blue_C.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_blue_C_cont = np.copy(tab_blue_C)
        tab_blue_C_cont['flux'] = tab_blue_C['flux'] + tab_blue_C['contam']
    else: tab_blue_C = None; tab_blue_C_cont = None
    if os.path.exists(specnameg2_C):
        tab_mid_C = Table.read(specnameg2_C, format="ascii")
        tab_mid_C.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_mid_C_cont = np.copy(tab_mid_C)
        tab_mid_C_cont['flux'] = tab_mid_C['flux'] + tab_mid_C['contam']
    else: tab_mid_C = None; tab_mid_C_cont = None
    if os.path.exists(specnameg3_C):
        tab_red_C = Table.read(specnameg3_C, format="ascii")
        tab_red_C.rename_columns(["wave", "ferr", "zeroth"], ["lambda", "ferror", "zero"])
        tab_red_C_cont = np.copy(tab_red_C)
        tab_red_C_cont['flux'] = tab_red_C['flux'] + tab_red_C['contam']
    else: tab_red_C = None; tab_red_C_cont = None 

    # =================== Show spec2d new (begin) =====================

    showSpec2D_PASSAGE(par, obj, path_to_data=path_to_data)

    # =================== Show spec2d new (end) =====================

    # pan full images to the new object
    # showDirectNEW(obj, par, g141zeros, path_to_wisp_data=path_to_wisp_data)

    # tbaines: pan full frame image to object of interest (new)

    #     if show_dispersed:
    #         showDispersed(obj, par, path_to_wisp_data = path_to_wisp_data)

    # define parameters for this object
    ra = objinfo["ra"]
    dec = objinfo["dec"]
    a_image = objinfo["a_image"]
    b_image = objinfo["b_image"]
    jmag = objinfo["jmag"]
    jerr = objinfo["jerr"]
    hmag = objinfo["hmag"]
    herr = objinfo["herr"]

    # pan to object in direct frames
    panDirect_PASSAGE(ra, dec)
    #panDirect_PASSAGE(x_pix, y_pix)
    ### Updated by KVN because this panning needs to be offset for each grism
    ### See new panDispersed_PASSAGE function in guis.py
    panDispersed_PASSAGE(obj, parno=par, path_to_data=path_to_data)

    # start with a fresh set of config pars
    config_pars = read_config(path_to_code+"/default.config", availgrism=availgrism)
    print('Reading in this config file: ', path_to_code+"/default.config")

    # Data have been successfully loaded for this object. If it has been inspected
    # previously, the original results will have been stored in the SQLLite database
    # and a retrieval obtion should be offered.
    #    databaseManager = WDBM(dbFileNamePrefix='Par{}'.format(par))
    databaseManager = WDBM(
        dbFileNamePrefix=os.path.join(outdir, "Par{}".format(par)))  
    # databaseManager = WDBM(dbFileNamePrefix=os.path.join(outdir,'Par{}'.format(par))) - M.D.R. - 10/08/2020
    mostRecentObjectData = databaseManager.getMostRecentObject()
    # if mostRecentObjectData is not None:
    #     print('Most recent object in database: Par {}, Obj {}, Date {}'.format(*mostRecentObjectData))
    # else:
    #     print('Database is empty.')
    catalogueEntryData = databaseManager.loadCatalogueEntry(parNumber=par, objectId=obj)

    ### # MDR 2022/05/17
    print("\nSearching for previous fits to object {}...\n".format(str(obj)))
    # print('file = ', outdir + '/done_%s'%user)

    if os.path.exists(outdir + "/done_%s" % user):  # needed for the first run before file exists
        with open(outdir + "/done_%s" % user) as f:
            for index, line in enumerate(f):
                # print("Line {}: {}".format(index, line.strip()))
                if int(line.strip()) == obj:
                    catalogueEntryData = 1
                    print("Match found...\n")
                    break

        if catalogueEntryData != 1:
            print("No match found...\n")

    print("Done searching for previous fits...")
    ### # MDR 2022/05/17

    rejectPrevFit = True

    if catalogueEntryData == 1:
        # nonFitResults, fitResults = catalogueEntryData
        # (par_db, obj_db, ra_db, dec_db, jmag_db, hmag_db, a_image_db, b_image_db, contamflag_db, entrytime_db) = nonFitResults
        # print('Found previous fit results for Pointing {}, Object {}.\nEnter "y" to accept the earlier fit.'.format(par_db, obj_db))
        print("You have already fit Obj {}. Refit? [y/N]".format(str(obj)))
        rejectPrevFit = input("> ").strip().lower() == "y"
        if rejectPrevFit:
            comment_out_obj(par, obj, linelistoutfile)

    # print('Accepting previous fit.' if acceptPrevFit else 'Re-fitting this object.')
    # else :
    # print('No previous fit results found. Fitting this object now.')

    # get line, fwhm, z estimate
    # choose the lamline that has the highest S/N estimate
    s = np.argsort(ston_found)
    # reverse s/n order

    ston_found = ston_found[s[::-1]]
    lamlines_found = lamlines_found[s[::-1]]
    index_of_strongest_line = 0
    lamline = lamlines_found[index_of_strongest_line]

    # MDR 2022/06/10 - Changed the first guess line to [O II] for MUSE.
    # KVN 2023/09/01 - Changed first guess back to Halpha for PASSAGE
    zguess = lamline / ha_6565_vac - 1
    # zguess = lamline / ((o2_3727_vac + o2_3730_vac) / 2.0) - 1.0
    # zguess = (lamline / o2_3730_vac) - 1.0
    # fwhm is defined for the red side, regardless of where line is
    fwhm_guess = 2.35 * a_image * config_pars["dispersion_red"]

    if stored_fits != False:
        first_stored_fit = stored_fits[0]
        users = [path.split("/")[-3].split("_")[-1] for path in stored_fits]
        fileObject = open(first_stored_fit, "r")
        alldata = pickle.load(fileObject)
        config_pars = alldata[10]
        fitresults_old = alldata[8]
        zguess = fitresults_old["redshift"]
        fwhm_guess = fitresults_old["fwhm_g141"]
        print("using stored fit from: " + users[0])
        print("available stored fits: ")
        print(users)
        ### also need to figure out what else to add?
        ### config pars for nodes can also be entered here.

    ### replace this with printouts from pickle files

    # print object info to screen
    if rejectPrevFit:
        print(" ")
        print_prompt("=" * 72)
        print_prompt(f"Par{par} Obj {obj:d}:")
        print_prompt("Initial redshift guess: z = %f" % (zguess))
        print_prompt(
            "\nWhat would you like to do with this object?\nSee the README for options, or type 'h' to print them all to the screen."
        )

    comment = ""

    # set this to loop over flags created at top of code.
    contamflags = {
        "la_1216": 0,
        "n5_1238": 0,
        "n5_1242": 0,
        "c4_1548": 0,
        "c4_1550": 0,
        "h2_1640": 0,
        "o3_1660": 0,
        "o3_1666": 0,
        "s3_1883": 0,
        "s3_1892": 0,
        "c3_1907": 0,
        "c3_1909": 0,
        "m2_2796": 0,
        "m2_2803": 0,
        "o2_3727": 0,
        "o2_3730": 0,
        "hg_4342": 0,
        "o3_4363": 0,
        "h2_4686": 0,
        "hb_4863": 0,
        "o3_4959": 0,
        "o3_5007": 0,
        "o1_6300": 0,
        "o1_6363": 0,
        "n2_6550": 0,
        "ha_6565": 0,
        "n2_6585": 0,
        "s2_6716": 0,
        "s2_6731": 0,
        "s3_9069": 0,
        "s3_9532": 0,
        "he10830": 0,
        "pg_10941": 0,
        "pb_12822": 0,
        "pa_18756": 0,
        "cont": 0,
    }  # MDR 2022/07/22

    # Skip if previous fit is to be accepted
    done = 0 if rejectPrevFit else 1
    fast_fit = False  # MDR 2022/06/30 - move to configuration file?
    orientation = None

    while done == 0:
        # KVN: creating spectra for each orientation ("T" for total/combined, "R" for row, "C" for column)
        # Note - do not move below outside the while loop. The masking (and probably other features) will not work properly
        spdata_T = trim_spec(tab_blue, tab_mid, tab_red, config_pars, mask_zeros=True, return_masks=True)

        spdata_R = trim_spec(tab_blue_R, tab_mid_R, tab_red_R, config_pars, mask_zeros=True, return_masks=True)

        spdata_C = trim_spec(tab_blue_C, tab_mid_C, tab_red_C, config_pars, mask_zeros=True, return_masks=True)

        spdata_T_contam = trim_spec(tab_blue_cont, tab_mid_cont, tab_red_cont, config_pars, mask_zeros=True, return_masks=True)

        spdata_R_contam = trim_spec(tab_blue_R_cont, tab_mid_R_cont, tab_red_R_cont, config_pars, mask_zeros=True, return_masks=True)

        spdata_C_contam = trim_spec(tab_blue_C_cont, tab_mid_C_cont, tab_red_C_cont, config_pars, mask_zeros=True, return_masks=True)

        # KVN: the default is to use the combined spectra as this will likely be the preferred option for most cases.
        # This will be updated when/if the user makes a different selection
        if orientation is None:
            spdata = copy(spdata_T)
        spec_lam = spdata[0]; spec_val = spdata[1]; spec_unc = spdata[2]; spec_con = spdata[3]; spec_zer = spdata[4]; mask_flg = spdata[5]
        
        # sticking with the while loop to determine whether user is finished with object
        # get spectrum for obj. do this every time because sometimes we
        # re-read with a mask or a different transition wavelength

        # KVN: trim_spec has been updated to take 3 grism filters
        # and the code now does this for the row & column (R & C)
        plot_chooseSpec(spdata_T, spdata_R, spdata_C, config_pars, plottitle, outdir)
        print_prompt("If you would like to change which spectrum is being fit, the options are: grismR, grismC, grismRcontam, grismCcontam, CombContam, or Comb to go back to the combined at any time. ")

        # Determine the largest extent of the object so broadening of the lines can be accounted for in the fitting. MDR 2022/06/30
        ab_image_max = np.max([objinfo["a_image"][0], objinfo["b_image"][0]])

        # apply the mask to the wavelength array
        masked_spec_lam = np.ma.masked_where(np.ma.getmask(spec_val), spec_lam)
        
        # compress the masked arrays for fitting
        
        print('Running the fit with the following settings: redshift = ',zguess,', fast_fit = ',fast_fit,', comp_fit = ',comp_fit,', polycont_fit = ',polycont_fit,', lincont_fit = ',lincont_fit)
        fit_inputs = [
            np.ma.compressed(masked_spec_lam),
            np.ma.compressed(spec_val),
            np.ma.compressed(spec_unc),
            config_pars,
            zguess,
            fwhm_guess,
            str(obj),
            ab_image_max,
            fast_fit,
            comp_fit,
            polycont_fit, 
            lincont_fit] 
        # parsing the input to facilitate parallel processing when fitting is done in batch mode.
        fitresults = fit_obj(fit_inputs)
        zfit = fitresults["redshift"]
        fitpars = fitresults["fit_parameters"]
        fitpars_nolines = cp.deepcopy(fitpars)
        fitpars_onlybroad = cp.deepcopy(fitpars)

        ############################################################################
        """
        The model parameters for the emission line amplitudes must be set to zero
        for the continuum fit, while the values for the line ratios must be non-zero
        to avoid division by zero. The 'get_fitpar_indices()' and 'get_ratio_indices()'
        functions below are defined in fitting.py to send these model parameter indices
        back to measure_z_interactive().
        """
        first_line_index, first_node_index = passage.get_fitpar_indices()
        ####### KVN -- line below doesn't work (idk why?), so hard-coding the broad line index. Will need to update. 
        # first_broad_line_index = passage.get_broad_indices()
        first_broad_line_index = 51
        
        fitpars_onlybroad[first_line_index:first_broad_line_index] = 0.0 
        fitpars_nolines[first_line_index:first_node_index] = 0.0

        for idx in get_ratio_indices():
            fitpars_nolines[idx] = 1.0
        for idx in get_ratio_indices():
            if idx < first_broad_line_index:
                fitpars_onlybroad[idx] = 1.0

        
        ############################################################################

        fitmodel = (
            emissionline_model(fitpars, np.ma.compressed(masked_spec_lam), comp_fit, polycont_fit, lincont_fit)
            * fitresults["scl_factor"])
        
        fitmodel_broad_gauss_fit = (
            emissionline_model(fitpars_onlybroad, np.ma.compressed(masked_spec_lam), comp_fit, polycont_fit, lincont_fit)
            * fitresults["scl_factor"])
        
        contmodel = (
            emissionline_model(fitpars_nolines, np.ma.compressed(masked_spec_lam), comp_fit, polycont_fit, lincont_fit)
            * fitresults["scl_factor"])

        # the fitting is done on compressed arrays, so we need to
        # create masked versions of the fit and continuum models
        full_fitmodel = np.zeros(spec_lam.shape, dtype=float)
        full_contmodel = np.zeros(spec_lam.shape, dtype=float)
        broad_fitmodel = np.zeros(spec_lam.shape, dtype=float)
        
        full_fitmodel[np.ma.nonzero(spec_val)] = fitmodel
        full_contmodel[np.ma.nonzero(spec_val)] = contmodel
        broad_fitmodel[np.ma.nonzero(spec_val)] = fitmodel_broad_gauss_fit
        
        full_fitmodel = np.ma.masked_where(np.ma.getmask(spec_val), full_fitmodel)
        full_contmodel = np.ma.masked_where(np.ma.getmask(spec_val), full_contmodel)
        broad_fitmodel = np.ma.masked_where(np.ma.getmask(spec_val), broad_fitmodel)

        # loop over the lines specified in 'flux_strings' and save the results to an array.
        # ---------
        # Added by KVN 13-Aug-2024 because double gaussian fit has extra lines. 
        if comp_fit == True: flux_strings = flux_strings_2gauss
        elif comp_fit == False: flux_strings = flux_strings_1gauss
        # ---------

        snr_meas_array = []
        for line in flux_strings:
            snr_meas_array.append(fitresults[line + "_flux"] / fitresults[line + "_error"])
        snr_meas_array = np.array(snr_meas_array)

        signal_lines = []
        for line in flux_strings:
            signal_lines.append(fitresults[line + "_flux"])
        signal_lines = np.array(signal_lines)

        err_lines = []
        for line in flux_strings:
            err_lines.append(fitresults[line + "_error"])
        err_lines = np.array(err_lines)

        # calculate a weighted s/n ratio for all lines and print on the plot.
        # signal_lines = np.array([fitresults['o2_total_flux'], fitresults['hg_4342_flux'], fitresults['hb_4863_flux'],fitresults['ha_total_flux'], fitresults['s2_total_flux']])
        # err_lines   = np.array([fitresults['o2_total_error'], fitresults['hg_4342_error'], fitresults['hb_4863_error'], fitresults['ha_total_error'], fitresults['s2_total_error']])

        w = np.where(signal_lines > 0)

        # MDR 2022/06/10 - updated the definition of weight SNR to be sum of SNRs that are  > 3 weighted by their contribution to the total flux.
        total_flux = 0
        for line in flux_strings:
            snr_line = fitresults[line + "_flux"] / fitresults[line + "_error"]
            if snr_line >= 3.0:
                total_flux = total_flux + fitresults[line + "_flux"]
        # MDR 2022/06/10
        snr_tot_others = []
        for line in flux_strings:
            snr_line = fitresults[line + "_flux"] / fitresults[line + "_error"]
            if snr_line >= 3.0:
                snr_tot_weight = fitresults[line + "_flux"] / total_flux
                if np.isfinite(snr_line * snr_tot_weight):
                    snr_tot_others.append(snr_line * snr_tot_weight)
                else:
                    snr_tot_others.append(0.0)
        # MDR 2022/06/10
        snr_tot_others = np.sum(snr_tot_others)

        # plot the whole darn thing
        plot_object(
            zguess,
            fitresults["redshift"],
            spdata,
            config_pars,
            snr_meas_array,
            snr_tot_others,
            full_fitmodel,
            full_contmodel,
            broad_fitmodel,
            lamline,
            lamlines_found,
            index_of_strongest_line,
            contmodel,
            plottitle,
            outdir, 
            orientation='Combined')
        
        #        print "    Guess Redshift: z = %f" % (zguess)
        print_prompt("    Fit Redshift:   z = %f\n" % (zfit))
        # input
        option = input("> ")

        # checking user's input. keeping this format the same as before
        # any time done is set to 1, the object is considered fit

        # reject object
        if option.strip().lower() == "r":
            zset = 0
            if len(comment) > 0:
                comment = "rejected" + ", " + comment
            else:
                comment = "rejected"
            done = 1

        # # accept object
        # elif option.strip().lower() == 'a':
        #     done = 1
        #     zset = 1
        #     flagcont = 1
        # accept object MDR 2022/06/30
        elif option.strip().lower() == "a":
            if comp_fit == True:
                if fast_fit == False:
                    done = 1
                    zset = 1
                    flagcont = 1
                elif fast_fit == True:
                    print("\nWARNING: Still using fast fit mode, type full for refined fit.")

            elif comp_fit == False:
                print("\n--- You are accepting the single Gaussian fit ---")
                if fast_fit == False:
                    done = 1
                    zset = 1
                    flagcont = 1
                elif fast_fit == True:
                        print("\nWARNING: Still using fast fit mode, type full for refined fit.")


        ### KVN 05-Aug-2024
        ### Adding option to fit double gaussian to emission lines:
        elif option.strip().lower() == "2gauss":
            print('Fitting emission lines as double gaussians.')
            comp_fit = True

        ### KVN 12-Aug-2024
        ### Adding option to go back to 1 gaussian fit after selecting 2 gaussian fit
        elif option.strip().lower() == "1gauss":
            comp_fit = False


        ### KVN 06-Aug-2024
        ### Adding option to fit continuum as a polynomial
        elif option.strip().lower() == "polycont":
            polycont_fit = True; lincont_fit = False

        ### KVN 06-Aug-2024
        ### Adding option to fit continuum as a line
        elif option.strip().lower() == "lincont":
            lincont_fit = True; polycont_fit = False

        ### KVN 06-Aug-2024
        ### Adding option to fit continuum as a line
        elif option.strip().lower() == "splinecont":
            lincont_fit = False; polycont_fit = False

        ### MDR 2022/06/30
        elif option.strip().lower() == "full":
            fast_fit = False
        elif option.strip().lower() == "fast":
            fast_fit = True
        ### MDR 2022/06/30

        # accept object and note contamination
        elif option.strip().lower() == "ac":
            done = 1
            zset = 1
            flagcont = 2
            # add to contamination flags
            for k, v in contamflags.items():
                contamflags[k] = contamflags[k] | 1

        # change redshift guess
        elif option.strip().lower() == "z":
            print_prompt("The current redshift guess is: %f\nEnter Redshift Guess:" % zguess)
            try: zguess = float(input("> "))
            except ValueError: print_prompt("Invalid Entry.")

        # change wavelength guess
        elif option.strip().lower() == "w":
            print_prompt("The current emission line wavelength is: %f\nEnter Wavelength Guess in Angstroms:" % lamline)
            # save previous line guess (if not Ha)
            old_rest_wave = lamline / (1.0 + zguess)
            try:
                newwave = float(input("> "))
            except ValueError:
                print_prompt("Invalid Entry.")
            else:
                zguess = newwave / old_rest_wave - 1.0
                lamline = newwave

        # change the fwhm guess
        elif option.strip().lower() == "fw":
            print_prompt("Enter a Guess for FWHM in pixels")
            print_prompt(
                "The current fwhm_fit is:  "
                + str(fitresults["fwhm_g141"] / config_pars["dispersion_red"])
                + " and 2*A_image is: "
                + str(2 * a_image[0])
            )
            try:
                fwhm_guess = config_pars["dispersion_red"] * float(input("> "))
                if fwhm_guess <= 0:
                    print_prompt("FWHM (in pixels) must be larger than zero. Keeping previous FWHM.")
                    fwhm_guess = fitresults["fwhm_g141"] / config_pars["dispersion_red"]
            except ValueError:
                print_prompt("Invalid Entry.")

        # change delta_z range
        elif option.strip().lower() == "dz":
            print_prompt("Enter a new delta_z limit")
            print_prompt("The current delta_z is:  " + str(config_pars["delta_z"]))
            try:
                config_pars["delta_z"] = float(input("> "))
            except ValueError:
                print_prompt("Invalid Entry.")

        # mask out up to 8 regions of the spectrum. m1 & m2 assigned to gaps between grism filters
        # but can be altered here
        elif option.strip().lower() == "m1":
            print_prompt('Masked Region 1 currently masks: '+str(config_pars["mask_region1"]))
            print_prompt("Enter wavelength window to mask out:  blue, red:")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region1"] = maskwave

        elif option.strip().lower() == "m2":
            print_prompt('Masked Region 2 currently masks: '+str(config_pars["mask_region2"]))
            print_prompt("Enter wavelength window to mask out:  blue, red:")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region2"] = maskwave

        elif option.strip().lower() == "m3":
            print_prompt('Masked Region 3 currently masks: '+str(config_pars["mask_region3"]))
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region3"] = maskwave

        elif option.strip().lower() == "m4":
            print_prompt('Masked Region 4 currently masks: '+str(config_pars["mask_region4"]))
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region4"] = maskwave

        elif option.strip().lower() == "m5":
            print_prompt('Masked Region 5 currently masks: '+str(config_pars["mask_region5"]))
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region5"] = maskwave

        elif option.strip().lower() == "m6":
            print_prompt('Masked Region 6 currently masks: '+str(config_pars["mask_region6"]))
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region6"] = maskwave

        elif option.strip().lower() == "m7":
            print_prompt('Masked Region 7 currently masks: '+str(config_pars["mask_region7"]))
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region7"] = maskwave

        elif option.strip().lower() == "m8":
            print_prompt('Masked Region 8 currently masks: '+str(config_pars["mask_region8"]))
            print_prompt("Enter wavelength window to mask out:  blue, red (Angstroms):")
            maskstr = input("> ")
            try:
                maskwave = [float(maskstr.split(",")[0]), float(maskstr.split(",")[1])]
            except (IndexError, ValueError):
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                config_pars["mask_region8"] = maskwave

        # change the transition wavelength between the grisms
        elif option.strip().lower() == "t1":
            print_prompt(
                "The current transition wavelength is: "
                + str(config_pars["transition_wave1"])
                + "\nEnter the wavelength for the G115 to G150 transition:"
            )
            try:
                config_pars["transition_wave1"] = float(input("> "))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelength of grism transition.")

        # change the transition wavelength between the grisms
        elif option.strip().lower() == "t2":
            print_prompt(
                "The current transition wavelength is: "
                + str(config_pars["transition_wave2"])
                + "\nEnter the wavelength for the G150 to G200 transition:"
            )
            try:
                config_pars["transition_wave2"] = float(input("> "))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelength of grism transition.")

        # change the nodes used for the continuum spline
        elif option.strip().lower() == "nodes":
            strnw = ",".join(str(nw) for nw in config_pars["node_wave"])
            print_prompt("Enter Wavelengths for Continuum Spline: w1, w2, w3, w4, ....")
            print_prompt("current node wavelengths are: %s" % strnw)
            nodestr = input("> ")
            nodesplit = nodestr.split(",")
            node_arr = []
            try:
                for nodelam in nodesplit:
                    node_arr.append(float(nodelam))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            else:
                node_arr = np.array(node_arr)
                # sort by wavelength
                node_arr = np.sort(node_arr)
                config_pars["node_wave"] = node_arr

        # added by KVN June 2024
        # REMOVE any of the nodes used for the continuum spline by wavelength
        elif option.strip().lower() == "rmnodes":
            strnw = ",".join(str(nw) for nw in config_pars["node_wave"])
            print_prompt(
                "Enter wavelength of the nodes you would like to REMOVE: w1, w2, w3, ..."
            )
            print_prompt("current node wavelengths are: %s" % strnw)
            nodestr = input("> ")
            nodesplit = nodestr.split(",")
            node_arr_rem = []
            try:
                for nodelam in nodesplit:
                    node_arr_rem.append(float(nodelam))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            node_arr_rem = np.array(node_arr_rem)
            # check if nodes to be removed are actually there to begin with
            unmatched_nodes = [
                k for k in node_arr_rem if k not in config_pars["node_wave"]
            ]
            if len(unmatched_nodes) != 0:
                print_prompt("At least one of the nodes at the wavelength entered does not exist. Returning same fit.")
            else:
                # remove the nodes at entered wavelenths
                node_arr = [k for k in config_pars["node_wave"] if k not in node_arr_rem]
                # sort by wavelength
                node_arr = np.sort(node_arr)
                config_pars["node_wave"] = node_arr

        # added by KVN June 2024
        # ADD nodes used for the continuum spline
        elif option.strip().lower() == "addnodes":
            strnw = ",".join(str(nw) for nw in config_pars["node_wave"])
            upper_wav_bound = config_pars["lambda_max"]
            lower_wav_bound = config_pars["lambda_min"]
            print_prompt(
                "Enter wavelength of the nodes you would like to ADD: w1, w2, w3, ..."
            )
            print_prompt("current node wavelengths are: %s" % strnw)
            nodestr = input("> ")
            nodesplit = nodestr.split(",")
            node_arr_add = []
            print
            try:
                for nodelam in nodesplit:
                    node_arr_add.append(float(nodelam))
            except ValueError:
                print_prompt(
                    "Invalid entry. Enter wavelengths separated by commas"
                )
            node_arr_add = np.array(node_arr_add)
            # remove any added nodes that are outside the blue/red cutoff
            node_goodadd = [
                k for k in node_arr_add if k <= upper_wav_bound and k >= lower_wav_bound
            ]
            if len(node_goodadd) != len(node_arr_add):
                print_prompt(
                    "At least one of the nodes is outside the wavelength range and has not been added."
                )
                node_arr_add = np.copy(node_goodadd)
            # add the new nodes simply by turning arrays into lists
            node_arr = np.array(list(config_pars["node_wave"]) + list(node_arr_add))
            # sort by wavelength
            node_arr = np.sort(node_arr)
            config_pars["node_wave"] = node_arr

        # added by KVN June 2024
        # SHIFT ALL nodes used for the continuum spline by some wavelength
        elif option.strip().lower() == "shiftallnodes":
            strnw = ",".join(str(nw) for nw in config_pars["node_wave"])
            upper_wav_bound = config_pars["lambda_max"]
            lower_wav_bound = config_pars["lambda_min"]
            print_prompt("Enter wavelength by which you would like to SHIFT ALL nodes (in Angstrom(?)): w1")
            print_prompt("current node wavelengths are: %s" % strnw)
            try:
                shift_wav = float(input("> "))
            except ValueError:
                print_prompt("Invalid entry. Enter single wavelength")
            node_arr = [k + shift_wav for k in config_pars["node_wave"]]
            # sort by wavelength
            node_arr = np.sort(node_arr)
            ### note if your shifted nodes are outside the blue/red cutoff,
            ### then the nodes gets placed at the cutoff,
            ### which is usually a bad place to measure the continuum.
            # Remove any nodes that fall outside the blue/red cutoff:
            node_good = [
                k for k in node_arr if k <= upper_wav_bound and k >= lower_wav_bound
            ]
            if len(node_good) != len(node_arr):
                print_prompt(
                    "At least one of the nodes is outside the wavelength range and is removed."
                )
                config_pars["node_wave"] = node_good
            else:
                config_pars["node_wave"] = node_arr

        # added by KVN June 2024
        # SHIFT INDIVIDUAL nodes used for the continuum spline by some wavelength
        elif option.strip().lower() == "shiftnodes":
            strnw = ",".join(str(nw) for nw in config_pars["node_wave"])
            upper_wav_bound = config_pars["lambda_max"]
            lower_wav_bound = config_pars["lambda_min"]
            print_prompt("Enter wavelength(s) of node(s) you would like to shift: w1, w2, w3...")
            print_prompt("current node wavelengths are: %s" % strnw)
            nodestr = input("> ")
            nodesplit = nodestr.split(",")
            node_arr_shift = []
            try:
                for nodelam in nodesplit:
                    node_arr_shift.append(float(nodelam))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelengths separated by commas")
            print_prompt("Enter wavelength by which you would like to SHIFT this node: w1")
            try:
                shift_wav = float(input("> "))
            except ValueError:
                print_prompt("Invalid entry. Enter single wavelength")
            # find the indices of the nodes to be shifted & shift nodes at those indices by input ('shift_wav')
            nodes_index = [k for k in range(len(config_pars["node_wave"])) if config_pars["node_wave"][k] in node_arr_shift]
            node_arr = np.copy(config_pars["node_wave"])
            node_arr[nodes_index] += shift_wav
            # sort by wavelength -- this is important as shift may cause order to change
            node_arr = np.sort(node_arr)
            ### note if your shifted nodes are outside the blue/red cutoff, then the nodes gets placed at the cutoff,
            ### which is usually a bad place to measure the continuum.
            # Remove any nodes that fall outside the blue/red cutoff:
            node_good = [
                k for k in node_arr if k <= upper_wav_bound and k >= lower_wav_bound
            ]
            if len(node_good) != len(node_arr):
                print_prompt(
                    "At least one of the nodes is outside the wavelength range and is removed."
                )
                config_pars["node_wave"] = node_good
            else:
                config_pars["node_wave"] = node_arr


        # added by KVN July 2024
        # Choose which spectrum should be used for the line fitting.
        elif option.strip().lower() == "grismr":
            spdata_lam_org = spdata[0]
            if len(spdata_lam_org) == len(spdata_R[0]):
                spdata = spdata_R
                orientation = "R"
            else: print("The selected R orientation appears to have missing data, keeping the original selection.")
                
        elif option.strip().lower() == "grismrcontam":
            spdata_lam_org = spdata[0]
            if len(spdata_lam_org) == len(spdata_R_contam[0]):
                spdata = spdata_R_contam
                orientation = "RContam"
            else: print("The selected R(with contamination) orientation appears to have missing data, keeping the original selection.")
                
        elif option.strip().lower() == "grismc":
            spdata_lam_org = spdata[0]
            if len(spdata_lam_org) == len(spdata_C[0]):
                spdata = spdata_C
                orientation = "C"
            else: print("The selected C orientation appears to have missing data, keeping the original selection.")
                
        elif option.strip().lower() == "grismccontam":
            spdata_lam_org = spdata[0]
            if len(spdata_lam_org) == len(spdata_C_contam[0]):
                spdata = spdata_C_contam
                orientation = "CContam"
            else: print("The selected C(with contamination) orientation appears to have missing data, keeping the original selection.")

        elif option.strip().lower() == "comb":
            spdata = spdata_T
            orientation = "Combined"

        elif option.strip().lower() == "combcontam":
            spdata = spdata_T_contam
            orientation = "CombinedContam"
        
        elif option.strip().lower() == "user":

            if stored_fits != False:
                print_prompt("Enter name of user for toggling between stored fits")
                user_input = input("> ")
                try:
                    w = users.index(user_input)
                except ValueError:
                    print_prompt("Invalid entry. Enter a valid user name.")

                different_stored_fit = stored_fits[w]
                fileObject = open(different_stored_fit, "r")
                alldata = pickle.load(fileObject)
                config_pars = alldata[10]
                fitresults_old = alldata[8]
                zguess = fitresults_old["redshift"]
                fwhm_guess = fitresults_old["fwhm_g141"]
            else:
                print("there are no stored fits!")

        # reset all options
        elif option.strip().lower() == "reset":
            print_prompt(
                "Reset configuration parameters, fwhm guess, and zguess to default values"
            )
            print('Resetting parameters to those specified in: ', path_to_code+"/default.config")
            config_pars = read_config(path_to_code+"/default.config", availgrism=availgrism)
            fwhm_guess = 2.35 * a_image * config_pars["dispersion_red"]
            # reset strongest line, too
            index_of_strongest_line = 0
            lamline = lamlines_found[index_of_strongest_line]
            # zguess = lamline / ha_6565_vac - 1
            zguess = lamline / ((o2_3727_vac + o2_3730_vac) / 2.0) - 1.0
            # zguess = (lamline / o2_3730_vac) - 1.0
            # reset contamflags
            for k, v in contamflags.items():
                contamflags[k] = contamflags[k] & 0
            ### if use stored = true, this should set us back to using the pickle file values

        # change the blue cutoff of G102 (or whichever grism is present?)
        elif option.strip().lower() == "bluecut":
            print_prompt(
                "The current blue cutoff is: "
                + str(config_pars["lambda_min"])
                + "\nChange the blue cutoff of the F115W grism filter:"
            )
            try:
                config_pars["lambda_min"] = float(input("> "))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelength of blue cutoff.")

        # change the red cutoff of G141
        elif option.strip().lower() == "redcut":
            print_prompt(
                "The current red cutoff is: "
                + str(config_pars["lambda_max"])
                + "\nChange the red cutoff of the F200W grism filter:"
            )
            try:
                config_pars["lambda_max"] = float(input("> "))
            except ValueError:
                print_prompt("Invalid entry. Enter wavelength of red cutoff.")

        # change to next brightest line
        elif option.strip().lower() == "n":
            nlines_found_cwt = np.size(lamlines_found)
            index_of_strongest_line = index_of_strongest_line + 1
            if index_of_strongest_line < (nlines_found_cwt):
                lamline = lamlines_found[index_of_strongest_line]
                # zguess = lamline / 6564.610 - 1.0
                zguess = lamline / ((o3_5007_vac)) - 1.0
                # zguess = (lamline / o2_3730_vac) - 1.0
            else:
                print_prompt(
                    "There are no other automatically identified peaks. Select another option."
                )
                # stay at current line
                index_of_strongest_line -= 1

        # change to previous brightest line
        # added KVN 30-Jan-2025
        elif option.strip().lower() == "p":
            nlines_found_cwt = np.size(lamlines_found)
            if index_of_strongest_line >= 1:
                index_of_strongest_line = index_of_strongest_line - 1
            else: print_prompt("Already at the strongest peak.")
            if nlines_found_cwt > 1 and index_of_strongest_line >= 0:
                lamline = lamlines_found[index_of_strongest_line]
                # zguess = lamline / 6564.610 - 1.0
                zguess = lamline / ((o3_5007_vac)) - 1.0
                # zguess = (lamline / o2_3730_vac) - 1.0
            else:
                print_prompt(
                    "There are no other automatically identified peaks. Select another option."
                )
                # stay at current line
                index_of_strongest_line -= 1

        # change to another line
        elif option.strip().lower() == "lya":
            zguess = (lamline / la_1216_vac) - 1.0
        elif option.strip().lower() == "c4":
            zguess = (lamline / c4_1548_vac) - 1.0
        elif option.strip().lower() == "o2":
            zguess = (lamline / o2_3730_vac) - 1.0
        elif option.strip().lower() == "hb":
            zguess = (lamline / hb_4863_vac) - 1.0
        elif option.strip().lower() == "hg":
            zguess = (lamline / hg_4342_vac) - 1.0
        elif option.strip().lower() == "o31":
            zguess = (lamline / o3_4959_vac) - 1.0
        elif option.strip().lower() == "o32":
            zguess = (lamline / o3_5007_vac) - 1.0
        elif option.strip().lower() == "ha":
            zguess = (lamline / ha_6565_vac) - 1.0
        elif option.strip().lower() == "s2":
            zguess = (lamline / s2_6716_vac) - 1.0
        elif option.strip().lower() == "s31":
            zguess = (lamline / s3_9069_vac) - 1.0
        elif option.strip().lower() == "s32":
            zguess = (lamline / s3_9532_vac) - 1.0
        elif option.strip().lower() == "pb":
            zguess = (lamline / pb_12822_vac) - 1.0
        elif option.strip().lower() == "pa":
            zguess = (lamline / pa_18756_vac) - 1.0
        elif option.strip().lower() == "mg":
            zguess = (lamline / m2_2796_vac) - 1.0


        # note contamination
        elif option.strip().lower() == "contam":
            print_prompt(
                "Specify contamination.\nEnter a comma-separated list of identifiers choosing from:\n"
                + contam_flags_string
            )
            cf = input("> ")
            cflags = [thing.strip().lower() for thing in cf.split(",")]
            # continuum contamination sets bit 1 for all lines and the continuum itself
            if "c" in cflags:
                for k, v in contamflags.items():
                    contamflags[k] = contamflags[k] | 2
            cflaglines = [thing for thing in cflags if thing != "c"]
            # specific line contamination sets bit 2 for all lines
            for contamflag in cflaglines:
                try:
                    contamflags[contamflag] = contamflags[contamflag] | 4
                except KeyError:
                    print_prompt("{} not known. Skipping".format(contamflag))
            # sqlite3 database support - automatically creates and initializes DB if required
            # databaseManager.setFlags(par, obj, [(flagName, flagValue) for flagName, flagValue in contamflags.iteritems()])

        # add a comment
        elif option.strip().lower() == "c":
            print_prompt("Enter your comment here:")
            if len(comment) > 0:
                comment = input("> ") + ", " + comment
            else:
                comment = input("> ")

            # sqlite3 database support - automatically creates and initializes DB if required
        # databaseManager.saveAnnotation((par, obj, comment.decode('utf-8')))

        # set or unset one or more flags
        elif option.strip().lower() == "flag":
            print_prompt(
                "Enter a comma-separated list of flag, value pairs e.g. CONTAM, 1, CONTIN, 2:"
            )
            print_prompt(
                "Valid flags are {}".format(
                    WISPLFDatabaseManager.WISPLFDatabaseManager.validMutableFlags
                )
            )
            flagList = input("> ")
            # sqlite3 database support - automatically creates and initializes DB if required
            # databaseManager.setFlagsFromString(par, obj, flagList.decode('utf-8'))

        # write object summary
        elif option.strip().lower() == "s":
            write_object_summary(
                par,
                obj,
                fitresults,
                snr_meas_array,
                contamflags,
                comp_fit,
                summary_type="working",
            )

        # print help message
        elif option.strip().lower() == "h":
            print_help_message()

        ### image/display options ###
        # change 2d stamp scaling to linear
        elif option.strip().lower() == "lin":
            if g102zeros is not None:
                show2dNEW(
                    "F115W",
                    par,
                    obj,
                    g102zeros,
                    user,
                    "linear",
                    path_to_data=path_to_data,
                )
            if g141zeros is not None:
                show2dNEW(
                    "F150W",
                    par,
                    obj,
                    g141zeros,
                    user,
                    "linear",
                    path_to_data=path_to_data,
                )

        # change 2d stamp scaling to log
        elif option.strip().lower() == "log":
            if g102zeros is not None:
                show2dNEW(
                    "F115W",
                    par,
                    obj,
                    g102zeros,
                    user,
                    "log",
                    path_to_data=path_to_data,
                )
            if g141zeros is not None:
                show2dNEW(
                    "F150W",
                    par,
                    obj,
                    g141zeros,
                    user,
                    "log",
                    path_to_data=path_to_data,
                )

        # change g102 2d stamp scaling to zscale
        elif option.strip().lower() == "zs102":
            print_prompt("Enter comma-separated range for G102 zscale: z1,z2")
            zscale = input("> ")
            zs = zscale.split(",")
            try:
                z1 = float(zs[0])
                z2 = float(zs[1])
            except ValueError:
                print_prompt("Invalid entry.")
            else:
                if g102zeros is not None:
                    show2dNEW(
                        "F115W",
                        par,
                        obj,
                        g102zeros,
                        user,
                        "linear",
                        zran1=z1,
                        zran2=z2,
                        path_to_data=path_to_data,
                    )

        # change g141 2d stamp scaling to zscale
        elif option.strip().lower() == "zs141":
            print_prompt("Enter comma-separated range for G141 zscale: z1,z2")
            zscale = input("> ")
            zs = zscale.split(",")
            try:
                z1 = float(zs[0])
                z2 = float(zs[1])
            except ValueError:
                print_prompt("Invalid entry.")
            else:
                if g141zeros is not None:
                    show2dNEW(
                        "F150W",
                        par,
                        obj,
                        g141zeros,
                        user,
                        "linear",
                        zran1=z1,
                        zran2=z2,
                        path_to_data=path_to_data,
                    )

        # recenter full images
        # Updated for PASSAGE KVN 19-Nov-2024
        elif option.strip().lower() == "dc":
            panDirect_PASSAGE(ra, dec)
            if show_dispersed:  # MB
                panDispersed_PASSAGE(obj, parno=par, path_to_data=path_to_data)

        # reload full iamges
        # Updated for PASSAGE KVN 19-Nov-2024
        elif option.strip().lower() == "reload":
            showDirect_PASSAGE(par, path_to_data=path_to_data)
            showSpec2D_PASSAGE(par, obj, path_to_data=path_to_data)

        # reload direct image region files
        elif option.strip().lower() == "dr":
            reloadReg()

        # new options dealing with iterating objects
        # can't actually go back or choose another object now,
        # but allow for them now just in case
        elif option.strip().lower() == "b":
            print_prompt("Please either reject or accept this object first.")
        elif "obj" in option:
            print_prompt("Please either reject or accept this object first.")
        # print remaining objects that have not yet been inspected
        elif option.strip().lower() == "left":
            print_prompt("Remaining objects:")
            print(remaining)
        # print all objects in line list
        elif option.strip().lower() == "list":
            print_prompt("All objects:")
            print(allobjects)

        # quit this object
        elif option.strip().lower() == "q":
            print_prompt("Quitting Obj %i. Nothing saved to file" % (obj))
            print_prompt("-" * 72)
            return 0

        # catch-all for everything else
        else:
            print_prompt("Invalid entry.  Try again.")

    # only re-save data if the previous fit was discarded
    if rejectPrevFit:
        # plot the whole darn thing
        plot_object(
            zguess,
            fitresults["redshift"],
            spdata,
            config_pars,
            snr_meas_array,
            snr_tot_others,
            full_fitmodel,
            full_contmodel,
            broad_fitmodel,
            lamline,
            lamlines_found,
            index_of_strongest_line,
            contmodel,
            plottitle,
            outdir,
            zset=zset,
        )

        # write to file if object was accepted
        if zset == 1:
            if np.any(spdata[1].mask):
                fitresults, snr_meas_array = check_masked_lines(
                    fitresults, snr_meas_array, spdata, flux_strings)

            # write object summary
            write_object_summary(par, obj, fitresults, snr_meas_array, contamflags, comp_fit)

            # sqlite3 database support - automatically creates and initializes DB if required
            # databaseManager.saveCatalogueEntry(databaseManager.layoutCatalogueData(par, obj, ra[0], dec[0], a_image[0],
            #                                                                       b_image[0], jmag[0], hmag[0], fitresults, flagcont))
            if comp_fit == False:
                writeToCatalog(
                    linelistoutfile,
                    par,
                    obj,
                    ra,
                    dec,
                    a_image,
                    b_image,
                    jmag,
                    hmag,
                    snr_tot_others,
                    fitresults,
                    contamflags,
                    comp_fit
                )
        
                writeFitdata(
                    fitdatafilename,
                    spec_lam,
                    spec_val,
                    spec_unc,
                    spec_con,
                    spec_zer,
                    full_fitmodel,
                    full_contmodel,
                    mask_flg,
                )
        
                fitspec_pickle = open(fitdatafilename + ".pickle", "wb")
                output_meta_data = [
                    par,
                    obj,
                    ra,
                    dec,
                    a_image,
                    b_image,
                    jmag,
                    hmag,
                    fitresults,
                    flagcont,
                    config_pars,
                ]
                pickle.dump(output_meta_data, fitspec_pickle)
                fitspec_pickle.close()
            elif comp_fit == True:
                writeToCatalog2gauss(
                    linelistoutfile, par, obj,ra, dec, a_image, b_image,
                    jmag, hmag, snr_tot_others, fitresults, contamflags, comp_fit)
    
                writeFitdata(
                    fitdatafilename,
                    spec_lam,
                    spec_val,
                    spec_unc,
                    spec_con,
                    spec_zer,
                    full_fitmodel,
                    full_contmodel,
                    mask_flg)
    
                fitspec_pickle = open(fitdatafilename + ".pickle", "wb")
                output_meta_data = [
                    par,
                    obj,
                    ra,
                    dec,
                    a_image,
                    b_image,
                    jmag,
                    hmag,
                    fitresults,
                    flagcont,
                    config_pars,
                ]
                pickle.dump(output_meta_data, fitspec_pickle)
                fitspec_pickle.close()
        
        # else:
        #     # done == 1, but zset == 0 => rejected
        #     databaseManager.saveCatalogueEntry(databaseManager.layoutCatalogueData(par, obj, ra[0], dec[0], a_image[0],
        #                                                                           b_image[0], jmag[0], hmag[0],
        #                                                                           None,
        #                                                                           None))
        #     databaseManager.setFlags(par, obj, [('REJECT', 1)])
        #     databaseManager.saveAnnotation((par, obj, 'REJECTED'))

        # write comments to file
        # if we go back to the previous objects, duplicate comments will still be
        # written
        writeComments(commentsfile, par, obj, comment)

        # write object to done file, incase process gets interrupted
        if not os.path.exists(outdir + "/done_%s" % user):
            f = open(outdir + "/done_%s" % user, "w")
        else:
            f = open(outdir + "/done_%s" % user, "a")
        f.write("%i\n" % obj)
        f.close()


def check_input_objid(objlist, objid, nextup):
    if verbose == True:
        print("\nRunning check_input_objid...\n")  # MDR 2022/05/17

    fulllist = ", ".join(["%i" % o for o in objlist])
    if objid not in objlist:
        print_prompt("Obj %i is not in the line list" % (objid), prompt_type="interim")
        print_prompt("Full line list: \n%s" % (fulllist), prompt_type="interim")
    while objid not in objlist:
        if nextup == 0:
            o = input("Please try again, or hit enter or 'q' to quit. > ")
        else:
            o = input(
                "Please try again, or hit enter to continue with Obj %s: > " % nextup
            )
        if o.strip().lower() == "q":
            return False
        elif isFloat(o.strip()):
            objid = int(o)
        elif "obj" in o:
            objid = int(re.search("\d+", o).group())
        elif o.strip() == "":
            return False
    return objid


def measure_z_interactive(
    linelistfile=" ",
    path_to_data=" ",
    path_to_code=" ",
    show_dispersed=True,
    path_to_stored_fits=" ",
    print_colors=True,
    parno=0):

    if verbose == True:
        print("\nRunning measure_z_interactive...\n")  # MDR 2022/05/17
    # turn off color printing to terminal if required
    if print_colors is False:
        global setcolors
        for k, v in setcolors.items():
            setcolors[k] = "\033[0m"

    if path_to_data == " ":
        ### running from the Spectra directory
        path_to_data = "../../"

    # if path_to_stored_fits == ' ':
    #    use_stored_fits  = False
    # elif os.path.exists(path_to_stored_fits) :
    #    use_stored_fits = True
    #    print 'looking for stored fit data'
    # else:
    #    use_stored_fits = False
    #    print 'not using stored fit data'

    #### STEP 0:   set ds9 window to tile mode ################################
    ###########################################################################
    # not the best way to do this, but matching the method in guis.py
    # cmd = "xpaset -p ds9 tile"
    # os.system(cmd)
    # if show_dispersed:
    #     cmd = "xpaset -p ds9 tile grid layout 3 2"
    # else:
    #     cmd = "xpaset -p ds9 tile grid layout 2 2"
    # os.system(cmd)

    #### STEP 1:   get linelist ###############################################
    ###########################################################################
    if linelistfile == " ":
        print('path_to_data= ', path_to_data, os.getcwd() )
        files = glob(path_to_data+"/linelist/Par"+str(parno)+"lines.dat")
        if len(files) == 0:
            print_prompt("No line list file found", prompt_type="interim")
            return 0
        else:
            linelistfile = files[0]
    if not os.path.exists(linelistfile):
        print_prompt(
            "Invalid path to line list file: %s" % (linelistfile), prompt_type="interim"
        )
        return 0
    else:
        print_prompt("Found line list file: %s" % (linelistfile), prompt_type="interim")

    #### STEP 1b:   read the list of candidate lines  ####################
    ###########################################################################

    llin = asciitable.read(
        linelistfile, names=["parnos", "grism", "objid", "wavelen", "npix", "ston"], converters={"parnos": str},
    )

    parnos = llin["parnos"]
    grism = llin["grism"]
    objid = llin["objid"]
    wavelen = llin["wavelen"]
    npix = llin["npix"]
    ston = llin["ston"]
    objid_unique = np.unique(objid)
    #par = parnos[0]  # MDR 2022/05/17
    par = parno # KVN 2024/07/31

    #### STEP 2:  set user name and output directory #########################
    ###########################################################################
    if verbose == True:
        print("")
        print("All outputs will be stored in...\n")  # MDR 2022/05/17
        print(os.getcwd())
        print("")

    tmp = glob(path_to_data + "Par" + str(par) + "/Spectra/*.dat")  # MDR 2022/05/17 and updated KVN 2024/07/31
    print_prompt(
        "You are about to inspect emission lines identified in parallel field {}".format(parno),
        prompt_type="interim",
    )
    print_prompt("Please enter your name or desired username", prompt_type="interim")
    while True:
        user = input("> ")
        if len(user) == 0:
            print_prompt("Username required", prompt_type="interim")
            continue
        else:
            break
    user = user.strip().lower()
    # create output directory
    outdir = "Par%s_output_%s" % (parno, user)
    if not os.path.isdir(outdir):
        os.mkdir(outdir)

    # if (verbose == True):
    #     print('tmp =', tmp) # MDR 2022/05/17
    #     print('user =', user) # MDR 2022/05/17
    #     print('outdir =', outdir) # MDR 2022/05/17

    #### STEP 3: define filenames and check for partially complete work #####
    #########################################################################
    if verbose == True:
        print("\nCreating figs and fitdata directories...\n")  # MDR 2022/05/17
    if not os.path.exists(os.path.join(outdir, "figs")):
        os.makedirs(os.path.join(outdir, "figs"))
    if not os.path.exists(os.path.join(outdir, "fitdata")):
        os.makedirs(os.path.join(outdir, "fitdata"))

    parts = os.path.splitext(os.path.basename(linelistfile))
    linelistoutfile = os.path.join(outdir, "%s_catalog_%s.dat" % (parts[0], user))
    commentsfile = os.path.join(outdir, "%s_comments_%s.dat" % (parts[0], user))
    # the file that will be used to determine which objects are "done"
    donefile = outdir + "/done_%s" % user

    # if (verbose == True):
    #     print('parts =', parts) # MDR 2022/05/17
    #     print('linelistoutfile =', linelistoutfile) # MDR 2022/05/17
    #     print('commentsfile =', commentsfile) # MDR 2022/05/17
    #     print('donefile =', donefile) # MDR 2022/05/17

    if os.path.isfile(linelistoutfile):
        print_prompt(
            "\nOutput file: \n  %s \nalready exists\n" % linelistoutfile,
            prompt_type="interim",
        )
        ask = input("Append? [Y/n] ")
        if ask.lower() == "n":
            os.unlink(linelistoutfile)
            os.unlink(commentsfile)
            # starting over, no objects have been done
            os.unlink(donefile)
            objid_done = np.array([])

            # sqlite3 database support - automatically creates and initializes DB if required
            # If the field has been previously examined, but those results are to be discarded,
            # then reset the database tables.
            # All Par numbers in the putative line list file should be the same, so the zeroth
            # element corresponds to the current field ID.
            # databaseManager = WDBM(dbFileNamePrefix=os.path.join(outdir,'Par{}'.format(parnos[0])))
            # databaseManager.resetDatabaseTables()
        else:
            # an object may be written to the comment file before it has
            # actually been inspected, so use donefile for a list
            # of the "done" objects
            objid_done = np.atleast_1d(np.genfromtxt(donefile, dtype=int))
    else:
        if os.path.exists(donefile):
            os.unlink(donefile)
        objid_done = np.array([], dtype=int)

    #### STEP 4: create trace.reg files ############################
    #########################################################################
    if verbose == True:
        print("Creating trace.reg files...\n")  # MDR 2022/05/17

    trace102 = open(
        path_to_data + "/Par" + str(parno) + "/Spectra/G102_trace.reg", "w"
    )
    trace102.write(
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
    )
    trace102.write("wcs;\n")
    # sensitivity drops below 25% of max at wave < 8250 and wave > 11540
    # so box should be 3290 angstroms wide and be centered at 9895.
    trace102.write("box(9895,0,3290,1,1.62844e-12)\n")
    trace102.close()
    trace141 = open(
        path_to_data + "/Par" + str(par) + "/Spectra/G141_trace.reg", "w"
    )
    trace141.write(
        'global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n'
    )
    trace141.write("wcs;\n")
    # sensitivity drops below 25% of max at wave < 10917 and wave > 16904
    # so box should be 5897 angstroms wide and be centered at 13910.5
    trace141.write("box(13910.5,0,5897,1,0)\n")
    trace141.close()

    #     #### STEP 5:  Get zero and first order positions; unpack them ###########
    #     #########################################################################
    #     g102zeroordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G102_0th.reg'
    #     g102firstordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G102_1st.reg'
    #     g141zeroordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G141_0th.reg'
    #     g141firstordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G141_1st.reg'

    #### STEP 6:  Get object information from SExtractor catalog ############
    #########################################################################
    # a_image will give an order of magnitude estimate on the FWHM of the line
    #   this determines the initial guess and sets an upper limit on how broad
    #   it can be.
    # ra/dec, b_image, jmag, jerr, hmag, herr will be carried forward into
    #   the output linelist.
    # find all available cats

    try:
        secats = glob(path_to_data + "/Par" + str(par) + "/DATA/DIRECT_GRISM/Par*phot*.fits")  # MDR 2022/05/17
    except:
        secats = glob(path_to_data + "/Par" + str(par) + "/Products/Par*phot*.fits")  # KVN allowing for different path structure (?)
    secats.sort()
    cat = Table.read(secats[0])

    if verbose == True:
        print("I found the following photometric catalogs...\n")  # MDR 2022/05/17
        print(secats)  # MDR 2022/05/17
        # print('\nThe catalog was read in as...\n') # MDR 2022/05/17
        # print(cat)
        # print('') # MDR 2022/05/17

    # beam = cat['col2'] - M.D.R. - 10/08/2020
    # a_image = cat['col5'] - M.D.R. - 10/08/2020
    # b_image = cat['col6'] - M.D.R. - 10/08/2020
    # ra = cat['col8'] - M.D.R. - 10/08/2020
    # dec = cat['col9'] - M.D.R. - 10/08/2020

    # Edited to call specific columns. - M.D.R. - 10/08/2020
    beam = cat["id"]
    a_image = cat["a_image"]
    b_image = cat["b_image"]
    ra = cat["ra"]
    dec = cat["dec"]

    # Edited to get positions from catalog instead of region files. - K.V.N. 09/28/2023
    #### STEP 5 now comes after STEP 6...
    #### STEP 5:  Get zero and first order positions; unpack them ###########
    #########################################################################
    # g102zeroordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G102_0th.reg'
    # g102firstordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G102_1st.reg'
    # g141zeroordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G141_0th.reg'
    # g141firstordreg = path_to_wisp_data + '/Par' + str(par) + '/DATA/DIRECT_GRISM/G141_1st.reg'

    g102zeroordreg = secats
    g141zeroordreg = secats
    f200zeroordreg = secats

    # which filter is this?
    # if os.path.basename(secats[0]) == 'fin_F110.cat': # MDR 2022/05/17
    #     jmag = cat['col13'] # MDR 2022/05/17
    #     jerr = cat['col14'] # MDR 2022/05/17
    #     # in case there is a F110-only field # MDR 2022/05/17
    #     hmag = np.ones(ra.shape, dtype=float) * 99. # MDR 2022/05/17
    #     herr = np.ones(ra.shape, dtype=float) * 99. # MDR 2022/05/17
    # else:
    jmag = np.ones(ra.shape, dtype=float) * 99.0
    jerr = np.ones(ra.shape, dtype=float) * 99.0
    hmag = cat["mag_auto"]  # hmag = cat['col13'] - M.D.R. - 10/08/2020
    herr = cat["magerr_auto"]  # herr = cat['col14'] - M.D.R. - 10/08/2020
    # read in second file if there are two
    # if len(secats) == 2: # MDR 2022/05/17
    #     cat2 = asciitable.read(secats[1]) # MDR 2022/05/17
    #     # second catalog should be hband # MDR 2022/05/17
    #     hmag = cat2['col13'] # MDR 2022/05/17
    #     herr = cat2['col14'] # MDR 2022/05/17
    objtable = Table(
        [beam, ra, dec, a_image, b_image, jmag, jerr, hmag, herr],
        names=(
            "obj",
            "ra",
            "dec",
            "a_image",
            "b_image",
            "jmag",
            "jerr",
            "hmag",
            "herr",
        ),
    )

    #### STEP 7:  Set up initial ds9 display ################################
    #########################################################################
    if os.path.exists(g102zeroordreg[0]):
        g102zeroarr = getzeroorders_from_cat(g102zeroordreg[0], g="G102")
        show2dNEW(
            "F115W",
            parnos[0],
            objid_unique[0],
            g102zeroarr,
            user,
            "linear",
            path_to_data=path_to_data,
        )
    else:
        g102zeroarr = None
        g102firstarr = None
    if os.path.exists(g141zeroordreg[0]):
        g141zeroarr = getzeroorders_from_cat(g141zeroordreg[0], g="G102")
        show2dNEW(
            "F150W",
            parnos[0],
            objid_unique[0],
            g141zeroarr,
            user,
            "linear",
            path_to_data=path_to_data,
        )
    else:
        g141zeroarr = None
        g141firstarr = None
    if os.path.exists(f200zeroordreg[0]):
        f200zeroarr = getzeroorders_from_cat(f200zeroordreg[0], g="G102")
        show2dNEW(
            "F200W",
            parnos[0],
            objid_unique[0],
            f200zeroarr,
            user,
            "linear",
            path_to_data=path_to_data,
        )
    else:
        f200zeroarr = None
        f200firstarr = None

    # showDirectNEW(
    #     objid_unique[0],
    #     parnos[0],
    #     g102zeroarr,
    #     load_image=True,
    #     path_to_data=path_to_data,
    # )

    # tbaines: show direct images of PASSAGE data
    showDirect_PASSAGE(parno=parnos[0], path_to_data=path_to_data)

    #     if show_dispersed:  # MB
    #         showDispersed(objid_unique[0], parnos[0], load_image=True, path_to_data  = path_to_data)

    #### STEP 8:  Loop through objects ############
    #########################################################################
    if verbose == True:
        print("\nStarting loop through objects...\n")  # MDR 2022/05/17

    remaining_objects = get_remaining_objects(objid_unique, objid_done)
    allobjects = [unique_obj for unique_obj in objid_unique]

    # if (verbose == True):
    #     print('remaining_objects =', remaining_objects) # MDR 2022/05/17
    #     print('allobjects =', allobjects) # MDR 2022/05/17

    print_prompt(
        "\nAs you loop through the objects, you can choose from the following\noptions at any time:\n\txxx = skip to object xxx\n\tb = revisit the previous object\n\tleft = list all remaining objects that need review\n\tlist = list all objects in line list\n\tany other key = continue with the next object\n\th = help/list interactive commands\n\tq = quit\n",
        prompt_type="interim",
    )

    while remaining_objects.shape[0] > 0:
        if path_to_stored_fits == " ":
            use_stored_fits = False
        elif os.path.exists(path_to_stored_fits):
            use_stored_fits = True
            print("looking for stored fit data")
        else:
            use_stored_fits = False
            print("not using stored fit data")

        ndone = len(np.unique(objid_done))
        progress = float(ndone) / float(len(objid_unique)) * 100.0
        print_prompt("\nProgress: %.1f percent" % (progress), prompt_type="interim")

        # do some things as long as there are still objects to inspect
        next_obj = remaining_objects[0]
        print_prompt("Next up: Obj %i" % (next_obj), prompt_type="interim")
        o = input("Enter 'xxx' to skip to Obj xxx or hit any key to continue. > ")

        if o.strip().lower() == "left":
            # remaining_list = ', '.join(['%i'%i for i in remaining_objects])
            print_prompt("Remaining objects:", prompt_type="interim")
            print(remaining_objects)
            o = input("> ")

        if o.strip().lower() == "list":
            print_prompt("All objects:", prompt_type="interim")
            print(allobjects)
            o = input("> ")

        if o.strip().lower() == "b":
            previous_obj = int(objid_done[-1])
            # need to figure out what object came before this one
            # w = np.where(objid_unique == remaining_objects[0])
            # if on first object, this will roll around to previous object
            next_obj = previous_obj
            print_prompt(
                "Going back to previous object: Obj %i" % (next_obj),
                prompt_type="interim",
            )

        elif o.strip().lower() == "q":
            print_prompt(
                "Quitting; saved through previously completed object.",
                prompt_type="interim",
            )
            return 0

        elif o.strip().lower() == "random":
            print_prompt("Hi Vihang!", prompt_type="random")
            next_obj_idx = np.random.randint(remaining_objects.shape[0])
            next_obj = remaining_objects[next_obj_idx]

        elif isFloat(o.strip()):
            next_obj = int(re.search("\d+", o).group())
            # confirm that requested object is in line list
            next_obj = check_input_objid(objid_unique, next_obj, remaining_objects[0])
            next_obj = next_obj if next_obj else remaining_objects[0]
        elif "obj" in o:
            next_obj = int(re.search("\d+", o).group())
            # confirm that requested object is in line list
            next_obj = check_input_objid(objid_unique, next_obj, remaining_objects[0])
            next_obj = next_obj if next_obj else remaining_objects[0]

        # pass the information for this object
        wlinelist = np.where(objid == next_obj)
        lamlines_found = wavelen[wlinelist]
        ston_found = ston[wlinelist]
        wcatalog = np.where(objtable["obj"] == next_obj)
        objinfo = objtable[wcatalog]
        # inspect_object(user, parnos[0], next_obj, objinfo, lamlines_found,
        #               ston_found, g102zeroarr, g141zeroarr, linelistoutfile,
        #               commentsfile, remaining_objects, allobjects,
        #               show_dispersed=show_dispersed)

        if use_stored_fits == True:
            ### get pickle files:
            inpickles = []
            path_pickle1 = (
                path_to_stored_fits
                + "/Par"
                + str(parnos[0])
                + "_output_a/fitdata/Par0_"
                + str(next_obj)
                + "_fitspec.pickle"
            )
            # path_pickle1 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_mbagley/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle2 = path_to_stored_fits + '/Par'  + str(parnos[0]) +    '_output_marc/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle3 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_ben/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle4 = path_to_stored_fits + '/Par'  + str(parnos[0]) +     '_output_claudia/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle5 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_vihang/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle6 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_ivano/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle7 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_mbeck/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle8 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_karlenoid/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle9 = path_to_stored_fits + '/Par'  + str(parnos[0]) +  '_output_mjr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle10 = path_to_stored_fits + '/Par'  + str(parnos[0]) + '_output_sophia/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle11 = '/Volumes/Thunderbay/wisps/mzr_refit/Par'  + str(parnos[0]) + '_output_marc-mzr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'
            # path_pickle12 = '/Volumes/Thunderbay/wisps/mzr_refit/Par'  + str(parnos[0]) + '_output_alaina-mzr/fitdata/Par' + str(parnos[0]) + '_BEAM_' + str(next_obj) + '_fitspec.pickle'

            ### put new fits first
            # if os.path.exists(path_pickle11):
            #     inpickles.append(path_pickle11)
            # if os.path.exists(path_pickle12):
            #     inpickles.append(path_pickle12)
            if os.path.exists(path_pickle1):
                inpickles.append(path_pickle1)
            # if os.path.exists(path_pickle2):
            #     inpickles.append(path_pickle2)
            # if os.path.exists(path_pickle3):
            #     inpickles.append(path_pickle3)
            # if os.path.exists(path_pickle4):
            #     inpickles.append(path_pickle4)
            # if os.path.exists(path_pickle5):
            #     inpickles.append(path_pickle5)
            # if os.path.exists(path_pickle6):
            #     inpickles.append(path_pickle6)
            # if os.path.exists(path_pickle7):
            #     inpickles.append(path_pickle7)
            # if os.path.exists(path_pickle8):
            #     inpickles.append(path_pickle8)
            # if os.path.exists(path_pickle9):
            #     inpickles.append(path_pickle9)
            # if os.path.exists(path_pickle10):
            #     inpickles.append(path_pickle10)

            if len(inpickles) == 0:
                use_stored_fits = False

        if use_stored_fits == True:
            inspect_object(
                user,
                parnos[0],
                next_obj,
                objinfo,
                lamlines_found,
                ston_found,
                g102zeroarr,
                g141zeroarr,
                linelistoutfile,
                commentsfile,
                remaining_objects,
                allobjects,
                show_dispersed=show_dispersed,
                stored_fits=inpickles,
                path_to_data=path_to_data, path_to_code=path_to_code)
        else:
            inspect_object(
                user,
                parnos[0],
                next_obj,
                objinfo,
                lamlines_found,
                ston_found,
                g102zeroarr,
                g141zeroarr,
                linelistoutfile,
                commentsfile,
                remaining_objects,
                allobjects,
                show_dispersed=show_dispersed,
                stored_fits=False,
                path_to_data=path_to_data, path_to_code=path_to_code)
            # if len(glob.glob(path_to_wisp_data +'Par'+ str(parnos[0])+ '/Spectra/Par' +str(parnos[0])+ '_' + str(next_obj).zfill(5)+'*_R.dat')) > 0:
            #     inspect_object(
            #         user,
            #         parnos[0],
            #         next_obj,
            #         objinfo,
            #         lamlines_found,
            #         ston_found,
            #         g102zeroarr,
            #         g141zeroarr,
            #         linelistoutfile,
            #         commentsfile,
            #         remaining_objects,
            #         allobjects,
            #         show_dispersed=show_dispersed,
            #         stored_fits=False,
            #         path_to_wisp_data=path_to_wisp_data, orientation='R')
            # if len(glob.glob(path_to_wisp_data +'Par'+ str(parnos[0])+ '/Spectra/Par' +str(parnos[0])+ '_' + str(next_obj).zfill(5)+'*_C.dat'))> 0:
            #     inspect_object(
            #         user,
            #         parnos[0],
            #         next_obj,
            #         objinfo,
            #         lamlines_found,
            #         ston_found,
            #         g102zeroarr,
            #         g141zeroarr,
            #         linelistoutfile,
            #         commentsfile,
            #         remaining_objects,
            #         allobjects,
            #         show_dispersed=show_dispersed,
            #         stored_fits=False,
            #         path_to_wisp_data=path_to_wisp_data, orientation='C')
                

        objid_done = np.append(objid_done, next_obj)
        remaining_objects = get_remaining_objects(objid_unique, objid_done)

    # outside the while loop, field is done
    redo = " "
    while redo != "q":
        print_prompt(
            "You've finished this field.\nEnter an object ID below to revisit a particular object.\nOtherwise enter 'q' to quit the field.",
            prompt_type="interim",
        )
        redo = input("> ").strip().lower()
        if redo != "q":
            try:
                next_obj = int(re.search("\d+", redo).group())
            except (ValueError, AttributeError):
                print_prompt(
                    "Invalid entry. Enter an object ID or enter 'q' to quit",
                    prompt_type="interim",
                )
            else:
                next_obj = check_input_objid(objid_unique, next_obj, 0)
                if next_obj:
                    # pass the information for this object
                    wlinelist = np.where(objid == next_obj)
                    lamlines_found = wavelen[wlinelist]
                    ston_found = ston[wlinelist]
                    wcatalog = np.where(objtable["obj"] == next_obj)
                    objinfo = objtable[wcatalog]

                    ### get pickle files:
                    ### Note from KVN. field is done so all obj should have pickle files?
                    inpickles = []
                    path_pickle1 = (
                        path_to_stored_fits
                        + "/Par"
                        + str(parnos[0])
                        + "_output_"+str(user)+"/fitdata/PASSAGE_"
                        + str(next_obj)
                        + "_fitspec.pickle")

                    ### put new fits first
                    # if os.path.exists(path_pickle11):
                    #     inpickles.append(path_pickle11)
                    # if os.path.exists(path_pickle12):
                    #     inpickles.append(path_pickle12)
                    if os.path.exists(path_pickle1):
                        inpickles.append(path_pickle1)
                    # if os.path.exists(path_pickle2):
                    #    inpickles.append(path_pickle2)
                    # if os.path.exists(path_pickle3):
                    #    inpickles.append(path_pickle3)
                    # if os.path.exists(path_pickle4):
                    #    inpickles.append(path_pickle4)
                    # if os.path.exists(path_pickle5):
                    #    inpickles.append(path_pickle5)
                    # if os.path.exists(path_pickle6):
                    #    inpickles.append(path_pickle6)
                    # if os.path.exists(path_pickle7):
                    #    inpickles.append(path_pickle7)
                    # if os.path.exists(path_pickle8):
                    #    inpickles.append(path_pickle8)
                    # if os.path.exists(path_pickle9):
                    #    inpickles.append(path_pickle9)
                    # if os.path.exists(path_pickle10):
                    #    inpickles.append(path_pickle10)

                    if len(inpickles) == 0:
                        use_stored_fits = False

                    if use_stored_fits == True:
                        inspect_object(
                            user,
                            parnos[0],
                            next_obj,
                            objinfo,
                            lamlines_found,
                            ston_found,
                            g102zeroarr,
                            g141zeroarr,
                            linelistoutfile,
                            commentsfile,
                            remaining_objects,
                            allobjects,
                            show_dispersed=show_dispersed,
                            stored_fits=inpickles,
                            path_to_data=path_to_data,path_to_code=path_to_code)
                    else:
                        inspect_object(
                            user,
                            parnos[0],
                            next_obj,
                            objinfo,
                            lamlines_found,
                            ston_found,
                            g102zeroarr,
                            g141zeroarr,
                            linelistoutfile,
                            commentsfile,
                            remaining_objects,
                            allobjects,
                            show_dispersed=show_dispersed,
                            stored_fits=False,
                            path_to_data=path_to_data, path_to_code=path_to_code)

                else:
                    break

    make_tarfile(outdir, path_to_code=path_to_code )
    print_prompt(
        "A tarfile of your outputs has been created: %s.tar.gz" % outdir,
        prompt_type="interim",
    )

    # Clean up temp files
    if os.path.exists("./tempcoo.dat") == 1:
        os.unlink("./tempcoo.dat")
    if os.path.exists("./temp_zero_coords.coo") == 1:
        os.unlink("./temp_zero_coords.coo")
    if os.path.exists("./temp110.fits") == 1:
        os.unlink("./temp110.fits")
    if os.path.exists("./temp140.fits") == 1:
        os.unlink("./temp140.fits")
    if os.path.exists("./temp160.fits") == 1:
        os.unlink("./temp160.fits")
    if os.path.exists("./temp_zero_coords.reg") == 1:
        os.unlink("./temp_zero_coords.reg")
    if os.path.exists("G102_trace.reg") == True:
        os.unlink("G102_trace.reg")
    if os.path.exists("G141_trace.reg") == True:
        os.unlink("G141_trace.reg")


# # parnos, objid are scalar not array.
def writeToCatalog(
    catalogname,
    parnos,
    objid,
    ra_obj,
    dec_obj,
    a_image_obj,
    b_image_obj,
    jmag_obj,
    hmag_obj,
    snr_tot_others,
    fitresults,
    contamflags,
    comp_fit):
    if not os.path.exists(catalogname):
        cat = open(catalogname, "w")
        cat.write("#1 objid \n")
        cat.write("#2 redshift \n")
        cat.write("#3 redshift_error \n")
        cat.write("#4 ra_obj \n")
        cat.write("#5 dec_obj \n")
        cat.write("#6 f140w_mag \n")
        cat.write("#7 a_image_obj \n")
        cat.write("#8 b_image_obj \n")
        cat.write("#9 snr_tot_others \n")
        cat.write("#10 chisq \n")
        cat.write("#11 fwhm_muse \n")
        cat.write("#12 fwhm_muse_error \n")
        cat.write("#13 fwhm_g141 \n")
        cat.write("#14 fwhm_g141_error \n")
        cat.write("#15 la_1216_dz \n")
        cat.write("#16 c4_1548_dz \n")
        cat.write("#17 uv_line_dz \n")
        cat.write("#18 m2_2796_dz \n")
        cat.write("#19 o2_3727_dz \n")
        cat.write("#20 o3_5007_dz \n")
        cat.write("#21 s3_he_dz \n")

        result_lines = [
            "la_1216",
            "la_wing",
            "la_1216_wing",
            "n5_1238",
            "n5_1242",
            "n5_1238_1242",
            "c4_1548",
            "c4_1550",
            "c4_1548_1550",
            "h2_1640",
            "o3_1660",
            "o3_1666",
            "o3_1660_1666",
            "s3_1883",
            "s3_1892",
            "s3_1883_1892",
            "c3_1907",
            "c3_1909",
            "c3_1907_1909",
            "m2_2796",
            "m2_2803",
            "m2_2796_2803",
            "o2_3727",
            "o2_3730",
            "o2_3727_3730",
            "hg_4342",
            "o3_4363",
            "h2_4686",
            "hb_4863",
            "o3_4959",
            "o3_5007",
            "o3_4959_5007",
            "o1_6300",
            "o1_6363",
            "o1_6300_6363",
            "n2_6550",
            "ha_6565",
            "n2_6585",
            "ha_6550_6565_6585",
            "s2_6716",
            "s2_6731",
            "s2_6716_6731",
            "s3_9069",
            "s3_9532",
            "s3_9069_9532",
            "he10830",
            "pg_10941",
            "pb_12822",
            "pa_18756",
        ]

        results_idx = 23

        for line in result_lines:
            cat.write("#" + str(results_idx + 0) + " " + line + "_flux \n")
            cat.write("#" + str(results_idx + 1) + " " + line + "_error \n")
            cat.write("#" + str(results_idx + 2) + " " + line + "_ew_obs \n")
            #cat.write("#" + str(results_idx + 3) + " " + line + "_ratio \n")
            cat.write("#" + str(results_idx + 3) + " " + line + "_contam \n")
            results_idx = results_idx + 4

        cat.close()
    # does not leave space before RA?

    # Added KVN 12/2024 - ratio of the broad line is always zero for 1 Gaussian component fit
    ratio = 0

    outstr = (
        "{:<6d}".format(objid)
        + "{:>8.5f}".format(fitresults["redshift"])
        + "{:>8.5f}".format(fitresults["redshift_error"])
        + "{:>12.6f}".format(ra_obj[0])
        + "{:>12.6f}".format(dec_obj[0])
        + "{:>8.2f}".format(hmag_obj[0])
        + "{:>8.3f}".format(a_image_obj[0])
        + "{:>8.3f}".format(b_image_obj[0])
        + "{:>10.2f}".format(snr_tot_others)
        + "{:>10.2f}".format(fitresults["chisq"])
        + "{:>10.2f}".format(fitresults["fwhm_muse"])
        + "{:>10.2f}".format(fitresults["fwhm_muse_error"])
        + "{:>10.2f}".format(fitresults["fwhm_g141"])
        + "{:>10.2f}".format(fitresults["fwhm_g141_error"])
        + "{:>s}".format(str(comp_fit))
        + "{:>10.5f}".format(fitresults["la_1216_dz"])
        + "{:>10.5f}".format(fitresults["c4_1548_dz"])
        + "{:>10.5f}".format(fitresults["uv_line_dz"])
        + "{:>10.5f}".format(fitresults["m2_2796_dz"])
        + "{:>10.5f}".format(fitresults["o2_3727_dz"])
        + "{:>10.5f}".format(fitresults["o3_5007_dz"])
        + "{:>10.5f}".format(fitresults["s3_he_dz"])
        + "{:>13.2e}".format(fitresults["la_1216_flux"])
        + "{:>13.2e}".format(fitresults["la_1216_error"])
        + "{:>13.2e}".format(fitresults["la_1216_ew_obs"])
        + "{:>6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["la_1216"])
        + "{:>13.2e}".format(fitresults["la_wing_flux"])
        + "{:>13.2e}".format(fitresults["la_wing_error"])
        + "{:>13.2e}".format(fitresults["la_wing_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["la_1216"])
        + "{:>13.2e}".format(fitresults["la_1216_wing_flux"])
        + "{:>13.2e}".format(fitresults["la_1216_wing_error"])
        + "{:>13.2e}".format(fitresults["la_1216_wing_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["la_1216"])
        + "{:>13.2e}".format(fitresults["n5_1238_flux"])
        + "{:>13.2e}".format(fitresults["n5_1238_error"])
        + "{:>13.2e}".format(fitresults["n5_1238_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["n5_1238"])
        + "{:>13.2e}".format(fitresults["n5_1242_flux"])
        + "{:>13.2e}".format(fitresults["n5_1242_error"])
        + "{:>13.2e}".format(fitresults["n5_1242_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["n5_1242"])
        + "{:>13.2e}".format(fitresults["n5_1238_1242_flux"])
        + "{:>13.2e}".format(fitresults["n5_1238_1242_error"])
        + "{:>13.2e}".format(fitresults["n5_1238_1242_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["n5_1238"], contamflags["n5_1242"]]))
        + "{:>13.2e}".format(fitresults["c4_1548_flux"])
        + "{:>13.2e}".format(fitresults["c4_1548_error"])
        + "{:>13.2e}".format(fitresults["c4_1548_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["c4_1548"])
        + "{:>13.2e}".format(fitresults["c4_1550_flux"])
        + "{:>13.2e}".format(fitresults["c4_1550_error"])
        + "{:>13.2e}".format(fitresults["c4_1550_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["c4_1550"])
        + "{:>13.2e}".format(fitresults["c4_1548_1550_flux"])
        + "{:>13.2e}".format(fitresults["c4_1548_1550_error"])
        + "{:>13.2e}".format(fitresults["c4_1548_1550_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["c4_1548"], contamflags["c4_1550"]]))
        + "{:>13.2e}".format(fitresults["h2_1640_flux"])
        + "{:>13.2e}".format(fitresults["h2_1640_error"])
        + "{:>13.2e}".format(fitresults["h2_1640_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["h2_1640"])
        + "{:>13.2e}".format(fitresults["o3_1660_flux"])
        + "{:>13.2e}".format(fitresults["o3_1660_error"])
        + "{:>13.2e}".format(fitresults["o3_1660_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o3_1660"])
        + "{:>13.2e}".format(fitresults["o3_1666_flux"])
        + "{:>13.2e}".format(fitresults["o3_1666_error"])
        + "{:>13.2e}".format(fitresults["o3_1666_ew_obs"])        
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o3_1666"])
        + "{:>13.2e}".format(fitresults["o3_1660_1666_flux"])
        + "{:>13.2e}".format(fitresults["o3_1660_1666_error"])
        + "{:>13.2e}".format(fitresults["o3_1660_1666_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["o3_1660"], contamflags["o3_1666"]]))
        + "{:>13.2e}".format(fitresults["s3_1883_flux"])
        + "{:>13.2e}".format(fitresults["s3_1883_error"])
        + "{:>13.2e}".format(fitresults["s3_1883_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["s3_1883"])
        + "{:>13.2e}".format(fitresults["s3_1892_flux"])
        + "{:>13.2e}".format(fitresults["s3_1892_error"])
        + "{:>13.2e}".format(fitresults["s3_1892_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["s3_1892"])
        + "{:>13.2e}".format(fitresults["s3_1883_1892_flux"])
        + "{:>13.2e}".format(fitresults["s3_1883_1892_error"])
        + "{:>13.2e}".format(fitresults["s3_1883_1892_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["s3_1883"], contamflags["s3_1892"]]))
        + "{:>13.2e}".format(fitresults["c3_1907_flux"])
        + "{:>13.2e}".format(fitresults["c3_1907_error"])
        + "{:>13.2e}".format(fitresults["c3_1907_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["c3_1907"])
        + "{:>13.2e}".format(fitresults["c3_1909_flux"])
        + "{:>13.2e}".format(fitresults["c3_1909_error"])
        + "{:>13.2e}".format(fitresults["c3_1909_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["c3_1909"])
        + "{:>13.2e}".format(fitresults["c3_1907_1909_flux"])
        + "{:>13.2e}".format(fitresults["c3_1907_1909_error"])
        + "{:>13.2e}".format(fitresults["c3_1907_1909_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["c3_1907"], contamflags["c3_1909"]]))
        + "{:>13.2e}".format(fitresults["m2_2796_flux"])
        + "{:>13.2e}".format(fitresults["m2_2796_error"])
        + "{:>13.2e}".format(fitresults["m2_2796_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["m2_2796"])
        + "{:>13.2e}".format(fitresults["m2_2803_flux"])
        + "{:>13.2e}".format(fitresults["m2_2803_error"])
        + "{:>13.2e}".format(fitresults["m2_2803_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["m2_2803"])
        + "{:>13.2e}".format(fitresults["m2_2796_2803_flux"])
        + "{:>13.2e}".format(fitresults["m2_2796_2803_error"])
        + "{:>13.2e}".format(fitresults["m2_2796_2803_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["m2_2796"], contamflags["m2_2803"]]))
        + "{:>13.2e}".format(fitresults["o2_3727_flux"])
        + "{:>13.2e}".format(fitresults["o2_3727_error"])
        + "{:>13.2e}".format(fitresults["o2_3727_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o2_3727"])
        + "{:>13.2e}".format(fitresults["o2_3730_flux"])
        + "{:>13.2e}".format(fitresults["o2_3730_error"])
        + "{:>13.2e}".format(fitresults["o2_3730_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o2_3730"])
        + "{:>13.2e}".format(fitresults["o2_3727_3730_flux"])
        + "{:>13.2e}".format(fitresults["o2_3727_3730_error"])
        + "{:>13.2e}".format(fitresults["o2_3727_3730_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["o2_3727"], contamflags["o2_3730"]]))
        + "{:>13.2e}".format(fitresults["hg_4342_flux"])
        + "{:>13.2e}".format(fitresults["hg_4342_error"])
        + "{:>13.2e}".format(fitresults["hg_4342_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["hg_4342"])
        + "{:>13.2e}".format(fitresults["o3_4363_flux"])
        + "{:>13.2e}".format(fitresults["o3_4363_error"])
        + "{:>13.2e}".format(fitresults["o3_4363_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o3_4363"])
        + "{:>13.2e}".format(fitresults["h2_4686_flux"])
        + "{:>13.2e}".format(fitresults["h2_4686_error"])
        + "{:>13.2e}".format(fitresults["h2_4686_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["h2_4686"])
        + "{:>13.2e}".format(fitresults["hb_4863_flux"])
        + "{:>13.2e}".format(fitresults["hb_4863_error"])
        + "{:>13.2e}".format(fitresults["hb_4863_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["hb_4863"])
        + "{:>13.2e}".format(fitresults["o3_4959_flux"])
        + "{:>13.2e}".format(fitresults["o3_4959_error"])
        + "{:>13.2e}".format(fitresults["o3_4959_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o3_4959"])
        + "{:>13.2e}".format(fitresults["o3_5007_flux"])
        + "{:>13.2e}".format(fitresults["o3_5007_error"])
        + "{:>13.2e}".format(fitresults["o3_5007_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o3_5007"])
        + "{:>13.2e}".format(fitresults["o3_4959_5007_flux"])
        + "{:>13.2e}".format(fitresults["o3_4959_5007_error"])
        + "{:>13.2e}".format(fitresults["o3_4959_5007_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["o3_4959"], contamflags["o3_5007"]]))
        + "{:>13.2e}".format(fitresults["o1_6300_flux"])
        + "{:>13.2e}".format(fitresults["o1_6300_error"])
        + "{:>13.2e}".format(fitresults["o1_6300_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o1_6300"])
        + "{:>13.2e}".format(fitresults["o1_6363_flux"])
        + "{:>13.2e}".format(fitresults["o1_6363_error"])
        + "{:>13.2e}".format(fitresults["o1_6363_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["o1_6363"])
        + "{:>13.2e}".format(fitresults["o1_6300_6363_flux"])
        + "{:>13.2e}".format(fitresults["o1_6300_6363_error"])
        + "{:>13.2e}".format(fitresults["o1_6300_6363_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["o1_6300"], contamflags["o1_6363"]]))
        + "{:>13.2e}".format(fitresults["n2_6550_flux"])
        + "{:>13.2e}".format(fitresults["n2_6550_error"])
        + "{:>13.2e}".format(fitresults["n2_6550_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["n2_6550"])
        + "{:>13.2e}".format(fitresults["ha_6565_flux"])
        + "{:>13.2e}".format(fitresults["ha_6565_error"])
        + "{:>13.2e}".format(fitresults["ha_6565_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["ha_6565"])
        + "{:>13.2e}".format(fitresults["n2_6585_flux"])
        + "{:>13.2e}".format(fitresults["n2_6585_error"])
        + "{:>13.2e}".format(fitresults["n2_6585_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["n2_6585"])
        + "{:>13.2e}".format(fitresults["ha_6550_6565_6585_flux"])
        + "{:>13.2e}".format(fitresults["ha_6550_6565_6585_error"])
        + "{:>13.2e}".format(fitresults["ha_6550_6565_6585_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(
            np.max(
                [contamflags["n2_6550"], contamflags["ha_6565"], contamflags["n2_6585"]]
            )
        )
        + "{:>13.2e}".format(fitresults["s2_6716_flux"])
        + "{:>13.2e}".format(fitresults["s2_6716_error"])
        + "{:>13.2e}".format(fitresults["s2_6716_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["s2_6716"])
        + "{:>13.2e}".format(fitresults["s2_6731_flux"])
        + "{:>13.2e}".format(fitresults["s2_6731_error"])
        + "{:>13.2e}".format(fitresults["s2_6731_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["s2_6731"])
        + "{:>13.2e}".format(fitresults["s2_6716_6731_flux"])
        + "{:>13.2e}".format(fitresults["s2_6716_6731_error"])
        + "{:>13.2e}".format(fitresults["s2_6716_6731_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["s2_6716"], contamflags["s2_6731"]]))
        + "{:>13.2e}".format(fitresults["s3_9069_flux"])
        + "{:>13.2e}".format(fitresults["s3_9069_error"])
        + "{:>13.2e}".format(fitresults["s3_9069_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["s3_9069"])
        + "{:>13.2e}".format(fitresults["s3_9532_flux"])
        + "{:>13.2e}".format(fitresults["s3_9532_error"])
        + "{:>13.2e}".format(fitresults["s3_9532_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["s3_9532"])
        + "{:>13.2e}".format(fitresults["s3_9069_9532_flux"])
        + "{:>13.2e}".format(fitresults["s3_9069_9532_error"])
        + "{:>13.2e}".format(fitresults["s3_9069_9532_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(np.max([contamflags["s3_9069"], contamflags["s3_9532"]]))
        + "{:>13.2e}".format(fitresults["he10830_flux"])
        + "{:>13.2e}".format(fitresults["he10830_error"])
        + "{:>13.2e}".format(fitresults["he10830_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["he10830"])
        + "{:>13.2e}".format(fitresults["pg_10941_flux"])
        + "{:>13.2e}".format(fitresults["pg_10941_error"])
        + "{:>13.2e}".format(fitresults["pg_10941_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["pg_10941"])
        + "{:>13.2e}".format(fitresults["pb_12822_flux"])
        + "{:>13.2e}".format(fitresults["pb_12822_error"])
        + "{:>13.2e}".format(fitresults["pb_12822_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["pb_12822"])
        + "{:>13.2e}".format(fitresults["pa_18756_flux"])
        + "{:>13.2e}".format(fitresults["pa_18756_error"])
        + "{:>13.2e}".format(fitresults["pa_18756_ew_obs"])
        + "{:<6d}".format(ratio)   # ratio = 0 ; added KVN 12/24
        + "{:>6d}".format(contamflags["pa_18756"])
        + "\n"
    )

    """
    # if a row already exists for this object, comment it out
    objstr = '{:<8d}'.format(parnos) + '{:<6d}'.format(objid)
    for line in fileinput.input(catalogname, inplace=True):
        if objstr in line:
            print "#%s" % line,
        else:
            print '%s' % line,
#    """

    cat = open(catalogname, "a")
    cat.write(outstr)
    cat.close()

# Added by KVN 21-Aug-2024
def writeToCatalog2gauss(
    catalogname,
    parnos,
    objid,
    ra_obj,
    dec_obj,
    a_image_obj,
    b_image_obj,
    jmag_obj,
    hmag_obj,
    snr_tot_others,
    fitresults,
    contamflags,
    comp_fit):

    if comp_fit == False:
        ratios = 0
    
    if not os.path.exists(catalogname):
        cat = open(catalogname, "w")
        cat.write("#1 objid \n")
        cat.write("#2 redshift \n")
        cat.write("#3 redshift_error \n")
        cat.write("#4 ra_obj \n")
        cat.write("#5 dec_obj \n")
        cat.write("#6 f140w_mag \n")
        cat.write("#7 a_image_obj \n")
        cat.write("#8 b_image_obj \n")
        cat.write("#9 snr_tot_others \n")
        cat.write("#10 chisq \n")
        cat.write("#11 fwhm_muse \n")
        cat.write("#12 fwhm_muse_error \n")
        cat.write("#13 fwhm_g141 \n")
        cat.write("#14 fwhm_g141_error \n")
        cat.write("#15 2Gauss_fit \n")
        cat.write("#16 la_1216_dz \n")
        cat.write("#17 c4_1548_dz \n")
        cat.write("#18 uv_line_dz \n")
        cat.write("#19 m2_2796_dz \n")
        cat.write("#20 o2_3727_dz \n")
        cat.write("#21 o3_5007_dz \n")
        cat.write("#22 s3_he_dz \n")

        
        # flux_strings_2gauss = [
        #     "la_1216_wing",
        #     "n5_1238_1242",
        #     "c4_1548_1550",
        #     "h2_1640tot", "h2_1640nar", "h2_1640bro", 
        #     "o3_1660_1666",
        #     "s3_1883_1892",
        #     "c3_1907_1909",
        #     "m2_2796_2803",
        #     "o2_3727_3730",
        #     "hg_4342tot", "hg_4342nar", "hg_4342bro",
        #     "o3_4363tot", "o3_4363nar", "o3_4363bro",
        #     "h2_4686tot", "h2_4686nar", "h2_4686bro",
        #     "hb_4863tot", "hb_4863nar", "hb_4863bro",
        #     "o3_4959_5007",
        #     "o1_6300_6363",
        #     "ha_6550_6565_6585",
        #     "s2_6716_6731",
        #     "s3_9069",
        #     "s3_9532",
        #     "he10830tot", "he10830nar", "he10830bro",
        #     "pg_10941tot", "pg_10941nar", "pg_10941bro",
        #     "pb_12822tot", "pb_12822nar", "pb_12822bro",
        #     "pa_18756tot", "pa_18756nar", "pa_18756bro",
        # ]
        result_lines = [
            "la_1216",
            "la_wing",
            "la_1216_wing",
            "n5_1238",
            "n5_1242",
            "n5_1238_1242",
            "c4_1548",
            "c4_1550",
            "c4_1548_1550",
            "h2_1640", 
            "o3_1660",
            "o3_1666",
            "o3_1660_1666",
            "s3_1883",
            "s3_1892",
            "s3_1883_1892",
            "c3_1907",
            "c3_1909",
            "c3_1907_1909",
            "m2_2796",
            "m2_2803",
            "m2_2796_2803",
            "o2_3727",
            "o2_3730",
            "o2_3727_3730",
            "hg_4342",
            "o3_4363",
            "h2_4686",
            "hb_4863",
            "o3_4959",
            "o3_5007",
            "o3_4959_5007",
            "o1_6300",
            "o1_6363",
            "o1_6300_6363",
            "n2_6550",
            "ha_6565",
            "n2_6585",
            "ha_6550_6565_6585",
            "s2_6716",
            "s2_6731",
            "s2_6716_6731",
            "s3_9069",
            "s3_9532",
            "s3_9069_9532",
            "he10830",
            "pg_10941",
            "pb_12822",
            "pa_18756"]

        results_idx = 23

        for line in result_lines:
            cat.write("#" + str(results_idx + 0) + " " + line + "_flux \n")
            cat.write("#" + str(results_idx + 1) + " " + line + "_error \n")
            cat.write("#" + str(results_idx + 2) + " " + line + "_ew_obs \n")
            cat.write("#" + str(results_idx + 3) + " " + line + "_ratio \n")
            cat.write("#" + str(results_idx + 4) + " " + line + "_contam \n")
            results_idx = results_idx + 5

        cat.close()
    # does not leave space before RA?

    outstr = (
        "{:<6d}".format(objid)
        + "{:>8.5f}".format(fitresults["redshift"])
        + "{:>8.5f}".format(fitresults["redshift_error"])
        + "{:>12.6f}".format(ra_obj[0])
        + "{:>12.6f}".format(dec_obj[0])
        + "{:>8.2f}".format(hmag_obj[0])
        + "{:>8.3f}".format(a_image_obj[0])
        + "{:>8.3f}".format(b_image_obj[0])
        + "{:>10.2f}".format(snr_tot_others)
        + "{:>10.2f}".format(fitresults["chisq"])
        + "{:>10.2f}".format(fitresults["fwhm_muse"])
        + "{:>10.2f}".format(fitresults["fwhm_muse_error"])
        + "{:>10.2f}".format(fitresults["fwhm_g141"])
        + "{:>10.2f}".format(fitresults["fwhm_g141_error"])
        + "{:>s}".format(str(comp_fit))
        + "{:>10.5f}".format(fitresults["la_1216_dz"])
        + "{:>10.5f}".format(fitresults["c4_1548_dz"])
        + "{:>10.5f}".format(fitresults["uv_line_dz"])
        + "{:>10.5f}".format(fitresults["m2_2796_dz"])
        + "{:>10.5f}".format(fitresults["o2_3727_dz"])
        + "{:>10.5f}".format(fitresults["o3_5007_dz"])
        + "{:>10.5f}".format(fitresults["s3_he_dz"])
        + "{:>13.2e}".format(fitresults["la_1216_flux"])
        + "{:>13.2e}".format(fitresults["la_1216_error"])
        + "{:>13.2e}".format(fitresults["la_1216_ew_obs"])
        + "{:>6d}".format(contamflags["la_1216"])
        + "{:>13.2e}".format(fitresults["la_wing_flux"])
        + "{:>13.2e}".format(fitresults["la_wing_error"])
        + "{:>13.2e}".format(fitresults["la_wing_ew_obs"])
        + "{:>6d}".format(contamflags["la_1216"])
        + "{:>13.2e}".format(fitresults["la_1216_wing_flux"])
        + "{:>13.2e}".format(fitresults["la_1216_wing_error"])
        + "{:>13.2e}".format(fitresults["la_1216_wing_ew_obs"])
        + "{:>6d}".format(contamflags["la_1216"])
        + "{:>13.2e}".format(fitresults["n5_1238_flux"])
        + "{:>13.2e}".format(fitresults["n5_1238_error"])
        + "{:>13.2e}".format(fitresults["n5_1238_ew_obs"])
        + "{:>6d}".format(contamflags["n5_1238"])
        + "{:>13.2e}".format(fitresults["n5_1242_flux"])
        + "{:>13.2e}".format(fitresults["n5_1242_error"])
        + "{:>13.2e}".format(fitresults["n5_1242_ew_obs"])
        + "{:>6d}".format(contamflags["n5_1242"])
        + "{:>13.2e}".format(fitresults["n5_1238_1242_flux"])
        + "{:>13.2e}".format(fitresults["n5_1238_1242_error"])
        + "{:>13.2e}".format(fitresults["n5_1238_1242_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["n5_1238"], contamflags["n5_1242"]]))
        + "{:>13.2e}".format(fitresults["c4_1548_flux"])
        + "{:>13.2e}".format(fitresults["c4_1548_error"])
        + "{:>13.2e}".format(fitresults["c4_1548_ew_obs"])
        + "{:>6d}".format(contamflags["c4_1548"])
        + "{:>13.2e}".format(fitresults["c4_1550_flux"])
        + "{:>13.2e}".format(fitresults["c4_1550_error"])
        + "{:>13.2e}".format(fitresults["c4_1550_ew_obs"])
        + "{:>6d}".format(contamflags["c4_1550"])
        + "{:>13.2e}".format(fitresults["c4_1548_1550_flux"])
        + "{:>13.2e}".format(fitresults["c4_1548_1550_error"])
        + "{:>13.2e}".format(fitresults["c4_1548_1550_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["c4_1548"], contamflags["c4_1550"]]))
        + "{:>13.2e}".format(fitresults["h2_1640tot_flux"])
        + "{:>13.2e}".format(fitresults["h2_1640tot_error"])
        + "{:>13.2e}".format(fitresults["h2_1640tot_ew_obs"])
        + "{:>6d}".format(contamflags["h2_1640"])
        + "{:>13.2e}".format(fitresults["o3_1660_flux"])
        + "{:>13.2e}".format(fitresults["o3_1660_error"])
        + "{:>13.2e}".format(fitresults["o3_1660_ew_obs"])
        + "{:>6d}".format(contamflags["o3_1660"])
        + "{:>13.2e}".format(fitresults["o3_1666_flux"])
        + "{:>13.2e}".format(fitresults["o3_1666_error"])
        + "{:>13.2e}".format(fitresults["o3_1666_ew_obs"])
        + "{:>6d}".format(contamflags["o3_1666"])
        + "{:>13.2e}".format(fitresults["o3_1660_1666_flux"])
        + "{:>13.2e}".format(fitresults["o3_1660_1666_error"])
        + "{:>13.2e}".format(fitresults["o3_1660_1666_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["o3_1660"], contamflags["o3_1666"]]))
        + "{:>13.2e}".format(fitresults["s3_1883_flux"])
        + "{:>13.2e}".format(fitresults["s3_1883_error"])
        + "{:>13.2e}".format(fitresults["s3_1883_ew_obs"])
        + "{:>6d}".format(contamflags["s3_1883"])
        + "{:>13.2e}".format(fitresults["s3_1892_flux"])
        + "{:>13.2e}".format(fitresults["s3_1892_error"])
        + "{:>13.2e}".format(fitresults["s3_1892_ew_obs"])
        + "{:>6d}".format(contamflags["s3_1892"])
        + "{:>13.2e}".format(fitresults["s3_1883_1892_flux"])
        + "{:>13.2e}".format(fitresults["s3_1883_1892_error"])
        + "{:>13.2e}".format(fitresults["s3_1883_1892_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["s3_1883"], contamflags["s3_1892"]]))
        + "{:>13.2e}".format(fitresults["c3_1907_flux"])
        + "{:>13.2e}".format(fitresults["c3_1907_error"])
        + "{:>13.2e}".format(fitresults["c3_1907_ew_obs"])
        + "{:>6d}".format(contamflags["c3_1907"])
        + "{:>13.2e}".format(fitresults["c3_1909_flux"])
        + "{:>13.2e}".format(fitresults["c3_1909_error"])
        + "{:>13.2e}".format(fitresults["c3_1909_ew_obs"])
        + "{:>6d}".format(contamflags["c3_1909"])
        + "{:>13.2e}".format(fitresults["c3_1907_1909_flux"])
        + "{:>13.2e}".format(fitresults["c3_1907_1909_error"])
        + "{:>13.2e}".format(fitresults["c3_1907_1909_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["c3_1907"], contamflags["c3_1909"]]))
        + "{:>13.2e}".format(fitresults["m2_2796_flux"])
        + "{:>13.2e}".format(fitresults["m2_2796_error"])
        + "{:>13.2e}".format(fitresults["m2_2796_ew_obs"])
        + "{:>6d}".format(contamflags["m2_2796"])
        + "{:>13.2e}".format(fitresults["m2_2803_flux"])
        + "{:>13.2e}".format(fitresults["m2_2803_error"])
        + "{:>13.2e}".format(fitresults["m2_2803_ew_obs"])
        + "{:>6d}".format(contamflags["m2_2803"])
        + "{:>13.2e}".format(fitresults["m2_2796_2803_flux"])
        + "{:>13.2e}".format(fitresults["m2_2796_2803_error"])
        + "{:>13.2e}".format(fitresults["m2_2796_2803_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["m2_2796"], contamflags["m2_2803"]]))
        + "{:>13.2e}".format(fitresults["o2_3727_flux"])
        + "{:>13.2e}".format(fitresults["o2_3727_error"])
        + "{:>13.2e}".format(fitresults["o2_3727_ew_obs"])
        + "{:>6d}".format(contamflags["o2_3727"])
        + "{:>13.2e}".format(fitresults["o2_3730_flux"])
        + "{:>13.2e}".format(fitresults["o2_3730_error"])
        + "{:>13.2e}".format(fitresults["o2_3730_ew_obs"])
        + "{:>6d}".format(contamflags["o2_3730"])
        + "{:>13.2e}".format(fitresults["o2_3727_3730_flux"])
        + "{:>13.2e}".format(fitresults["o2_3727_3730_error"])
        + "{:>13.2e}".format(fitresults["o2_3727_3730_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["o2_3727"], contamflags["o2_3730"]]))
        + "{:>13.2e}".format(fitresults["hg_4342tot_flux"])
        + "{:>13.2e}".format(fitresults["hg_4342tot_error"])
        + "{:>13.2e}".format(fitresults["hg_4342tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["hg_4342bro_flux"]/fitresults["hg_4342tot_flux"] if not np.isnan(fitresults["hg_4342bro_flux"]/fitresults["hg_4342tot_flux"]) else 0) # ratio of lines; added KVN 12/2024 
        + "{:>6d}".format(contamflags["hg_4342"])
        + "{:>13.2e}".format(fitresults["o3_4363tot_flux"])
        + "{:>13.2e}".format(fitresults["o3_4363tot_error"])
        + "{:>13.2e}".format(fitresults["o3_4363tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["o3_4363bro_flux"]/fitresults["o3_4363tot_flux"] if not np.isnan(fitresults["o3_4363bro_flux"]/fitresults["o3_4363tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["o3_4363"])
        + "{:>13.2e}".format(fitresults["h2_4686tot_flux"])
        + "{:>13.2e}".format(fitresults["h2_4686tot_error"])
        + "{:>13.2e}".format(fitresults["h2_4686tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["h2_4686bro_flux"]/fitresults["h2_4686tot_flux"] if not np.isnan(fitresults["h2_4686bro_flux"]/fitresults["h2_4686tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["h2_4686"])
        + "{:>13.2e}".format(fitresults["hb_4863tot_flux"])
        + "{:>13.2e}".format(fitresults["hb_4863tot_error"])
        + "{:>13.2e}".format(fitresults["hb_4863tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["hb_4863bro_flux"]/fitresults["hb_4863tot_flux"] if not np.isnan(fitresults["hb_4863bro_flux"]/fitresults["hb_4863tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["hb_4863"])
        + "{:>13.2e}".format(fitresults["o3_4959_flux"])
        + "{:>13.2e}".format(fitresults["o3_4959_error"])
        + "{:>13.2e}".format(fitresults["o3_4959_ew_obs"])
        + "{:>6d}".format(contamflags["o3_4959"])
        + "{:>13.2e}".format(fitresults["o3_5007_flux"])
        + "{:>13.2e}".format(fitresults["o3_5007_error"])
        + "{:>13.2e}".format(fitresults["o3_5007_ew_obs"])
        + "{:>6d}".format(contamflags["o3_5007"])
        + "{:>13.2e}".format(fitresults["o3_4959_5007_flux"])
        + "{:>13.2e}".format(fitresults["o3_4959_5007_error"])
        + "{:>13.2e}".format(fitresults["o3_4959_5007_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["o3_4959"], contamflags["o3_5007"]]))
        + "{:>13.2e}".format(fitresults["o1_6300_flux"])
        + "{:>13.2e}".format(fitresults["o1_6300_error"])
        + "{:>13.2e}".format(fitresults["o1_6300_ew_obs"])
        + "{:>6d}".format(contamflags["o1_6300"])
        + "{:>13.2e}".format(fitresults["o1_6363_flux"])
        + "{:>13.2e}".format(fitresults["o1_6363_error"])
        + "{:>13.2e}".format(fitresults["o1_6363_ew_obs"])
        + "{:>6d}".format(contamflags["o1_6363"])
        + "{:>13.2e}".format(fitresults["o1_6300_6363_flux"])
        + "{:>13.2e}".format(fitresults["o1_6300_6363_error"])
        + "{:>13.2e}".format(fitresults["o1_6300_6363_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["o1_6300"], contamflags["o1_6363"]]))
        + "{:>13.2e}".format(fitresults["n2_6550tot_flux"])
        + "{:>13.2e}".format(fitresults["n2_6550tot_error"])
        + "{:>13.2e}".format(fitresults["n2_6550tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["n2_6550bro_flux"]/fitresults["n2_6550tot_flux"] if not np.isnan(fitresults["n2_6550bro_flux"]/fitresults["n2_6550tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["n2_6550"])
        + "{:>13.2e}".format(fitresults["ha_6565tot_flux"])
        + "{:>13.2e}".format(fitresults["ha_6565tot_error"])
        + "{:>13.2e}".format(fitresults["ha_6565tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["ha_6565bro_flux"]/fitresults["ha_6565tot_flux"] if not np.isnan(fitresults["ha_6565bro_flux"]/fitresults["ha_6565tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["ha_6565"])
        + "{:>13.2e}".format(fitresults["n2_6585tot_flux"])
        + "{:>13.2e}".format(fitresults["n2_6585tot_error"])
        + "{:>13.2e}".format(fitresults["n2_6585tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["n2_6585bro_flux"]/fitresults["n2_6585tot_flux"] if not np.isnan(fitresults["n2_6585bro_flux"]/fitresults["n2_6585tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["n2_6585"])
        + "{:>13.2e}".format(fitresults["ha_6550_6565_6585_flux"])
        + "{:>13.2e}".format(fitresults["ha_6550_6565_6585_error"])
        + "{:>13.2e}".format(fitresults["ha_6550_6565_6585_ew_obs"])
        + "{:>6d}".format(
            np.max(
                [contamflags["n2_6550"], contamflags["ha_6565"], contamflags["n2_6585"]]
            )
        )
        + "{:>13.2e}".format(fitresults["s2_6716_flux"])
        + "{:>13.2e}".format(fitresults["s2_6716_error"])
        + "{:>13.2e}".format(fitresults["s2_6716_ew_obs"])
        + "{:>6d}".format(contamflags["s2_6716"])
        + "{:>13.2e}".format(fitresults["s2_6731_flux"])
        + "{:>13.2e}".format(fitresults["s2_6731_error"])
        + "{:>13.2e}".format(fitresults["s2_6731_ew_obs"])
        + "{:>6d}".format(contamflags["s2_6731"])
        + "{:>13.2e}".format(fitresults["s2_6716_6731_flux"])
        + "{:>13.2e}".format(fitresults["s2_6716_6731_error"])
        + "{:>13.2e}".format(fitresults["s2_6716_6731_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["s2_6716"], contamflags["s2_6731"]]))
        + "{:>13.2e}".format(fitresults["s3_9069_flux"])
        + "{:>13.2e}".format(fitresults["s3_9069_error"])
        + "{:>13.2e}".format(fitresults["s3_9069_ew_obs"])
        + "{:>6d}".format(contamflags["s3_9069"])
        + "{:>13.2e}".format(fitresults["s3_9532_flux"])
        + "{:>13.2e}".format(fitresults["s3_9532_error"])
        + "{:>13.2e}".format(fitresults["s3_9532_ew_obs"])
        + "{:>6d}".format(contamflags["s3_9532"])
        + "{:>13.2e}".format(fitresults["s3_9069_9532_flux"])
        + "{:>13.2e}".format(fitresults["s3_9069_9532_error"])
        + "{:>13.2e}".format(fitresults["s3_9069_9532_ew_obs"])
        + "{:>6d}".format(np.max([contamflags["s3_9069"], contamflags["s3_9532"]]))
        + "{:>13.2e}".format(fitresults["he10830tot_flux"])
        + "{:>13.2e}".format(fitresults["he10830tot_error"])
        + "{:>13.2e}".format(fitresults["he10830tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["he10830bro_flux"]/fitresults["he10830tot_flux"] if not np.isnan(fitresults["he10830bro_flux"]/fitresults["he10830tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["he10830"])
        + "{:>13.2e}".format(fitresults["pg_10941tot_flux"])
        + "{:>13.2e}".format(fitresults["pg_10941tot_error"])
        + "{:>13.2e}".format(fitresults["pg_10941tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["pg_10941bro_flux"]/fitresults["pg_10941tot_flux"] if not np.isnan(fitresults["pg_10941bro_flux"]/fitresults["pg_10941tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["pg_10941"])
        + "{:>13.2e}".format(fitresults["pb_12822tot_flux"])
        + "{:>13.2e}".format(fitresults["pb_12822tot_error"])
        + "{:>13.2e}".format(fitresults["pb_12822tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["pb_12822bro_flux"]/fitresults["pb_12822tot_flux"] if not np.isnan(fitresults["pb_12822bro_flux"]/fitresults["pb_12822tot_flux"]) else 0) # ratio of lines; added KVN 12/2024 
        + "{:>6d}".format(contamflags["pb_12822"])
        + "{:>13.2e}".format(fitresults["pa_18756tot_flux"])
        + "{:>13.2e}".format(fitresults["pa_18756tot_error"])
        + "{:>13.2e}".format(fitresults["pa_18756tot_ew_obs"])
        + "{:>13.2e}".format(fitresults["pa_18756bro_flux"]/fitresults["pa_18756tot_flux"] if not np.isnan(fitresults["pa_18756bro_flux"]/fitresults["pa_18756tot_flux"]) else 0) # ratio of lines; added KVN 12/2024
        + "{:>6d}".format(contamflags["pa_18756"])
        + "\n"
    )

    """
    # if a row already exists for this object, comment it out
    objstr = '{:<8d}'.format(parnos) + '{:<6d}'.format(objid)
    for line in fileinput.input(catalogname, inplace=True):
        if objstr in line:
            print "#%s" % line,
        else:
            print '%s' % line,
#    """

    print('What is this in the catalogs? ', "{:>13.2e}".format(fitresults["he10830tot_flux"] / fitresults["he10830bro_flux"] if not np.isnan(fitresults["he10830tot_flux"] / fitresults["he10830bro_flux"]) else 0)) # ratio of lines; added KVN 12/2024)
    cat = open(catalogname, "a")
    cat.write(outstr)
    cat.close()

def writeFitdata(filename, lam, flux, eflux, contam, zero, fit, continuum, masks):
    if verbose == True:
        print("\nRunning writeFitdata...\n")  # MDR 2022/05/17
    """ """
    fitspec_file = filename + ".dat"
    t = Table(
        [lam, flux, eflux, contam, zero, fit, continuum, masks],
        names=(
            "Lam",
            "Flam",
            "Flam_err",
            "Contam",
            "Zero",
            "Fitmodel",
            "Contmodel",
            "Masked",
        ),
    )
    t["Lam"].format = "{:.1f}"
    t["Flam"].format = "{:.5e}"
    t["Flam_err"].format = "{:.5e}"
    t["Contam"].format = "{:.5e}"
    t["Zero"].format = "{:.0f}"
    t["Fitmodel"].format = "{:.5e}"
    t["Contmodel"].format = "{:.5e}"
    t["Masked"].format = "{:.0f}"
    asciitable.write(
        t,
        fitspec_file,
        fill_values=[(asciitable.masked, "--")],
        overwrite=True,
        format="fixed_width_two_line",
        position_char="#",
        delimiter_pad=" ")


def writeComments(filename, parnos, objid, comment):
    if verbose == True:
        print("Running writeComments...\n")  # MDR 2022/05/17
    if os.path.exists(filename) == False:
        cat = open(filename, "w")
    else:
        cat = open(filename, "a")

    outstr = f"Par: {parnos} ID: {objid:<6d}" + comment + "\n"

    cat.write(outstr)
    cat.close()


def UpdateCatalog(linelistoutfile):
    if verbose == True:
        print("Running UpdateCatalog...\n")  # MDR 2022/05/17
    allDirectoryFiles = os.listdir("./fitdata/")
    objid_list = []
    for obj in allDirectoryFiles:
        x = obj.split("_")[2]
        objid_list.append(int(x))
        Parno = obj.split("_")[0]  # this is inefficient, but I don't care.
    objid_list = np.sort(np.unique(np.array(objid_list)))
    for obj in objid_list:
        print(obj)
        inpickle = "./fitdata/" + Parno + "_BEAM_" + str(obj) + "_fitspec.pickle"
        fileObject = open(inpickle, "r")
        alldata = pickle.load(fileObject)
        # definition from above
        #                      0          1                 2      3        4           5            6         7         8           9         10
        # output_meta_data = [parnos[0], objid_unique[i], ra_obj, dec_obj, a_image_obj, b_image_obj, jmag_obj, hmag_obj, fitresults, flagcont, config_pars]
        parnos = alldata[0]
        objid_unique = alldata[1]
        ra_obj = alldata[2]
        dec_obj = alldata[3]
        a_image_obj = alldata[4]
        b_image_obj = alldata[5]
        jmag_obj = alldata[6]
        hmag_obj = alldata[7]
        fitresults = alldata[8]
        flagcont = alldata[9]
        # config_pars = alldata[10] ## not used here.

        if verbose == True:
            print("Writing to catalog...\n")  # MDR 2022/05/17

        if comp_fit == True:
            writeToCatalog2gauss(
                linelistoutfile,
                parnos,
                objid_unique,
                ra_obj,
                dec_obj,
                a_image_obj,
                b_image_obj,
                jmag_obj,
                hmag_obj,
                fitresults,
                flagcont,
                comp_fit)

        if comp_fit == False:
            WriteToCatalog(
                linelistoutfile,
                parnos,
                objid_unique,
                ra_obj,
                dec_obj,
                a_image_obj,
                b_image_obj,
                jmag_obj,
                hmag_obj,
                fitresults,
                flagcont, 
                comp_fit)







