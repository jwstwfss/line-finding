# Documentation for line-finding code

This software is used to identify line-emitting objects and measure emission line properties in JWST NIRISS WFSS Grism spectra, based on the pure-parallel survey PASSAGE (PID#1571).

After installation, to run the software, simply go into the git cloned directory and run 

*python mainPASSAGE.py*

The user will be asked to enter the number of a parallel field and a username of choice. The fitting is then performed on an object-by-object basis. 

The following is a list of commands for this software.

**OBJECT SPECIFIC OPTIONS:**  

a = accept object fit  

ac = accept object fit, noting contamination  

r = reject object  

c = add comment  

user = toggle between previously saved fits  

contam = specify contamination to line flux and/or continuum  

reset = reset interactive options back to default for this object  

s = print the (in progress) object summary


**EMISSION LINE SPECIFIC OPTIONS:**  

z = enter a different z guess  

w = enter a different emission line wavelength guess

dz = change the allowable redshift difference between lines  

n = skip to next brightest line found in this object

2gauss = double gaussian profile for the line being fitted

1gauss = option to go back to 1 gaussian fit after selecting 2 gaussian fit

ha, hb, hg, o31, o32, o2, s2, s31, s32, lya, c4, pb, pa, pg = change strongest emission line

The full list of commands and corresponding lines are as follows

| **Command** | **Line**       | **Vacuum Wavelength (Å)** |
| ----------- | -------------- | ------------------------- |
| lya         | Ly-alpha 1215  | 1215.670                  |
| c4          | CIV 1548       | 1548.203                  |
| o2          | [OII] 3730     | 3729.875                  |
| hg          | H-gamma 4342   | 4341.684                  |
| hb          | H-beta 4863    | 4862.683                  |
| o31         | [OIII] 4959    | 4960.295                  |
| o32         | [OIII] 5007    | 5008.240                  |
| ha          | H-alpha 6563   | 6564.610                  |
| s2          | [SII] 6716     | 6718.290                  |
| s31         | [SIII] 9069    | 9071.100                  |
| s32         | [SIII] 9532    | 9533.200                  |
| he          | HeI 10830      | 10832.86                  |
| pg          | Pa-gamma 10941 | 10941.1                   |
| pb          | Pa-beta 12822  | 12821.6                   |
| pa          | Pa-alpha 18756 | 18756.1                   |


**SPECTRUM SPECIFIC OPTIONS:**  

fw = change the fwhm guess in pixels  

t1, t2 = change transition wavelength between F115W and F150W (t1) and F150W and F200W (t2)  

m1, m2, or m3 = mask up to three discontinuous wavelength regions  

nodes = change the wavelengths for the continuum spline nodes  

addnodes = add wavelengths for the continuum spline nodes

rmnodes = remove wavelengths from the continuum spline nodes

shiftallnodes = SHIFT ALL nodes used for the continuum spline by some wavelength   

bluecut = change the blue cutoff of the F115W grism  

redcut  = change the red cutoff of the F200W grism

lincont = fit continuum as a line

polycont = fit continuum as a higher-order polynomial

splinecont = fit continuum as a spline (piecewise) polynomial

grismr = use only Grism-R spectrum for line-fitting

grismrcontam = use only Grism-R spectrum (with contamination) for line-fitting

grismc = use only Grism-C spectrum for line-fitting

grismccontam = use only Grism-C spectrum (with contamination) for line-fitting

comb = Use combined spectrum (default)

combcontam = Use combined spectrum with contamination



**DS9 SPECIFIC OPTIONS:**  

lin = linear z-scale  

log = logarithmic z-scale

zs102 = z1,z2 comma-separated range for G102 z-scale  

zs141 = z1,z2 comma-separated range for G141 z-scale  

dc = recenter direct images  

reload = reload direct images  

dr = reload direct image reg files

**SOFTWARE SPECIFIC OPTIONS:**  

h = print this message  
q = quit
