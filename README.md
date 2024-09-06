This package merges line fitting code from AH with the interactive redshift 
catalog from NR.  Several helper functions can be found in the following modules.  
The key code is "measure_z_interactive", which operates on a candidate line list produced by
"detect_grism_lines.py."   

The key improvements over "get_z_interactive" are the following: 


1. The plots and ds9 views are displayed in the same way as previously in "get_z_interactive." 
However, now, the line and continuum fits are shown as well.   

2. The interactive options have changed to allow the user to get a good fit. new options include: 

-- masking regions of the spectrum
-- change the g102-g141 transition wavelength
-- modify the continuum spline by changing the nodes 
-- trying a different "guess" fwhm 
-- or entering an arbitrary redshift. 

When the fit is satisfactory, the user types "a" (accept) or "ac" (accept contaminatd), and all the flux measurements
for all possible lines are appended to a catalog.  

3. The catalog is no longer a line list, but rather one line per object.  It has many columns, but is 
readable with astropy.io.ascii.  The new catalog is called ParNNNlines_catalog.dat.  

4. Along the same strategy, we are carrying out the interactive fitting on an object-by-object basis, rather than
line-by-line.  Because all the lines are fit simultaneously, this is the only sensible thing to do. 

5. The user can still choose to leave comments, and these are stored separately since they are not machine readable. 
The file contains an entry for each object in the candidate list, whether comments are entered or not.  This serves as 
a way to mark our progress. 

6. When we stop the code and start again, the comments file is read to find the last object that was completed. 

7. Both the comments and catalog file are appended, rather than deleted and re-written in their entirty every time
an object is completed.  This makes the business of saving our temporary work, as was done with "get_z_interactive" no longer necessary.   Every time we 
accept/reject an object, the relevant files are opened, written, and closed.  

8.  The code saves fit figures in Spectra/figs. These are presently ugly. 
(The catalogs are still stored in Spectra/linelists).    The bottom panel (s/n vs lam) in get_z_interactive was 
judged to be low-information and removed.  It is not the job of a person to determine if the line has sufficient s/n 
to be included or not. 

9.  Some variables which should not be hard coded (i.e. the G102-G141 transition wavelength) are stored in a configuration file.  
At the moment, this file needs to be in the directory where the code is being run.

10. To run: 

--unpack the code and put the wisp directory in your python path 
I have the directory wisp_analysis in /Users/ahenry/ahenry_py 

--then my  bash file: 

export PYTHONPATH="/Users/ahenry/ahenry_py:$PYTHONPATH"
(if you have not defined $PYTHONPATH previously, then simply remove that part of the line.)

--- go to the Spectra subdirectory, after "detect_grism_lines.py" has been run and the "linelist" directory exists. 

--- copy the default.config file to this directory 

--- ipython 

--- import wisp_analysis as wisp 

--- wisp.measure_z_interactive() 

have fun! 


11.  the plots are ugly... this is what they're showing:

black = spectrum

magenta  = zero orders 

red = contam 

red  = model  
(These are really easy to tell apart even with the same color, but I should still change it.) 

blue-dashed = continuum model 

red vertical = line from candidate list; if multiple lines in the candidate list exist, then the one with the highest preliminary s/n is chosen and guessed to be ha. 

blue vertical = expected location of other lines; list taken from get_z_interactive 
lam_Oii,lam_Hbeta,lam_Oiii_1,lam_Oiii_2,lam_Halpha,lam_Sii,lam_Siii_1,lam_Siii_2,lam_He,lam_Pag,lam_Fe,lam_Pab]
note that the last three are not yet included in the model fit, but they can be added. 


cyan dotted vertical = G102/G141 transition. this line will move if you change the transition wavelength in the interactive fitter. 


black dots: spline nodes. The default for these is set in default.config, and they can also be changed in the interactiev fitter 

blue dots: locations of all peaks found by "detect_grism_lines.py" 





