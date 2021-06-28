# GITM_Python_Plotting_SampleCodes

Some Sample Plotting Codes (both 1-D) and 3-D.

1.  GITM_1D_Plotter_Multiple.py:  This reads in GITM binary files (1-D), accesses the types of variables to plot and then asks the user to select from the available keys.  This will then produce plots until the user wants to stop.

  This code has some key features:  (1) accessing and getting inputs from the user based upon the variables in the GitmBin files.  (2) using logarithmic / non-log plotting, (3) saving files to a .pdf file.
  
  2.  GITM_ConstantAlt_Viewer.py/GITM_ConstantLon_Viewer.py:  These use 3D files to slice the binary gitm files and plot contours of key variables.  These codes are fully automated and will simply plot a contour for each altitude or longitude, and they don't ask the user for any inputs.  These are probably only useful for seeing how to slice the GITM data and how I have use contour plots in the past.

3. Sample TimeSeriesPlot.eps:  This is for Gargi to see how we have produced time series data from the 1-D model in the past.  The bottom panel is the contour plot with time on the horizontal and altitude on the vertical.  The contours are Temperatures and top panels are merely constant altitude slices to show the user how the temperatures are changing over time.  The two lines in the top panel correspond to what is happening in the contour plot if you were to slice along the horizontal solid and dashed lines  Don't worry about trying to recreate this exact type of plot, but merely get a type of time-series contour going in python.

