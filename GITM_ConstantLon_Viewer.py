# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import spacepy
import gitm 
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
import cmocean
import glob, os
from cartopy import config
import cartopy.crs as ccrs
from matplotlib.backends.backend_pdf import PdfPages
# Pull in the rc from matplotlib to make pretty fonts (LaTeX Fonts)
from matplotlib import rc
import matplotlib as mpl

# First, we find all 3DALL binary files in the directory and then
# we simply access the first one.

#------------------------------------------------------
# Set some global parameters for Text
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
rc('text', usetex=True)
mpl.rcParams["lines.linewidth"] = 1.5
mpl.rcParams["contour.negative_linestyle"]='dashed'
# Turn off Interactive MODE (Suppresses the Figures when writing files)
plt.ioff()


#-----------------------------------------------------
# ACCESS GITM DATA FILES
# FILE 1
print('===> Searching Working Directory for GITM 3DALL Files. \n')
filelist = []  # initialize a filelist array
for file in glob.glob("*.bin"):
    filelist.append(file)
    
if len(filelist) > 0:
   print('===> Identified ', len(filelist), ' Gitm Binary Files.\n')
   testfile = []
   for file in filelist:
       print(file)
       testfile.append(file)
else:
    print("!!! Couldn't find any 3DALL*.bin GITM files !!!!")
    raise Exception()
    
gdata = gitm.GitmBin(testfile[0])

PI = np.pi

OrignLons = gdata.attrs['nLon']
OrignLats = gdata.attrs['nLat']
OrignAlts = gdata.attrs['nAlt']

# Pull out actual data
nGCs = 2
alts = gdata['Altitude' ][0     ,0     ,nGCs:OrignAlts-nGCs-1]/1000.0      # Alts in km
lats = gdata['Latitude' ][0     ,nGCs:OrignLats-nGCs-1,0     ]*(180.0/PI)  # Lats in deg
lons = gdata['Longitude'][nGCs:OrignLons-nGCs-1,0     ,0     ]*(180.0/PI)  # Lons in deg

nLons = len(lons)
nLats = len(lats)
nAlts = len(alts)

# These are the indices of the Altitudes that we want to pull out.
LonString = list()

#-----------------------------------------------------
# BEGIN OUTPUT OF PDF FILES
#-----------------------------------------------------
axis_text_size = 18
title_text_size = 28
tick_font_size = 13
subpanel_text_size = 14

pdfile = 'ConstantLongitude_TemperatureContourComparisons.pdf'
with PdfPages(pdfile) as pdf:
    #for i in AltitudeSlices:
    for i in range(nLons):
        lonpick = lons[i]
        newstring = str(lonpick)
        splitstring = newstring.split('.')
        LonString.append(splitstring[0])

    LonSliceIndex = 0
#    for iAlt in AltitudeSlices:
    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
        Temperature = gdata['Temperature'   ][iLon,:nLats,:nAlts]   # Pull out Temperature
#        U           = gdata['V!Dn!N (east)' ][iLon,:nLats,:nAlts]
#        V           = gdata['V!Dn!N (north)'][iLon,:nLats,:nAlts]
#        W           = gdata['V!Dn!N (up)'][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= Temperature[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.thermal)
        # Set axis text size and the labels
        TitleString = 'Temperature (K) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')
        
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)

        
        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        #fig.savefig('GITM_Sample_Contour_Plasma.eps')
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)
LonSliceIndex = 0

pdfile = 'ConstantLongitude_ZonalWindContourComparisons.pdf'
with PdfPages(pdfile) as pdf:

    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
#        Temperature = gdata['Temperature'   ][iLon,:nLats,:nAlts]   # Pull out Temperature
        U           = gdata['V!Dn!N (east)' ][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= U[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        #fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Zonal Winds (m/s) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')
        
        
        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)


LonSliceIndex = 0

pdfile = 'ConstantLongitude_MeridionalWindContourComparisons.pdf'
with PdfPages(pdfile) as pdf:

    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
        V           = gdata['V!Dn!N (north)'][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= V[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        #fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Meridional Winds (m/s) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')
        
        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        #fig.savefig('GITM_Sample_Contour_Plasma.eps')
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)


LonSliceIndex = 0

pdfile = 'ConstantLongitude_VerticalWindContourComparisons.pdf'
with PdfPages(pdfile) as pdf:

    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
        W           = gdata['V!Dn!N (up)'][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= W[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        #fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Verical Winds (m/s) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')

        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        #fig.savefig('GITM_Sample_Contour_Plasma.eps')
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)




LonSliceIndex = 0
pdfile = 'ConstantLongitude_IonZonalWindContourComparisons.pdf'
with PdfPages(pdfile) as pdf:

    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
#        Temperature = gdata['Temperature'   ][iLon,:nLats,:nAlts]   # Pull out Temperature
        U           = gdata['V!Di!N (east)' ][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= U[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        #fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Zonal Ion Winds (m/s) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')
        
        
        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)


LonSliceIndex = 0

pdfile = 'ConstantLongitude_IonMeridionalWindContourComparisons.pdf'
with PdfPages(pdfile) as pdf:

    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
        V           = gdata['V!Di!N (north)'][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= V[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        #fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Meridional Ion Winds (m/s) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')
        
        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        #fig.savefig('GITM_Sample_Contour_Plasma.eps')
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)


LonSliceIndex = 0

pdfile = 'ConstantLongitude_IonVerticalWindContourComparisons.pdf'
with PdfPages(pdfile) as pdf:

    for iLon in range(nLons):
        print('Accessing iLon %i',iLon)
        iLon = min(iLon,nLons-1)
        W           = gdata['V!Di!N (up)'][iLon,:nLats,:nAlts]
        # Create local numpy arrays for plotting
        X = np.linspace(0,1,len(lats))
        Y = np.linspace(0,1,len(alts))
        for i in range(0,len(lats)):
            X[i] = lats[i]
        for j in range(0,len(alts)):
            Y[j] = alts[j]
        
        # Normally, I would use the Meshgrid function here, but it creates the wrong
        # array size (want two 2-D Arrays X2D, Y2D that are (lons,lats) in size)
        array_shape = (len(lats), len(alts))
        X2D = np.zeros(array_shape,order='F')
        Y2D = np.zeros(array_shape,order='F')
        Z2D = np.zeros(array_shape,order='F')
        for i in range(0,nLats):
            for j in range(0,nAlts):
                X2D[i,j] = X[i]
                Y2D[i,j] = Y[j]
        
        for i in range(0,len(lats)):
            for j in range (0, len(alts)):
                Z2D[i,j]= W[i,j]

        projection = ccrs.PlateCarree()       
        # Figsize = (width, height) in inches
        # subplots (nrows, ncolumns, figure_settings)
        # --
        # Instantiate a figure container.  Fig is the whole figure
        # ax is each sub figure.
        #fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,11),sharex=True)
        fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(8.5,5.0),sharex=True)
        bw1 = ax.contour(X2D,Y2D,Z2D,10,colors='k',linewidths=1.5)
        cpf = ax.contourf(X2D,Y2D,Z2D,20,cmap=cmocean.cm.balance)
        # Set axis text size and the labels
        TitleString = 'Verical Ion Winds (m/s) at Lon = ' + LonString[LonSliceIndex] + ' deg'
        ax.annotate(TitleString,xy=(0.1,0.95),xycoords='figure fraction',fontsize=title_text_size)
        ax.tick_params(axis='both',labelsize=tick_font_size)
        ax.set_ylabel('Altitude (km)',size = axis_text_size)
        ax.set_xlabel('Latitude (deg)',size = axis_text_size)
        
        # Tell the colorbar where to exist
        # fig.colorbar(subplotname, axis name)
        cbar = fig.colorbar(cpf,ax=ax)   
        cbar.ax.tick_params(labelsize='12')

        plt.subplots_adjust(hspace=0.1, right = 1.05)
        pdf.savefig()
        plt.close()
        #fig.savefig('GITM_Sample_Contour_Plasma.eps')
        LonSliceIndex = LonSliceIndex + 1

print('Finished Output to File: '+ pdfile)
