import numpy as np
from netCDF4 import Dataset
from datetime import date, datetime, timedelta
import os
import shutil
import urllib.request
from contextlib import closing
from math import sqrt, sin, cos, pi, floor, isnan
import matplotlib
matplotlib.use('Agg')
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from flask import Flask, request, jsonify, send_file 

monthNames = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
ftpFolder = 'ftp://ftp.awi.de/sea_ice/product/cryosat2_smos/v206/nh/'

app = Flask(__name__)
@app.route("/")
def hello():
	return "Generate a cryosat-smos sea ice thickness map for a given date"

@app.get("/plot-thickness")
def plotThicknessMap():
	return plotDate(False)
	
@app.get("/plot-anomaly")
def plotAnomaly():
	return plotDate(True)
	
def plotDate(isAnomaly):
	dateString = request.args.get('date')
	if dateString == None:
		return "Please insert valid date (YYYY-MM-DD)", 400
	print(dateString)
	
	try:
		date = datetime.strptime(dateString, '%Y-%m-%d')
	except ValueError:
		return "Invalid  date. It needs to have the format YYYY-MM-DD"
	if date < datetime(2010,11,4):
		return "Date cannot be earlier than 2010-11-04", 400
	if (date.month, date.day) > (4,12) and (date.month, date.day) < (10,18):
		return "Date must lie between 15 October and 15 April", 400
	
	plotType = 'anomaly' if isAnomaly else 'thickness'
	storedFileName = 'data/' + plotType + '/' + str(date.year) + '/cryosat-smos-thickness' + ('-anomaly-' if isAnomaly else '-') + dateString.replace('-', '') + '.png'
	if os.path.isfile(storedFileName):
		copiedFileName = 'tmp/storedfile.png'
		shutil.copyfile(storedFileName, copiedFileName)
		return send_file(copiedFileName, mimetype='image/png')
	
	try:
		downloadedFileName = download(date)
	except:
		if date.year > 2023 or date.year == 2023 and date.month > 6:
			downloadedFileName = download(date, True)
		else:
			raise
	
	lat = np.loadtxt(open("lat.csv", "rb"), delimiter=",", skiprows=0)
	lon = np.loadtxt(open("lon.csv", "rb"), delimiter=",", skiprows=0)
	dummyvalue=10
	
	thicknessmax = 5.0 # 5 meters
	anomalymax = 1.0 # 1 meter
	anomyears = 10 # 10 years in anomaly base

	griddedThickness = getGriddedThickness(downloadedFileName)
	
	dayOfYear = date.timetuple().tm_yday
		
	if not isAnomaly:
		multiplier = 1
		landmask = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, date.year, dummyvalue, lat, lon)
		landmask = interpolate(landmask, dummyvalue, False, thicknessmax, anomalymax)

		plotTitle = "CryoSat-SMOS sea ice thickness " + str(date.day) + " " + monthNames[date.month-1] + " " + str(date.year)
		filename = 'tmp/cryosat-smos-thickness.png'
		plotThickness(landmask, plotTitle, filename, thicknessmax)
		
	else:
		startyear = 2014 if date.month <= 4 else 2013
		averageFileName = './data/avg/cryosat-smos-avg-' + str(startyear) + '-to-' + str(startyear + 9) + '-' + padzeros(date.month) + padzeros(date.day) + '.csv'
		
		multiplier = -1.0*0.001 #/anomyears 
		landmask = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, date.year, dummyvalue, lat, lon)
		maskbis = np.loadtxt(averageFileName, delimiter = ',', skiprows=0)
		landmask = addMasks(landmask, maskbis, multiplier, dummyvalue)
		print('plotting cryosat anomaly', date)
		# for k in range(anomyears):
			# compyear = 2023-k-1 # year-k-1
			# filename = download(datetime(compyear, date.month, date.day))
			# print('comp file',filename)
			# griddedThickness = getGriddedThickness(filename)
			# maskbis = insertCryosatDataInNsidcMask(griddedThickness, dayOfYear, compyear, dummyvalue, lat, lon)
			# landmask = addMasks(landmask, maskbis, multiplier, dummyvalue)

		landmask = interpolate(landmask, dummyvalue, True, thicknessmax, anomalymax)

		print(date.day)
		plotTitle = "CryoSat-SMOS thickness anomaly " + str(date.day) + " " + monthNames[date.month-1] + " " + str(date.year) + " vs 2013-2022"
		filename = 'tmp/cryosat-smos-thickness-anomaly.png'
		plotAnomaly(landmask, plotTitle, filename, anomalymax)
	
	return send_file(filename, mimetype='image/png')
	
def padzeros(n):
	"""
	Left pad a number with zeros. 
    """
	return str(n) if n >= 10 else '0'+str(n)

def getFileName(date, useRevised):
	startDate = date - timedelta(days = 3)
	endDate = date + timedelta(days = 3)
	return 'W_XX-ESA,SMOS_CS2,NH_25KM_EASE2_' + str(startDate.year) + padzeros(startDate.month) + padzeros(startDate.day)+ '_' + str(endDate.year) + padzeros(endDate.month) + padzeros(endDate.day) + '_' + ('r' if useRevised or date.year < 2023 or date.year == 2023 and date.month < 6 else 'o') + '_v206_01_l4sit.nc'
	
def getGriddedThickness(filename):
	f = Dataset(filename, 'r', format="NETCDF4")
	thicknessData = f.variables['analysis_sea_ice_thickness'][:]
	f.close()
	return thicknessData
	
def download(date, useRevised = False):
	"""
	Download Cryosat-SMOS ftp file. 
    """
	filename = getFileName(date, useRevised)
	ftpSubfolder = (str(date.year) + "/" + padzeros(date.month)) if (useRevised or date.year < 2023 or date.year == 2023 and date.month < 6) else 'LATEST'
	fullFtpPath = ftpFolder + ftpSubfolder + "/" + filename.replace(',','%2C')
	localpath = 'tmp/gridded-file.nc'
	print('downloading file ', fullFtpPath)
	with closing(urllib.request.urlopen(fullFtpPath)) as r:
		with open(localpath, 'wb') as f:
			shutil.copyfileobj(r, f)
	print('downloaded file ', fullFtpPath)
	return localpath	

def getNsidcLandMask():
	landmask = np.genfromtxt('landmask_nsidc.csv', delimiter=',')	
	print('nsidc land mask',landmask.shape)
	n = landmask.shape[0]  # n is assumed to be an odd number
	n2 = 1
	landmask = landmask[n2:-n2,n2:-n2]
	n = landmask.shape[0]  # n is assumed to be an odd number
	print('nsidc land mask cropped',landmask.shape)
	global landmaskcenter
	landmaskcenter = (n-1)/2
	return landmask
	
def insertCryosatDataInNsidcMask(cryosatData, day, year, dummyvalue, lat, lon):
	
	landmask = getNsidcLandMask()
	landmask = landmask*dummyvalue
	landmaskSize = landmask.shape[0]
	counter=1000.0
	rand = (year+day-2000)/2000000.0 #random.randrange(1,100000)/10000000.0
	print('random', rand)
	
	for i in range(0,432):
		for j in range(0,432):
			v = cryosatData[0,i,j]
			latitude = lat[i,j]
			longitude = lon[i,j]
			rad = 360*sqrt(2)*sin(pi*(90-latitude)/360)
			y = int(round(landmaskcenter+rad*sin(pi*longitude/180.0)))
			x = int(round(landmaskcenter+rad*cos(pi*longitude/180.0)))
			
			if(x < 0 or y < 0 or x >= landmaskSize or y >= landmaskSize or landmask[x,y] == 0): #or mask[i,j] == -1
				continue				
			if not type(v) is np.float64:
				#if isnan(refmask[i,j]):
				#	continue
				v = 0
				#continue
			if(landmask[x,y] == dummyvalue):
				landmask[x,y] = 0
			landmask[x,y] += counter + max(v, rand)

	for x in range(0,landmask.shape[0]):
		for y in range(0,landmask.shape[1]):
			if(landmask[x,y] == 0 or landmask[x,y] == dummyvalue):
				continue
			n = floor(landmask[x,y]/counter)
			landmask[x,y] = (landmask[x,y] - counter*n)/n
			
	return landmask

def getInterpolatedValue(x, y, landmask, dummyvalue):
	radius = 1
	while(radius < 10):
		for k in range(radius):
			other = landmask[x+k,y+radius-k] if x+k < 359 and y+radius-k < 359 else 0
			if(other != dummyvalue and other != 0):
				return other
			other = landmask[x+radius-k,y-k] if x+radius-k < 359 and y-k >= 0 else 0
			if(other != dummyvalue and other != 0):
				return other
			other = landmask[x-k,y-radius+k] if x-k >= 0 and y-radius+k >= 0 else 0
			if(other != dummyvalue and other != 0):
				return other
			other = landmask[x-radius+k,y+k] if x-radius+k >= 0 and y+k < 359 else 0
			if(other != dummyvalue and other != 0):
				return other
		radius += 1
	return dummyvalue

def interpolate(landmask, dummyvalue, anomalyplot, thicknessmax, anomalymax):
	mask = landmask.copy()
	
	for x in range(0,landmask.shape[0]):
		for y in range(0,landmask.shape[1]):
			if(mask[x,y] == 0): # land
				if(anomalyplot):
					mask[x,y] = -anomalymax
				continue
			if(mask[x,y] == dummyvalue): # value to be interpolated
				mask[x,y] = getInterpolatedValue(x, y, landmask, dummyvalue)
			if((anomalyplot and abs(mask[x,y]) < 0.001) or ((not anomalyplot) and mask[x,y] < 0.05)): # hide in maps
				mask[x,y] = dummyvalue
			if(anomalyplot and mask[x,y] != dummyvalue and mask[x,y] > anomalymax * 0.99):
				mask[x,y] = anomalymax * 0.99
			if(anomalyplot and mask[x,y] < -anomalymax * 0.99):
				mask[x,y] = -anomalymax * 0.99
			if((not anomalyplot) and mask[x,y] != dummyvalue and mask[x,y] > thicknessmax - 0.06):
				mask[x,y] = thicknessmax - 0.06				
					
	return mask

def plotThickness(landmask, plotTitle, filename, thicknessmax):
	cdict = {'red':   ((0.0,  0.5, 0.5),
					   (0.001, 0.5, 0.0),
		           	   (0.05, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.15, 0.0, 0.0),					   
					   (0.2, 0.0, 0.2),
				   	   (0.25, 0.2, 0.4),					  
					   (0.3, 0.4, 0.6),
				   	   (0.35, 0.6, 0.8),
					   (0.4, 0.8, 1.0),
					   (0.45, 1.0, 1.0),
				   	   (0.5, 1.0, 1.0),
					   (0.55, 1.0, 1.0),
					   (0.6, 1.0, 0.95),
					   (0.65, 0.95, 0.9),
					   (0.7, 0.9, 0.85),
					   (0.75, 0.85, 0.8),
					   (0.8, 0.8, 0.75),
					   (0.85, 0.75, 0.7),
					   (0.9, 0.7, 0.65),
					   (0.95, 0.65, 0.6),
				       (0.999,  0.6, 1),
                       (1.0,  1, 1)),

         'green':      ((0.0,  0.5, 0.5),
         	           (0.001, 0.5, 0.0),
         	           (0.05, 0.0, 0.1),
					   (0.1, 0.1, 0.25),
					   (0.15, 0.25, 0.4),
					   (0.2, 0.4, 0.55),
					   (0.25, 0.55, 0.7),
					   (0.3, 0.7, 0.85),
					   (0.35, 0.85, 1.0),
					   (0.4, 1.0, 1.0),		
					   (0.45, 1.0, 0.9),				   					   
				   	   (0.5, 0.9, 0.8),
					   (0.55, 0.8, 0.75),
					   (0.6, 0.75, 0.7),
					   (0.65, 0.7, 0.6),
					   (0.7, 0.6, 0.5),
					   (0.75, 0.5, 0.4),
					   (0.8, 0.4, 0.3),
					   (0.85, 0.3, 0.2),
					   (0.9, 0.2, 0.1),
					   (0.95, 0.1, 0.0),
         	           (0.999,  0.0, 1),
                       (1.0,  1, 1)),

         'blue':       ((0.0,  0.5, 0.5),
         	           (0.001, 0.4, 0.4),
         	           (0.05, 0.4, 0.55),
					   (0.1, 0.55, 0.7),
					   (0.15, 0.7, 0.85),
				       (0.2, 0.85, 1.0),	
					   (0.25, 1.0, 0.8),						     
				   	   (0.3, 0.8, 0.6),		
					   (0.35, 0.6, 0.4),		   	   
					   (0.4, 0.4, 0.2),		
					   (0.45, 0.2, 0.0),		   				   
					   (0.5, 0.0, 0.0),
					   (0.55, 0.0, 0.0),
					   (0.6, 0.0, 0.0),
					   (0.7, 0.0, 0.0),
					   (0.8, 0.0, 0.0),
					   (0.9, 0.0, 0.0),
         	           (0.999,  0.0, 0.0),
                       (1.0,  1, 1))}
	kleur = LinearSegmentedColormap('BlueRed1', cdict)
	#plt.register_cmap(cmap=kleur)
	# next part only serves to get a nicer colormap in the plot
	cbrol = {'red':   ((0.0,  0.0, 0.0),
					   (0.05, 0.0, 0.0),
					   (0.1, 0.0, 0.0),
					   (0.15, 0.0, 0.0),					   
					   (0.2, 0.0, 0.2),
				   	   (0.25, 0.2, 0.4),					  
					   (0.3, 0.4, 0.6),
				   	   (0.35, 0.6, 0.8),
					   (0.4, 0.8, 1.0),
					   (0.45, 1.0, 1.0),
				   	   (0.5, 1.0, 1.0),
					   (0.55, 1.0, 1.0),
					   (0.6, 1.0, 0.95),
					   (0.65, 0.95, 0.9),
					   (0.7, 0.9, 0.85),
					   (0.75, 0.85, 0.8),
					   (0.8, 0.8, 0.75),
					   (0.85, 0.75, 0.7),
					   (0.9, 0.7, 0.65),
					   (0.95, 0.65, 0.6),		
                       (1.0,  0.6, 0.6)),


			 'green': ((0.0,  0.0, 0.0),
					   (0.05, 0.0, 0.1),
					   (0.1, 0.1, 0.25),
					   (0.15, 0.25, 0.4),
					   (0.2, 0.4, 0.55),
					   (0.25, 0.55, 0.7),
					   (0.3, 0.7, 0.85),
					   (0.35, 0.85, 1.0),
					   (0.4, 1.0, 1.0),		
					   (0.45, 1.0, 0.9),				   					   
				   	   (0.5, 0.9, 0.8),
					   (0.55, 0.8, 0.75),
					   (0.6, 0.75, 0.7),
					   (0.65, 0.7, 0.6),
					   (0.7, 0.6, 0.5),
					   (0.75, 0.5, 0.4),
					   (0.8, 0.4, 0.3),
					   (0.85, 0.3, 0.2),
					   (0.9, 0.2, 0.1),
					   (0.95, 0.1, 0.0),
					   (1.0,  0.0, 0.0)),

			 'blue':  ((0.0, 0.4, 0.4),
					   (0.05, 0.4, 0.55),
					   (0.1, 0.55, 0.7),
					   (0.15, 0.7, 0.85),
				       (0.2, 0.85, 1.0),	
					   (0.25, 1.0, 0.8),						     
				   	   (0.3, 0.8, 0.6),		
					   (0.35, 0.6, 0.4),		   	   
					   (0.4, 0.4, 0.2),		
					   (0.45, 0.2, 0.0),		   				   
					   (0.5, 0.0, 0.0),
					   (0.55, 0.0, 0.0),
					   (0.6, 0.0, 0.0),
					   (0.7, 0.0, 0.0),
					   (0.8, 0.0, 0.0),
					   (0.9, 0.0, 0.0),
					   (1.0,  0.0, 0.0))}
	kleurbrol = LinearSegmentedColormap('BlueRed2', cbrol)
	#plt.register_cmap(cmap=kleurbrol)
	mask = landmask[50:-90,80:-90]#landmask[85:-100,95:-100]#landmask[30:-70,10:-70]
	n = landmask.shape[0]
	try:
		plt.colorbar().remove()
	except:
		print('error remove color bar thickness')
	plt.clf()
	plt.cla()
	figbrol = plt.imshow(mask, extent=(0,n,0,n), vmin= 0, vmax=thicknessmax,
			   interpolation='nearest', cmap=kleurbrol)

	#plot the relevant map:
	fig2 = plt.imshow(mask, extent=(0,n,0,n), vmin= 0, vmax=thicknessmax,
           interpolation='nearest', cmap=kleur)
	
	plt.title(plotTitle)
	cb = plt.colorbar(figbrol)
	plt.xticks([])
	plt.yticks([])
	cb.set_label("meters")
	#plt.show()
	plt.savefig(filename)
		
def plotAnomaly(landmask, plotTitle, filename, anomalymax):
	cdict = {'red': ((0.0,  0.4, 0.4),
         	       (0.001, 0.0, 0.0),
         	       #(0.4, 0.8, 0.8),
         	       (0.5, 1.0, 1.0),
         	       (0.999,  0.0, 0.0),
                   (1.0,  1, 1)),

         'green':   ((0.0,  0.4, 0.4),
		           (0.001, 0.0, 0.0),
		           (0.5, 1.0, 1.0),
				   (0.999,  1.0, 1.0),
                   (1.0,  1, 1)),

         'blue':  ((0.0,  0.4, 0.4),
         	       (0.001, 0.4, 0.4),
         	       #(0.4, 1, 0),
         	       (0.5, 1, 0.5),
         	       (0.999,  0.0, 0.0),
                   (1.0,  1, 1))}
	kleur = LinearSegmentedColormap('BlueRed3', cdict)
	#plt.register_cmap(cmap=kleur)
	# next part only serves to get a nicer colormap in the plot
	cbrol = {'red':   ((0.0,  0.0, 0.0),
					   (0.5, 1, 1),
					   (1.0,  0.0, 0.0)),
			 'green': ((0.0,  0, 0),
					   (0.5, 1, 1),
					   (1.0,  1, 1)),

			 'blue':  ((0.0, 0.4, 0.4),
					   #(0.4, 1, 0),
					   (0.5, 1, 0.5),
					   (1.0,  0.0, 0.0))}
	kleurbrol = LinearSegmentedColormap('BlueRed4', cbrol)
	#plt.register_cmap(cmap=kleurbrol)
	mask = landmask[50:-90,80:-90]#landmask[85:-100,95:-100]#landmask[30:-70,10:-70]
	n = mask.shape[0]
	m = mask.shape[1]
	try:
		plt.colorbar().remove()
	except:
		print('error remove color bar anomaly')
	plt.clf()
	plt.cla()
	figbrol = plt.imshow(mask, extent=(0,n,0,n), vmin= -anomalymax, vmax=anomalymax,
			   interpolation='nearest', cmap=kleurbrol)

	#plot the relevant map:
	fig2 = plt.imshow(mask, extent=(0,n,0,n), vmin= -anomalymax, vmax=anomalymax,
           interpolation='nearest', cmap=kleur)
	
	plt.title(plotTitle)
	cb = plt.colorbar(figbrol)
	plt.xticks([])
	plt.yticks([])
	cb.set_label("meters")
	#plt.show()
	plt.savefig(filename)
	
def addMasks(landmask, mask, multiplier, dummyvalue):
	for x in range(0,landmask.shape[0]):
		for y in range(0,landmask.shape[1]):
			if(landmask[x,y] == 0 or landmask[x,y] == dummyvalue):
				continue
			#if(mask[x,y] == 0 or mask[x,y] == dummyvalue):
			#	raise ValueError(str(mask[x,y]) + 'Mask adding not allowed')
			landmask[x,y] = landmask[x,y] + multiplier*mask[x,y]

	return landmask