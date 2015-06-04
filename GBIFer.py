#!/bin/python
'''
GBIFer Copyright (C) 2014 Jose Sergio Hleap

This a program that given a GBIF zip dataset(s) will extract the taxonomic and geographic info
to be analyzed. This program supercedes biogeographer, since the GBIF API change

TODO: Include GBIF API to do it directly
      Develop a list base exclude to exclude more than one query
	  Develop statistics
	  Include method to add csv files

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
Author: Jose Sergio Hleap
email: jshleap@dal.ca
'''

# importing bit###################################################################################
import zipfile, datetime
import numpy as np
import cPickle as pickl
import matplotlib.pyplot as plt
from glob import glob as G
from math import ceil,floor
from mpl_toolkits.basemap import Basemap, shiftgrid, cm

# END importing bit###############################################################################

# Info bit #######################################################################################
__author__   = 'Jose Sergio Hleap'
__email__    = 'jshleap (at) squalus.org'
__created__  = datetime.datetime(2014, 12, 15)
__modified__ = datetime.datetime(2015, 01, 21)
__version__  = "0.1"
__status__   = "Development"
# End Info bit ###################################################################################

# Functions bit###################################################################################

def isfloat(value):
	'''stupid isdigit for floats'''
	try:
		float(value)
		return True
	except:
		return False

class record:
	'''Class to hold individual GBIF records'''
	def __init__(self,prefix,line):
		'''
		This function initializes the record class.
		
		:param self: intance of the class record
		:type self: :class `record`
		:param prefix: Prefix for any of the outputs for the class
		:type prefix: string
		:param line: A comma separated line containing GBIF occurrences as in occurrence.txt
		:type line: string
		'''
		self.prefix = prefix
		self.parse_gbif_line(line)

	def parse_gbif_line(self,line):
		'''
		This function parse individual gbif records given a comma-delimited string, 
		populating the class.
		
		:param self: intance of the class record
		:type self: :class `record`
		:param line: A comma separated line containing GBIF occurrences as in occurrence.txt
		:type line: string
		'''
		codes = pickl.load(open('countrycodes.pickl','rb'))
		bl = line.strip().split('\t')
		self.gbifID = bl[0]
		self.identifier = bl[26]
		self.basisOfRecord = bl[62]
		self.catalogNumber=bl[65]
		self.phylum=bl[163]
		self.taxclass=bl[66]
		self.family=bl[92]
		self.genus=bl[99]
		self.order=bl[156]
		self.species=bl[218]
		self.scientificName=bl[172]
		self.taxonRank=bl[181]
		self.collectionCode=bl[67]
		self.collectionID=bl[68]
		self.countryCode=bl[70]
		self.country=codes[bl[70]]
		self.Lat=float(bl[77])
		self.Long=float(bl[78])
		self.institutionCode=bl[123]
		self.occurrenceID=bl[153]
		self.remarks=bl[154]
		if isfloat(bl[201]):
			self.depth = float(bl[201])
		else: 
			self.depth = bl[201]

	def write_record2CSV(self):
		'''
		This function will write (append) a single record to csv file.
		
		:param self: intance of the class record
		:type self: :class `record`
		'''
		title=['Species','Genus','Family','Order','Class','Phylum','Latitude',
		       'Longitude','country']
		title=','.join(title)
		line = [self.species,self.genus,self.family,self.order,self.taxclass,
		        self.phylum,self.Lat,self.Long,self.country]
		with open(self.prefix+'.csv', mode='a+') as F:
			F.seek(0)
			first_char = F.read(1)
			if not first_char:
				F.write(title+'\n')
			F.write(','.join([str(l) for l in line])+'\n')


class GBIFer:
	'''class to manipulate zip files downloaded from GBIF'''
	def __init__(self,prefix='gbif', filename='',exclude=None):
		'''
		This function initializes the GBFer class.
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param prefix: Prefix for any of the outputs for the class. (gbif by default)
		:type prefix: string
		:param filename: Name of a giving zip file downloaded from the GBIF database. If
		not provided, it will iterate over all zipfiles in the current working directory
		:type filename: string
		:param exclude: genus, species or both that you'd like to exclude from analyses.
		:type exclude: None or string
		'''		
		self.prefix = prefix
		self.filename = filename
		self.specciesData={}
		self.latglobmin=100
		self.longlobmin=300
		self.latglobmax=-100
		self.longlobmax=-300
		self.exclude=exclude

		if self.filename == '':
			zips = G("*.zip")
			self.get_zips(zips)
		else:
			self.get_zips([self.filename])

	def get_zips(self,lists):
		'''
		Iterate over files and extract occurrences
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param lists: List of string with the names of the zipfiles with the data
		:type lists: list
		'''
		for z in lists:
			zf = zipfile.ZipFile(z)
			data = zf.read('occurrence.txt')
			self.parse_occurrences(data)

	def get_globals(self,lat,lon):
		'''
		Populate the global latitude and longitude
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param lat: The latitude value of a particular occurrence
		:type lat: float
		:param lon: The longitude value of a particular occurrence
		:type lon: float
		'''
		lat=float(lat); lon=float(lon)
		if lat < self.latglobmin:
			self.latglobmin = float(lat)		
		elif lat > self.latglobmax:
			self.latglobmax = float(lat)

		if lon < self.longlobmin:
			self.longlobmin = lon		
		elif lon > self.longlobmax:
			self.longlobmax = float(lon)

	def parse_occurrences(self,filecontent):
		'''
		Parse GBIF occurrence file, and populate globals and speciesData
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param filecontent: The unzipped content of an occurrence file
		:type filecontent: string
		'''
		F = filecontent.split('\n')
		for line in F:
			if line.strip().startswith("gbifID") or line == '':
				continue
			else:
				r = record(self.prefix,line)
				r.write_record2CSV()
				sp = r.species
				lat= r.Lat
				lon= r.Long
				self.get_globals(lat, lon)
				if sp not in self.specciesData:
					self.specciesData[sp]=[(lat,lon)]
				else:
					self.specciesData[sp].append((lat,lon))

	def point_in_square(self, up,lo,point):
		'''
		Helper function to check if an occurrence is on the looking grid
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param up: Upper corner of the square/rectangle in (x,y)/(lat,lon)
		:type up: tuple
		:param lo: Lower corner of the square/rectangle in (x,y)/(lat,lon)
		:type lo: tuple
		:param point: Point to be checked
		:type point: tuple
		:return: boolean
		'''
		x = point[0]; y = point[1]; xu = up[0]; yu = up[1]; xl = lo[0]; yl = lo[1]

		if (x >= xl and x <= xu) and (y >= yl and y <= yu):
			return True
		else:
			return False

	def sp_per_unit(self,squp,sqlo,exclude=None):
		'''
		Count the number of species per grid cell.
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param squp: Upper corner of the square/rectangle in (x,y)/(lat,lon)
		:type squp: tuple
		:param sqlo: Lower corner of the square/rectangle in (x,y)/(lat,lon)
		:type sqlo: tuple
		'''
		spcount=0
		sps = []
		for k,v in self.specciesData.iteritems():
			if k == '':
				continue
			for e in v:
				if self.point_in_square(squp, sqlo, e) and (k not in sps) and (str(exclude) not in k):
					spcount+=1
					sps.append(k)
		return spcount

	def richness_grid_by_point(self,unit=5):
		'''
		Compute the richness in the grid given the maximum and minimum occurrences of the species
		(as opposed to by point). This will affect the counts if antitropical species are found.
		It most likely over estimate the counts.
			
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param unit: The bin unit for the grid in degrees.
		:type unit: integer
		:returns: A tuple with 3 numpy arrays (richness, latitudes, longitudes)
		'''
		rich=[]
		lats = np.arange(ceil(self.latglobmin),floor(self.latglobmax),unit)
		lons = np.arange(ceil(self.longlobmin),floor(self.longlobmax),unit)
		
		for la in lats:
			row=[]
			for lo in lons:
				lower = (la,lo) ; upper = (la+unit,lo+unit)
				row.append(self.sp_per_unit(upper, lower,exclude=self.exclude))
			rich.append(row)		
		return np.array(rich), lats, lons		
	
	def species_range(self):
		'''
		Compute the range of the species using its min and max coordinates. This will affect the 
		counts if antitropical species are found. It most likely over estimate the counts. It will
		populate the richnessRange dictionary.
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		'''
		self.richnessRange={}
		for k,v in self.specciesData.iteritems():
			if k == '':
				continue 
			maxlat = max(v, key = lambda x:x[0])[0]
			minlat = min(v, key = lambda x:x[0])[0]
			maxlon = max(v, key = lambda x:x[1])[1]
			minlon = min(v, key = lambda x:x[1])[1]
			self.richnessRange[k]=[(maxlat,maxlon),(minlat,minlon)]
			
	def richnes_by_range(self,unit=5):
		'''
		Compute the richness as an occurrence count (as oppose to range count). This underestimate counts
		and have the sampling bias.
		
		:param self: intance of the class GBIFer
		:type self: :class `GBIFer`
		:param unit: The bin unit for the grid in degrees.
		:type unit: integer
		:returns: A tuple with 3 numpy arrays (richness, latitudes, longitudes)
		'''
		self.species_range()
		rich =[]
		lats = np.arange(ceil(self.latglobmin),floor(self.latglobmax),unit)
		lons = np.arange(ceil(self.longlobmin),floor(self.longlobmax),unit)
		for la in lats:
			row=[]
			
			for lo in lons:
				midgrid = (((la + (la+unit))/2), ((lo + (lo+unit))/2))
				count = 0
				for k,v in self.richnessRange.iteritems():
					if self.point_in_square(v[0], v[1], midgrid):
						count+=1
				row.append(count)
			rich.append(row)
		return np.array(rich), lats, lons
		
	
	def plot_richness(self,proj='merc',grid=5,rang=False):
		'''
		Plot richness (counts of species per grid in degrees)
		
		:param self: Intance of the class GBIFer
		:type self: :class `GBIFer`
		:param proj: Which projection to use.The option has to be pass as the value 
		for the Basemap class (check 
		http://matplotlib.org/basemap/api/basemap_api.html#module-mpl_toolkits.basemap)
		:type proj: string
		:param grid: The bin unit for the grid in degrees
		:type grid: integer
		:param rang: Wheter or not to use species range or just species points (the default)
		:type rang: boolean
		'''
		#degree = 111 km
		#Data
		if rang:
			data,lats,lons = self.richnes_by_range(unit=grid)
			title = self.prefix+' species range richness'
		else:
			data,lats,lons = self.richness_grid_by_point(grid)
			title = title = self.prefix+' species richness by occurrence'
		#data,lons = shiftgrid(180.,data,lons,start=False)
		# create the figure and axes instances.
		fig = plt.figure()
		ax = fig.add_axes([0.05,0.05,0.9,0.9])
		meanlat = np.mean(lats)
		meanlon = np.mean(lons)
		sdlat = grid#ceil(np.std(lats))
		sdlon = grid#ceil(np.std(lons))
		#setup of basemap ('lcc' = lambert conformal conic).
		# use major and minor sphere radii from WGS84 ellipsoid.
		
		urcrnrlon=self.longlobmax + (sdlon)
		if urcrnrlon < -180: urcrnrlon=-180
		if urcrnrlon >  180: urcrnrlon= 180
		
		#if np.sign(self.longlobmin) == 1:
		llcrnrlon=self.longlobmin - (sdlon)
		if llcrnrlon < -180: llcrnrlon=-180
		if llcrnrlon >  180: llcrnrlon = 180
		
		#if np.sign(self.latglobmin) == 1:
		llcrnrlat=self.latglobmin - (sdlat)
		#else:
		#	llcrnrlat=self.latglobmin + (sdlat)
		if llcrnrlat < -85: llcrnrlat=-85
		if llcrnrlat >  85: llcrnrlat= 85
		
		
		urcrnrlat=self.latglobmax + (sdlat)
			
		if urcrnrlat < -85: urcrnrlat=-85
		if urcrnrlat >  85: urcrnrlat= 85
		
		if sdlat != 0:
			lons = np.append(llcrnrlon, lons)
			data = np.insert(data,0,0,axis=1)
			lats = np.append(llcrnrlat,lats)
			data = np.insert(data,0,0,axis=0)
		if sdlon != 0:
			lons = np.append(lons,urcrnrlon)
			data = np.insert(data,-1,0,axis=1)
			lats = np.append(lats,urcrnrlat)
			data = np.insert(data,-1,0,axis=0)
		
		m = Basemap(llcrnrlon=llcrnrlon,urcrnrlon=urcrnrlon,urcrnrlat=urcrnrlat,llcrnrlat=llcrnrlat,\
		            resolution='i',area_thresh=100.,projection=proj,lat_0=meanlat,lon_0=meanlon,ax=ax) 
		# transform to nx x ny regularly spaced 5km native projection grid
		#lons,lats = m.shiftdata(lons, lats)
		nx = int((m.xmax-m.xmin)/5000.)+1; ny = int((m.ymax-m.ymin)/5000.)+1
		dat = m.transform_scalar(data,lons,lats,nx,ny)		
		# plot image over map with imshow.
		im = m.imshow(dat,cm.GMT_haxby,interpolation='catrom')
		# draw coastlines and political boundaries.
		m.fillcontinents(color='0.8')
		m.drawcoastlines()
		m.drawcountries()
		#m.drawstates()
		# draw parallels and meridians.
		# label on left and bottom of map.
		m.drawparallels(np.arange(-80,81,20),labels=[1,0,0,0])
		m.drawmeridians(np.arange(0,360,40),labels=[0,0,0,1])
		# add colorbar
		cb = m.colorbar(im,"right", size="5%", pad='2%')
		ax.set_title(title)
		plt.savefig(self.prefix+'.pdf',dpi=600)		

# END Functions bit ###############################################################################

if __name__ == "__main__":
	import os,sys,optparse
	# Command line input #############################################################
	opts = optparse.OptionParser(usage='%prog <prefix> [options]')
	opts.add_option('-f','--filename',dest='filename', action="store", default='', 
	                help='Filename of the GBIF zip to use. If none passed, all zips in'\
	                ' the current working directory will be used.')
	opts.add_option('-p','--prefix',dest='prefix',action="store",default='GBIF',
	                help='Prefix for output files.')
	opts.add_option('-g','--unitgrid',dest='unitgrid',action="store",default=5, type=int,
	                help='Define the size of the unit of a grid.')
	opts.add_option('-q','--projection',dest='proj',action="store",default='lcc',
	                help='Define the projection type. By default is Lambert'\
	                ' Conformal (lcc). The option has to be pass as the value'\
	                ' for the Basemap class (check http://matplotlib.org/basemap'\
	                '/api/basemap_api.html#module-mpl_toolkits.basemap) for options.')	
	opts.add_option('-r','--range',dest='rang',action="store_true",default=False, 
		                help='Use species range or points. Points by default.')	

	options, args = opts.parse_args()

	filename = options.filename
	prefix = options.prefix
	grid = options.unitgrid
	proj = options.proj

	I = GBIFer(prefix=prefix, filename=filename)
	I.plot_richness(grid=grid,proj=proj,rang=options.rang)