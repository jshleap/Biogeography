#!/usr/bin/python
"""
biogeographer Copyright (C) 2012 Jose Sergio Hleap

This Script will fetch results from GBIF for macoecological and biogeografical studies.
It is Based on the BioGeography is a module under development by Nick Matzke for a Google
Summer of Code 2009 project. It is run through NESCENT's Phyloinformatics Summer of Code 
2009. See the project proposal at: 
Biogeographical Phylogenetics for BioPython. 
The mentors are Stephen Smith (primary), Brad Chapman, and David Kidd. 
The source code is in the Bio/Geography directory of the Geography fork of the nmatzke branch
on GitHub, and you can see a timeline and other info about ongoing development of the module
here. The new module is being documented on the BioPython wiki as BioGeography.

The CAS (Eschmeyer, W. N. (ed.) Catalog of Fishes electronic version.
http://research.calacademy.org/research/ichthyology/catalog/fishcatmain.asp) bit code was 
written by Mike Porter.

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

E-mail: jshleap@squalus.org

"""
################################################################################################
##Importing Bit
################################################################################################

import unicodedata , string , types , time , re , htmlentitydefs, sys, urllib, httplib
import warnings, os.path, sys
from glob import glob as G
import pickle as P
import matplotlib.pyplot as plt
from numpy import array, reshape, NaN, median
from urllib import urlencode
from Bio import File
from xml.etree import ElementTree as ET
import datetime as dt

################################################################################################
'''
Definitions bit: I couln't got to work module Bio.Geography, so manually imported the
fuctions. Is also useful to make a pseuso-stand alone script without having
biopython
'''
################################################################################################
# ############################################################################################ #
# Start definitions ########################################################################## #
# ############################################################################################ #
# Start Bio.Geography ######################################################################## #
"""
Functions for accessing GBIF, downloading records, processing them into a class, and extracting information 
from the xmltree in that class.
"""
class GbifObservationRecord(Exception): pass
class GbifObservationRecord():
	"""GbifObservationRecord is a class for holding an individual observation at an individual lat/long point."""
	# Info prints out the informative description of this class, if a user wants to see it.
	#info = "Information about this class: GbifObservationRecord is a class for holding an individual observation at an individual lat/long point."
	def __init__(self):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		self.gbifkey = None
		self.catalognum = None
		self.country = None
		self.lat = None
		self.lon = None
		self.earlydate = None
		self.taxonconceptkey = None
		self.taxonname = None
		self.namecomplete = None
		self.taxon = None
		self.genus = None
		self.species = None
		self.scientific = None
		self.latedate = None
		self.area = None # This is not taken from XML, it is determined by classify_point_into_area

	def latlong_to_obj(self, line):
		"""
		Read in a string which consists of a simple table, read species/lat/long to GbifObservationRecord object
		# This can be slow, e.g. 10 seconds for even just ~1000 records.
		"""
		print line
		words = line.split("\t")
		numwords = len(words)
		self.gbifkey = int(words[0])
		self.lat = float(words[1])
		self.lon = float(words[2])
		self.taxon = str(words[3])
		self.latedate = str(words[4])
		temptaxon_words = self.taxon.split()
		if len(temptaxon_words) == 2:
			self.genus = str(temptaxon_words[0])
			self.species = str(temptaxon_words[1])
		if len(temptaxon_words) == 1:
			self.genus = self.taxon
		return

	def classify_point_into_area(self, poly, polyname):
		"""
		Fill in the self.area attribute with polygon name "poly" if it falls within the polygon "poly".  
		Otherwise, don't change self.area (which will be "None" or a previously-determined value. Uses GeogUtils library.
		"""
		x = self.lon
		y = self.lat
		if poly:
			inside = point_inside_polygon(x,y,poly) # modified by Hleap does not uses GeogUtils
		else:
			inside =  True
		if inside == True:
			self.area = polyname

		return self.area

	def parse_occurrence_element(self, element):
		"""
		Parse a TaxonOccurrence element, store in OccurrenceRecord
		"""
		# Get the GBIF record key for the occurrence
		try:
			self.gbifkey = element.attrib['gbifKey']
		except:
			pass
		# Get the catalog number
		self.catalognum = self.fill_occ_attribute(element, 'catalogNumber', 'str')
		self.country = self.fill_occ_attribute(element, 'country', 'str')
		self.lat = self.fill_occ_attribute(element, 'decimalLatitude', 'float')
		self.lon = self.fill_occ_attribute(element, 'decimalLongitude', 'float')
		self.earlydate = self.fill_occ_attribute(element, 'earliestDateCollected', 'str')
		try:
			matching_subel = self.find_1st_matching_subelement(element, 'TaxonConcept', None)
			if matching_subel is not None:
				self.taxonconceptkey = matching_subel.attrib['gbifKey']
		except:
			pass
		self.taxonname = self.fill_occ_attribute(element, 'TaxonName', 'str')			
		self.namecomplete = self.fill_occ_attribute(element, 'nameComplete', 'str')
		self.genus = self.fill_occ_attribute(element, 'genusPart', 'str')
		self.species = self.fill_occ_attribute(element, 'specificEpithet', 'str')
		self.scientific = self.fill_occ_attribute(element, 'scientific', 'str')
		self.latedate = self.fill_occ_attribute(element, 'latestDateCollected', 'str')
		return

	def fill_occ_attribute(self, element, el_tag, format='str'):
		"""
		Return the text found in matching element matching_el.text.
		"""
		return_element = None
		matching_el = self.find_1st_matching_subelement(element, el_tag, return_element)
		#print matching_el
		# Typing these variables makes it go much faster.
		if matching_el is not None:
			result = matching_el.text
			if result == '':
				return None
			elif result == 'true':
				return bool(True)
			elif result == 'false':
				return bool(False)
			elif result == None:
				#print "None check #1"
				return None
			else:
				if format == 'str':
					return str(fix_ASCII_line(result))
				else:
					textstr = format + '(' + result + ')'
					try: 
						return eval(textstr)
					except:
						textstr = "Expected type (" + format + ") not matched. String version: " + str(fix_ASCII_line(result))
						return textstr
		else:
			return None
		return None

	def find_1st_matching_subelement(self, element, el_tag, return_element):
		"""
		Burrow down into the XML tree, retrieve the first element with the matching tag.
		"""
		if element.tag.endswith(el_tag):
			return_element = element
			return return_element
		else:
			children = element.getchildren()
			if len(children) > 0:
				for child in children:
					return_element = self.find_1st_matching_subelement(child, el_tag, return_element)

					# Check if it was found
					if return_element is not None:
						#print return_element
						return return_element
			# If it wasn't found, return the empty element
			#print return_element
			return return_element

	def record_to_string(self):
		"""
		Print the attributes of a record to a string
		"""

		temp_attributes_list = [self.gbifkey, self.catalognum, self.country, self.lat, self.lon, 
				                self.earlydate, self.taxonconceptkey, self.taxonname, self.namecomplete, 
				                self.taxon, self.genus, self.species, self.scientific, self.latedate, 
				                self.area]

		str_attributes_list = []
		for item in temp_attributes_list:
			str_attributes_list.append(str(item))

		return_string = '\t'.join(str_attributes_list)

		return return_string




class GbifDarwincoreXmlString(Exception): pass
class GbifDarwincoreXmlString(str):
	"""
	GbifDarwincoreXmlString is a class for holding the xmlstring returned by a GBIF search, & processing 
	it to plain text, then an xmltree (an ElementTree). GbifDarwincoreXmlString inherits string methods 
	from str (class String).
	"""
	# Info prints out the informative description of this class, if a user wants to see it.
	info = "Information about this class: GbifDarwincoreXmlString is a class for holding the xmlstring"\
		"returned by a GBIF search, & processing it to plain text, then an xmltree (an ElementTree)."

	def __init__(self, rawstring=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		if rawstring:
			self.plainstring = rawstring
		else:
			self.plainstring = None
		#return self.plainstring



class GbifXmlTreeError(Exception): pass
class GbifXmlTree():
	"""gbifxml is a class for holding and processing xmltrees of GBIF records."""
	# Info prints out the informative description of this class, if a user wants to see it.
	info = "Information about this class: gbifxml is a class for holding and processing xmltrees"\
		"of GBIF records."
	# Used in _open function to determine 3 second waiting time for server access

	def __init__(self, xmltree=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		self.data = []
		if xmltree:
			self.xmltree = xmltree
			self.root = xmltree.getroot()
		else:
			print "GbifXmlTree.__init__(): No xmltree, so creating empty object."
			self.xmltree = None
			self.root = None

	def print_xmltree(self):
		"""
		Prints all the elements & subelements of the xmltree to screen (may require 
		fix_ASCII to input file to succeed)
		"""
		xmltree = self.xmltree

		for element in xmltree.getroot():
			print element
			self.print_subelements(element)

	def print_subelements(self, element):
		"""
		Takes an element from an XML tree and prints the subelements tag & text, and
		the within-tag items (key/value or whatnot)
		"""	
		if element.__len__() == 0:
			# Print the text beginning the tag header, and between the tags
			print fix_ASCII_line(element.tag), fix_ASCII_line(element.text)

			# Check for any key/value pairs included in the tag header,
			# and print them if they exist
			if len(element.items()) > 0:
				print "Encoded items: ", fix_ASCII_line(repr(element.items()))
			return
		elif element.__len__() > 0:
			print fix_ASCII_line(element.tag), fix_ASCII_line(element.text), "#subelements =", element.__len__()
			if len(element.items()) > 0:
				print "Encoded items: ", fix_ASCII_line(repr(element.items()))

			for subelement in element.getchildren():
				self.print_subelements(subelement)
			return

	def _element_items_to_dictionary(self, element_items):
		"""
		If the XML tree element has items encoded in the tag, e.g. key/value or
		whatever, this function puts them in a python dictionary and returns 
		them.
		"""
		if len(element_items) < 2:
			print "_element_items_to_dictionary error: < 2 items in element"
			return None
		else:
			temp_dict = {}
			for item in element_items:
				temp_dict[item[0][1]] = item[1][1]
			return temp_dict

	def extract_latlongs(self, element):
		"""
		Create a temporary pseudofile, extract lat longs to it, 
		return results as string.
		Inspired by: http://www.skymind.com/~ocrow/python_string/
		(Method 5: Write to a pseudo file)
		"""
		from cStringIO import StringIO
		file_str = StringIO()
		self._extract_latlong_datum(element, file_str)
		return file_str.getvalue()


	def _extract_latlong_datum(self, element, file_str):
		"""
		Searches an element in an XML tree for lat/long information, and the 
		complete name. Searches recursively, if there are subelements.
		file_str is a string created by StringIO in extract_latlongs() (i.e., a temp filestr)
		"""
		if element.__len__() == 0:
			return (element.tag, element.text)
		elif element.__len__() > 0:
			#print element.tag, element.text, "#subelements =", element.__len__()
			for subelement in element.getchildren():
				# Get the TaxonOccurrence
				if subelement.tag.endswith('TaxonOccurrence'):
					for item in subelement.items():
						if item[0] == 'gbifKey':
							#print item[0], item[1]
							file_str.write(item[1] + '\t')

				# Get the lat/long/name
				temptuple = self._extract_latlong_datum(subelement, file_str)
				if temptuple[0].endswith('decimalLatitude'):
					file_str.write(temptuple[1] + '\t')
				elif temptuple[0].endswith('decimalLongitude'):
					file_str.write(temptuple[1] + '\t')
				elif temptuple[0].endswith('taxonName'):
					#if temptuple[1] != '':
					file_str.write(temptuple[1] + '\t')
				elif temptuple[0].endswith('latestDateCollected'):
					#if temptuple[1] != '':
					file_str.write(str(temptuple[1]) + '\n')
			return ('tag: parent subelement', 'text: multiple subelements')


	def extract_all_matching_elements(self, start_element, el_to_match):
		"""
		Returns a list of the elements, picking elements by TaxonOccurrence; this should 
		return a list of elements equal to the number of hits.
		"""
		print ''
		print 'Running extract_all_matching_elements(start_element, ' + el_to_match + ')'
		output_list = []
		self._recursive_el_match(start_element, el_to_match, output_list)
		return output_list


	def _recursive_el_match(self, element, el_to_match, output_list):
		"""
		Search recursively through xmltree, starting with element, recording all instances of el_to_match.
		"""
		# Does the element match? If so, add to list
		if element.tag.endswith(el_to_match):
			output_list.append(element)

		# If there are NO subelements, then return...
		if element.__len__() == 0:
			return output_list
		# If there are some subelements, search them...
		elif element.__len__() > 0:
			for subelement in element.getchildren():
				output_list = self._recursive_el_match(subelement, el_to_match, output_list)
			return output_list			
		# If nothing is found ever, return the below
		return ('Error: no output_list of XML elements returned')


	def find_to_elements_w_ancs(self, el_tag, anc_el_tag):
		"""
		Burrow into XML to get an element with tag el_tag, return only those el_tags underneath a 
		particular parent element parent_el_tag
		"""
		xmltree = self.xmltree
		match_el_list = []
		for element in xmltree.getroot():
			match_el_list = self.xml_recursive_search_w_anc(element, el_tag, anc_el_tag, match_el_list)

		# How many different ancs found?
		list_ancs = []
		for tupitem in match_el_list:
			list_ancs.append(tupitem[0])

		unique_ancs = unique(list_ancs)
		print unique_ancs

		print ""
		print "Number of elements found: " + str(len(list_ancs)) + " in " + str(len(unique_ancs)) + " XML 'ancestor' (higher group) categories."

		for unique_anc in unique_ancs:
			print ""
			print "Anc: " + str(unique_anc)
			for tupitem in match_el_list:
				if tupitem[0] == unique_anc:
					print "	el: " + str(tupitem[1])

		return match_el_list


	def xml_recursive_search_w_anc(self, element, el_tag, anc_el_tag, match_el_list):
		"""
		Recursively burrows down to find whatever elements with el_tag exist inside a parent_el_tag.
		"""
		xmltree = self.xmltree

		# If the element matches the tag you are looking for...
		if element.tag == el_tag:

			# Then check if the ancestor matches
			found_anc = None
			ancestor = self._xml_burrow_up(element, anc_el_tag, found_anc)

			if ancestor == None:
				pass
			else:
				if ancestor.tag == anc_el_tag:
					match_el = element
					match_el_list.append((ancestor, match_el))

		else:
			for child in element.getchildren():
				match_el_list = self.xml_recursive_search_w_anc(child, el_tag, anc_el_tag, match_el_list)

		return match_el_list


	def create_sub_xmltree(self, element):
		"""
		Create a subset xmltree (to avoid going back to irrelevant parents)
		"""

		xmltree = ET.ElementTree(element)

		return xmltree


	def _xml_burrow_up(self, element, anc_el_tag, found_anc):
		"""
		Burrow up xml to find anc_el_tag
		"""
		xmltree = self.xmltree

		if found_anc == None:
			# Just get the direct parent of child_to_search_for
			child_to_search_for = element
			parent_element = self._return_parent_in_xmltree(child_to_search_for)

			if parent_element == None:
				return found_anc

			# Does the parent match the searched-for ancestor?		
			if parent_element.tag == anc_el_tag:
				found_anc = parent_element
			else:
				# Move a level up and search again, return if found

				found_anc = self._xml_burrow_up(parent_element, anc_el_tag, found_anc)

			return found_anc

		else:
			return found_anc



	def _xml_burrow_up_cousin(element, cousin_el_tag, found_cousin):
		"""
		Burrow up from element of interest, until a cousin is found with cousin_el_tag
		"""
		xmltree = self.xmltree

		if found_cousin == None:
			# Just get the direct parent of child_to_search_for
			child_to_search_for = element
			parent_element = self._return_parent_in_xmltree(xmltree, child_to_search_for)

			if parent_element == None:
				return found_cousin

			grandparent_element = self._return_parent_in_xmltree(parent_element)
			if grandparent_element == None:
				return found_cousin

			# Does the parent or any cousins match the searched-for ancestor?		
			for aunt in grandparent_element.getchildren():
				if aunt.tag == cousin_el_tag:
					found_cousin = aunt
					return found_cousin

			if found_cousin == None:
				# Move a level up and search again, return if found
				found_cousin = self._xml_burrow_up_cousin(parent_element, cousin_el_tag, found_cousin)

			return found_cousin

		else:
			return found_cousin



	def _return_parent_in_xmltree(self, child_to_search_for):
		"""
		Search through an xmltree to get the parent of child_to_search_for
		"""
		xmltree = self.xmltree
		returned_parent = None
		for element in xmltree.getroot():
			potential_parent = element
			returned_parent = self._return_parent_in_element(potential_parent, 
						                                     child_to_search_for, returned_parent)
			return returned_parent



	def _return_parent_in_element(self, potential_parent, child_to_search_for, returned_parent):
		"""
		Search through an XML element to return parent of child_to_search_for
		"""
		if returned_parent == None:
			children = potential_parent.getchildren()
			if len(children) > 0:
				for child in potential_parent.getchildren():
					if child == child_to_search_for:
						returned_parent = potential_parent

					# If not found at this level, go down a level
					else:
						returned_parent = self._return_parent_in_element(child, child_to_search_for, returned_parent)

			return returned_parent
		else:
			return returned_parent


	def find_1st_matching_element(self, element, el_tag, return_element):
		"""
		Burrow down into the XML tree, retrieve the first element with the matching tag
		"""
		if return_element == None:
			if element.tag == el_tag:
				return_element = element
			else:
				children = element.getchildren()
				if len(children) > 0:
					for child in children:
						return_element = self.find_1st_matching_element(child, el_tag, return_element)
						if return_element is not None:
							return return_element
			return return_element	
		else:
			return return_element


	def extract_numhits(self, element):
		"""
		# Search an element of a parsed XML string and find the 
		# number of hits, if it exists.  Recursively searches, 
		# if there are subelements.
		# 
		"""
		# print "Running extract_numhits(element)..."
		if element.__len__() == 0:
			if len(element.items()) > 0:
				for item in element.items():
					for index, tupleitem in enumerate(item):
						if tupleitem == 'totalMatched':
							#print item
							#print int(item[index+1])
							return int(item[index+1])
						else:
							temp_return_item = None
		elif element.__len__() > 0:
			#print element.tag, element.text, "#subelements =", element.__len__()
			for subelement in element.getchildren():
				temp_return_item = self.extract_numhits(subelement)
			if temp_return_item != None:
				return_item = temp_return_item
				return return_item
			else:
				return return_item



class GbifSearchResults(Exception): pass
class GbifSearchResults():	
	"""
	GbifSearchResults is a class for holding a series of GbifObservationRecord records, and processing
	them e.g. into classified areas.
	"""
	# Info prints out the informative description of this class, if a user wants to see it.
	info = "GbifSearchResults is a class for holding a series of GbifObservationRecord records,"\
		" and processing them	e.g. into classified areas."

	def __init__(self, gbif_recs_xmltree=None):
		"""
		This is an instantiation class for setting up new objects of this class.
		"""
		# These are used by the _open function, e.g. to set wait times for 
		# accessing GBIF
		self.time_last_gbif_access = 0
		self.email = None
		# xmltree of GBIF records (GbifXmlTree class)	
		if gbif_recs_xmltree:
			self.gbif_recs_xmltree = gbif_recs_xmltree
		else:
			self.gbif_recs_xmltree = GbifXmlTree()
		# list of xmltrees returned by get_all_records_by_increment
		self.gbif_xmltree_list = []
		# xmltree of GBIF count query
		self.gbif_count_xmltree = None
		# String (tempfile-ish) list of records
		# An GbifDarwincoreXmlString object, with methods for cleaning the string
		# This could a list if multiple downloads are required
		self.obs_recs_xmlstring = None
		self.count_recs_xmlstring = None
		# Summary statistics for search
		self.gbif_hits_count = None
		self.gbif_count_params = None
		self.count_params_str = None
		# List of GbifObservationRecord objects
		self.obs_recs_list = []

	def print_records(self):
		"""
		Print all records in tab-delimited format to screen.
		"""

		for index, record in enumerate(self.obs_recs_list):
			# Get string for record
			print_string = record.record_to_string()

			# Fix ASCII
			print_string2 = fix_ASCII_line(print_string)

			# Print it
			print str(index+1) + '\t' + print_string

		return


	def print_records_to_file(self, fn):
		"""
		Print the attributes of a record to a file with filename fn
		"""
		# Open the file
		fh = open(fn, 'w')
		for index, record in enumerate(self.obs_recs_list):
			print_string = record.record_to_string()
			# Write to the file
			fh.write(str(index+1) + '\t' + print_string + '\n')
		# Close the file
		fh.close()
		return fn


	def classify_records_into_area(self, poly, polyname):
		"""
		Take all of the records in the GbifSearchResults object, fill in their area attribute if they fall within the polygon poly.
		"""
		if self.obs_recs_list == []:
			print 'Error: No records stored in self.obs_recs_list.'
			return
		matching_count = 0
		for record in self.obs_recs_list:
			record.classify_point_into_area(poly, polyname)
			if record.area == polyname:
				matching_count = matching_count + 1
		print str(matching_count),"of",str(len(self.obs_recs_list)),
		'fell inside area"' + polyname + '".'

		return


	def latlongs_to_obj(self):
		"""
		Takes the string from extract_latlongs, puts each line into a
		GbifObservationRecord object.
		Return a list of the objects
		"""
		if self.obs_recs_xmlstring == None:
			# If there is no string of records yet, extract the string,
			# using the gbif_recs_xmltree root (getroot from xmltree) as the input
			# element
			print "Filling in self.obs_recs_xmlstring"
			self.obs_recs_xmlstring = self.gbif_recs_xmltree.extract_latlongs(self.gbif_recs_xmltree.root)
			#print type(self.obs_recs_xmlstring)
			# Print the header
			for index, line in enumerate(self.obs_recs_xmlstring.splitlines()):
				print "Line" + str(index) + ": " + line
				if index > 4:
					print '...'
					break

		for line in self.obs_recs_xmlstring.splitlines():
			temprec = GbifObservationRecord()
			temprec.latlong_to_obj(line)
			self.obs_recs_list.append(temprec)
			#del temprec
		return


	# Functions devoted to accessing/downloading GBIF records
	def access_gbif(self, url, params):
		"""
		# Helper function to access various GBIF services
		# 
		# choose the URL ("url") from here:
		# http://data.gbif.org/ws/rest/occurrence
		#
		# params are a dictionary of key/value pairs
		#
		# "self._open" is from Bio.Entrez.self._open, online here: 
		# http://www.biopython.org/DIST/docs/api/Bio.Entrez-pysrc.html#self._open
		#
		# Get the handle of results
		# (looks like e.g.: <addinfourl at 75575128 whose fp = <socket._fileobject object at 0x48117f0>> )

		# (open with results_handle.read() )
		"""
		print 'Accessing GBIF with access_gbif...'
		results_handle = self._open(url, params)
		return results_handle


	def _get_hits(self, params):
		"""
		Get the actual hits that are be returned by a given search
		(this allows parsing & gradual downloading of searches larger 
		than e.g. 1000 records)

		It will return the LAST non-none instance (in a standard search result there
		should be only one, anyway).
		"""
		print '  Running _get_hits(params)...'

		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/list'

		cmd = url + self._paramsdict_to_string(params)
		results_handle = self.access_gbif(url, params)

		self.obs_recs_xmlstring = GbifDarwincoreXmlString(results_handle.read())

		print 'XML search results stored via: "self.obs_recs_xmlstring = GbifDarwincoreXmlString(results_handle.read())"'

		return self.obs_recs_xmlstring


	def get_xml_hits(self, params):
		"""
		Returns hits like _get_hits, but returns a parsed XML tree.
		"""
		print ''
		print 'Running get_xml_hits(params)...'
		# Returns GbifDarwincoreXmlString object
		# (stored in self.obs_recs_xmlstring)
		self._get_hits(params)
		# Fix the xmlstring:
		# returns a plain string
		plain_xmlstring_fixed = self.obs_recs_xmlstring.plainstring
		#plain_xmlstring_fixed = fix_ASCII_lines(self.obs_recs_xmlstring)
		# (temp_xmlstring.splitlines() works because GbifDarwincoreXmlString inherits from string class)
		# Store the plain string in an GbifDarwincoreXmlString object, store as method of
		# GbifXmlTree object
		self.obs_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)
		# make sure the self.obs_recs_xmlstring (an GbifDarwincoreXmlString object) parses into an ElementTree
		try:
			self.gbif_recs_xmltree = GbifXmlTree(ET.ElementTree(ET.fromstring(self.obs_recs_xmlstring)))
			print_it = 0
			if print_it == 1:
				print ''
				print "Printing xmlstring_fixed..."
				fixed_string = fix_ASCII_lines(self.obs_recs_xmlstring)
				print fixed_string
				print ''
		except Exception, inst:
			print "Unexpected error opening %s: %s" % ('"self.obs_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)"', inst)


		return 

	def get_record(self, key):
		"""
		Given the key, get a single record, return xmltree for it.
		"""
		print ''
		print 'Running get_record(params)...'
		# URL for the record utility
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/get'
		params = {'format': 'darwin', 'key' : key}
		cmd = url + self._paramsdict_to_string(params)
		results_handle = self.access_gbif(url, params)
		#print results_handle.read()
		xmlstring = GbifDarwincoreXmlString(results_handle.read())
		# returns plain string, with linebreaks for parsing with ET.fromstring
		xmlstring2 = fix_ASCII_lines(xmlstring.plainstring)
		xmltree = ET.ElementTree(ET.fromstring(xmlstring))
		temp = GbifXmlTree(xmltree)

		return xmltree


	def get_numhits(self, params):
		"""
		Get the number of hits that will be returned by a given search
		(this allows parsing & gradual downloading of searches larger 
		than e.g. 1000 records)
		It will return the LAST non-none instance (in a standard search result there
		should be only one, anyway).
		"""
		self.gbif_count_params = params
		print ''
		print 'Running get_numhits(params)...'
		# URL for the count utility
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/count'
		self.count_params_str = self._paramsdict_to_string(params)
		cmd = url + self.count_params_str
		results_handle = self.access_gbif(url, params)
		self.count_recs_xmlstring = GbifDarwincoreXmlString(results_handle.read())
		# Fix xmlstring
		plain_xmlstring_fixed = self.count_recs_xmlstring
		# Store the plain string in an GbifDarwincoreXmlString object, store as method of
		# GbifXmlTree object
		self.count_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)
		# make sure the self.obs_recs_xmlstring (an GbifDarwincoreXmlString object) parses into an ElementTree
		try:
			self.gbif_count_xmltree = GbifXmlTree(ET.ElementTree(ET.fromstring(self.count_recs_xmlstring)))
		except Exception, inst:
			print "Unexpected error opening %s: %s" % ('"self.count_recs_xmlstring = GbifDarwincoreXmlString(plain_xmlstring_fixed)"', inst)
		# Get the element from the xmltree attribute of the GbifXmlTree object stored in self
		for element in self.gbif_count_xmltree.xmltree.getroot():
			temp_numhits = self.gbif_count_xmltree.extract_numhits(element)
			if temp_numhits != None:
				self.gbif_hits_count = int(temp_numhits)

		print '# hits on "' + self.count_params_str + '" = ' + str(self.gbif_hits_count)

		return self.gbif_hits_count



	def xmlstring_to_xmltree(self, xmlstring):
		"""
		Take the text string returned by GBIF and parse to an XML tree using ElementTree.  
		Requires the intermediate step of saving to a temporary file (required to make
		ElementTree.parse work, apparently)
		"""
		# instructions for ElementTree:
		# http://docs.python.org/library/xml.etree.elementtree.html
		# make sure the file is findable
		try:
			xmltree = ET.fromstring(xmlstring)
		except Exception, inst:
			print "Unexpected error opening %s: %s" % ('file_str', inst)

		# make sure the file is parsable
		try:
			xmltree.getroot()
		except Exception, inst:
			print "Unexpected error running getroot() on text in file %s: %s" % (tempfn, inst)

		return xmltree


	def get_all_records_by_increment(self, params, inc):
		"""
		Download all of the records in stages, store in list of elements.
		Increments of e.g. 100 to not overload server
		"""
		print ''
		print "Running get_all_records_by_increment(params, inc)"
		# Set up the list of chunks to get
		numhits = self.get_numhits(params)
		print "#hits = ", numhits
		list_of_chunks = range(0, numhits-1, inc)
		# Set up list of xmlstring results
		self.gbif_xmltree_list = []
		# Set up url to parse
		# instructions: http://data.gbif.org/ws/rest/occurrence
		url = 'http://data.gbif.org/ws/rest/occurrence/list'
		# download them by increment	
		for index, startindex in enumerate(list_of_chunks):
			if startindex + inc  > numhits:
				print "Downloading records #" + str(startindex) + "-" + str(numhits-1) , 
				"of " + str(numhits) +" (-1)."
				params['startindex'] = str(startindex)
				params['maxresults'] = str(numhits - startindex)
			else:
				print "Downloading records #" + str(startindex) + "-" + str(startindex+inc-1),
				"of " + str(numhits) +" (-1)."
				params['startindex'] = str(startindex)
				params['maxresults'] = str(inc)
			print params
			results_handle = self.access_gbif(url, params)
			# returns GbifDarwincoreXmlString
			temp_xmlstring = GbifDarwincoreXmlString(results_handle.read())
			# returns plain string
			#xmlstring2 = temp_xmlstring.fix_ASCII_lines('\n')
			xmlstring2 = temp_xmlstring
			#Bug-checking; make True if you want to see the xmlstring as a text file
			xmlstring_tofile = False
			if xmlstring_tofile == True:
				error_fh = open('error2.txt', 'w')
				error_fh.write(xmlstring2)
				error_fh.close()
			# Close results_handle
			results_handle.close()
			#try:				
			# Add these latlongs to filestr
			gbif_xmltree = GbifXmlTree(ET.ElementTree(ET.fromstring(xmlstring2)))
			# Append the xmltree
			self.gbif_xmltree_list.append(gbif_xmltree)
			# Search through it to get the occurrences to add
			self.extract_occurrences_from_gbif_xmltree(gbif_xmltree)
		return self.gbif_xmltree_list



	def extract_occurrences_from_gbif_xmltree(self, gbif_xmltree):
		"""
		Extract all of the 'TaxonOccurrence' elements to a list, store them in a GbifObservationRecord.
		"""
		start_element = gbif_xmltree.xmltree.getroot()
		# Element to pull out:
		el_to_match = 'TaxonOccurrence'
		# Get all occurrences:
		occurrences_list = gbif_xmltree.extract_all_matching_elements(start_element, el_to_match)

		#print occurrences_list

		# For each one, extract info
		for element in occurrences_list:
			# Make a temporary observation record
			temp_observation = GbifObservationRecord()

			# Populate it from the element
			temp_observation.parse_occurrence_element(element)

			#print temp_observation.record_to_string()

			# Add the observation to the list
			self.obs_recs_list.append(temp_observation)

		return


	def _paramsdict_to_string(self, params):
		"""
		# Converts the python dictionary of search parameters into a text 
		# string for submission to GBIF
		"""
		temp_outstring_list = []
		for key in params.keys():
			temp_outstring_list.append(str(key) + '=' + str(params[key]))
		outstring = '&'.join(temp_outstring_list)
		return outstring


	def _open(self, cgi, params={}):
		"""
		Function for accessing online databases.

		Modified from: 
		http://www.biopython.org/DIST/docs/api/Bio.Entrez-module.html

		Helper function to build the URL and open a handle to it (PRIVATE).

		Open a handle to GBIF.  cgi is the URL for the cgi script to access.
		params is a dictionary with the options to pass to it.  Does some
		simple error checking, and will raise an IOError if it encounters one.

		This function also enforces the "three second rule" to avoid abusing
		the GBIF servers (modified after NCBI requirement).
		"""
		# NCBI requirement: At least three seconds between queries
		delay = 3.0
		current = time.time()

		wait = self.time_last_gbif_access + delay - current
		if wait > 0:
			time.sleep(wait)
			self.time_last_gbif_access = current + wait
		else:
			self.time_last_gbif_access = current
		# Remove None values from the parameters
		for key, value in params.items():
			if value is None:
				del params[key]

		# Eliminating this bit; irrelevant in GBIF
		# Tell Entrez that we are using Biopython
		"""
		if not "tool" in params:
			params["tool"] = "biopython"
		"""
		# Tell Entrez who we are
		if not "email" in params:
			if self.email != None:
				params["email"] = email
		# Open a handle to Entrez.
		options = urllib.urlencode(params, doseq=True)
		cgi += "?" + options
		handle = urllib.urlopen(cgi)

		# Wrap the handle inside an UndoHandle.
		uhandle = File.UndoHandle(handle)

		# Check for errors in the first 5 lines.
		# This is kind of ugly.
		lines = []
		for i in range(5):
			lines.append(uhandle.readline())
		for i in range(4, -1, -1):
			uhandle.saveline(lines[i])
		data = ''.join(lines)

		if "500 Proxy Error" in data:
			# Sometimes Entrez returns a Proxy Error instead of results
			raise IOError("500 Proxy Error (NCBI busy?)")
		elif "502 Proxy Error" in data:
			raise IOError("502 Proxy Error (NCBI busy?)")
		elif "WWW Error 500 Diagnostic" in data:
			raise IOError("WWW Error 500 Diagnostic (NCBI busy?)")
		elif data.startswith("Error:") :
			#e.g. 'Error: Your session has expired. Please repeat your search.\n'
			raise IOError(data.strip())
		elif data.startswith("The resource is temporarily unavailable") :
			#This can occur with an invalid query_key
			#Perhaps this should be a ValueError?
			raise IOError("The resource is temporarily unavailable")
		elif data.startswith("download dataset is empty") :
			#This can occur when omit the identifier, or the WebEnv and query_key
			#Perhaps this should be a ValueError?
			raise IOError("download dataset is empty")
		elif data[:5] == "ERROR":
			# XXX Possible bug here, because I don't know whether this really
			# occurs on the first line.  I need to check this!
			raise IOError("ERROR, possibly because id not available?")
		# Should I check for 404?  timeout?  etc?
		return uhandle


# END Bio.Geography ############################################################################################

"""
Some simple functions using free libraries for reading shapefiles, performing point-in-polygon 
operations, etc.
"""

def readshpfile(fn):

	# Get the filename
	# fn = shapefilename
	# Try to open the file, print error if error	
	try:
		f = open(fn)
	# some code here
	except IOError, reason:
		print "couldn't open file, because of", reason
		return
	# load the shapefile, populating a list of dictionaries
	shpRecords = shpUtils.loadShapefile(fn)
	return shpRecords


def summarize_shapefile(fn, output_option, outfn):
	shpRecords = readshpfile(fn)

	if output_option == 'toscreen':
		print "Feature#	Name	X	Y"

		for index, record in enumerate(shpRecords):

			# Make an output string

			# consecutive ount
			pt1 = str(index+1)

			# Feature name
			pt2 = str(record['dbf_data']['NAME'])

			# X & Y coords
			pt3 = str(record['shp_data']['x'])
			pt4 = str(record['shp_data']['y'])

			printstr = '	'.join([pt1, pt2, pt3, pt4])
			print printstr

		return "printed string to screen"


	elif output_option == 'tofile':
		outfile = open(outfn, 'w')

		header_str = '	'.join(['Feature#', 'Name', 'X', 'Y'])
		outfile.write(header_str)


		for index, record in enumerate(shpRecords):

			# Make an output string

			# consecutive ount
			pt1 = str(index+1)

			# Feature name
			pt2 = str(record['dbf_data']['NAME'])

			# X & Y coords
			pt3 = str(record['shp_data']['x'])
			pt4 = str(record['shp_data']['y'])

			printstr = '	'.join([pt1, pt2, pt3, pt4])
			outfile.write(printstr)

		outfile.close()
		return outfn


	elif output_option == 'tolist':

		# Create blank list
		outlist = []

		# Determine header row
		header_row = ['Feature#', 'Name', 'X', 'Y']
		outlist.append(header_row)

		for index, record in enumerate(shpRecords):

			# Make an output string

			# consecutive ount
			pt1 = str(index+1)

			# Feature name
			pt2 = str(record['dbf_data']['NAME'])

			# X & Y coords
			pt3 = str(record['shp_data']['x'])
			pt4 = str(record['shp_data']['y'])

			line = [pt1, pt2, pt3, pt4]
			outlist.append(line)

		return outlist

	return


def point_inside_polygon(x,y,poly):
	"""
	# Written by Paul Bourke in C (http://local.wasp.uwa.edu.au/~pbourke/geometry/insidepoly/), 
	# translated to python by Patrick Jordan (http://www.ariel.com.au/a/python-point-int-poly.html)
	# NOTE: does not presently deal with the case of polygons that
	# cross the international dateline!
	# determine if a point is inside a given polygon or not
	# Polygon is a list of (x,y) pairs.
	"""
	n = len(poly)
	inside = False

	p1x,p1y = poly[0]
	for i in range(n+1):
		p2x,p2y = poly[i % n]
		if y > min(p1y,p2y):
			if y <= max(p1y,p2y):
				if x <= max(p1x,p2x):
					if p1y != p2y:
						xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
					if p1x == p2x or x <= xinters:
						inside = not inside
		p1x,p1y = p2x,p2y

	return inside


"""
Generic functions for simple operations.
"""
# Incorporates:
# Public domain module to convert diverse characters to ASCII (dammit!)
# http://newsbruiser.tigris.org/source/browse/~checkout~/newsbruiser/nb/lib/AsciiDammit.py
# from AsciiDammit import asciiDammit
"""ASCII, Dammit

Stupid library to turn MS chars (like smart quotes) and ISO-Latin
chars into ASCII, dammit. Will do plain text approximations, or more
accurate HTML representations. Can also be jiggered to just fix the
smart quotes and leave the rest of ISO-Latin alone.

Sources:
 http://www.cs.tut.fi/~jkorpela/latin1/all.html
 http://www.webreference.com/html/reference/character/isolat1.html

1.0 Initial Release (2004-11-28)

The author hereby irrevocably places this work in the public domain.
To the extent that this statement does not divest the copyright,
the copyright holder hereby grants irrevocably to every recipient
all rights in this work otherwise reserved under copyright.
"""

__author__ = "Leonard Richardson (leonardr@segfault.org)"
__version__ = "$Revision: 1.3 $"
__date__ = "$Date: 2009/04/28 10:45:03 $"
__license__ = "Public domain"

CHARS = { '\x80' : ('EUR', 'euro'), '\x81' : ' ', '\x82' : (',', 'sbquo'), '\x83' : ('f', 'fnof'),
          '\x84' : (',,', 'bdquo'), '\x85' : ('...', 'hellip'), '\x86' : ('+', 'dagger'),
          '\x87' : ('++', 'Dagger'),'\x88' : ('^', 'caret'),'\x89' : '%','\x8A' : ('S', 'Scaron'),
          '\x8B' : ('<', 'lt;'), '\x8C' : ('OE', 'OElig'), '\x8D' : '?', '\x8E' : 'Z', '\x8F' : '?',
          '\x90' : '?', '\x91' : ("'", 'lsquo'), '\x92' : ("'", 'rsquo'), '\x93' : ('"', 'ldquo'),
          '\x94' : ('"', 'rdquo'), '\x95' : ('*', 'bull'), '\x96' : ('-', 'ndash'), '\x97' : ('--', 'mdash'),
          '\x98' : ('~', 'tilde'), '\x99' : ('(TM)', 'trade'), '\x9a' : ('s', 'scaron'),
          '\x9b' : ('>', 'gt'), '\x9c' : ('oe', 'oelig'), '\x9d' : '?', '\x9e' : 'z', '\x9f' : ('Y', 'Yuml'),
          '\xa0' : (' ', 'nbsp'), '\xa1' : ('!', 'iexcl'), '\xa2' : ('c', 'cent'), '\xa3' : ('GBP', 'pound'),
          '\xa4' : ('$', 'curren'), #This approximation is especially lame.
          '\xa5' : ('YEN', 'yen'), '\xa6' : ('|', 'brvbar'), '\xa7' : ('S', 'sect'), '\xa8' : ('..', 'uml'),
          '\xa9' : ('', 'copy'),'\xaa' : ('(th)', 'ordf'),'\xab' : ('<<', 'laquo'), '\xac' : ('!', 'not'),
          '\xad' : (' ', 'shy'), '\xae' : ('(R)', 'reg'), '\xaf' : ('-', 'macr'), '\xb0' : ('o', 'deg'),
          '\xb1' : ('+-', 'plusmm'),'\xb2' : ('2', 'sup2'),'\xb3' : ('3', 'sup3'), '\xb4' : ("'", 'acute'),
          '\xb5' : ('u', 'micro'), '\xb6' : ('P', 'para'), '\xb7' : ('*', 'middot'), '\xb8' : (',', 'cedil'),
          '\xb9' : ('1', 'sup1'), '\xba' : ('(th)', 'ordm'), '\xbb' : ('>>', 'raquo'), '\xbc' : ('1/4', 'frac14'),
          '\xbd' : ('1/2', 'frac12'), '\xbe' : ('3/4', 'frac34'), '\xbf' : ('?', 'iquest'),          
          '\xc0' : ('A', "Agrave"), '\xc1' : ('A', "Aacute"), '\xc2' : ('A', "Acirc"),
          '\xc3' : ('A', "Atilde"), '\xc4' : ('A', "Auml"), '\xc5' : ('A', "Aring"), '\xc6' : ('AE', "Aelig"), 
          '\xc7' : ('C', "Ccedil"), '\xc8' : ('E', "Egrave"), '\xc9' : ('E', "Eacute"), '\xca' : ('E', "Ecirc"), 
          '\xcb' : ('E', "Euml"), '\xcc' : ('I', "Igrave"), '\xcd' : ('I', "Iacute"), '\xce' : ('I', "Icirc"), 
          '\xcf' : ('I', "Iuml"), '\xd0' : ('D', "Eth"), '\xd1' : ('N', "Ntilde"), '\xd2' : ('O', "Ograve"), 
          '\xd3' : ('O', "Oacute"), '\xd4' : ('O', "Ocirc"), '\xd5' : ('O', "Otilde"), '\xd6' : ('O', "Ouml"), 
          '\xd7' : ('*', "times"), '\xd8' : ('O', "Oslash"), '\xd9' : ('U', "Ugrave"), '\xda' : ('U', "Uacute"), 
          '\xdb' : ('U', "Ucirc"), '\xdc' : ('U', "Uuml"), '\xdd' : ('Y', "Yacute"), '\xde' : ('b', "Thorn"), 
          '\xdf' : ('B', "szlig"), '\xe0' : ('a', "agrave"), '\xe1' : ('a', "aacute"), '\xe2' : ('a', "acirc"), 
          '\xe3' : ('a', "atilde"), '\xe4' : ('a', "auml"), '\xe5' : ('a', "aring"), '\xe6' : ('ae', "aelig"), 
          '\xe7' : ('c', "ccedil"), '\xe8' : ('e', "egrave"), '\xe9' : ('e', "eacute"), '\xea' : ('e', "ecirc"), 
          '\xeb' : ('e', "euml"), '\xec' : ('i', "igrave"), '\xed' : ('i', "iacute"), '\xee' : ('i', "icirc"), 
          '\xef' : ('i', "iuml"), '\xf0' : ('o', "eth"), '\xf1' : ('n', "ntilde"), '\xf2' : ('o', "ograve"), 
          '\xf3' : ('o', "oacute"), '\xf4' : ('o', "ocirc"), '\xf5' : ('o', "otilde"), '\xf6' : ('o', "ouml"), 
          '\xf7' : ('/', "divide"), '\xf8' : ('o', "oslash"), '\xf9' : ('u', "ugrave"), '\xfa' : ('u', "uacute"), 
          '\xfb' : ('u', "ucirc"), '\xfc' : ('u', "uuml"), '\xfd' : ('y', "yacute"), '\xfe' : ('b', "thorn"), 
          '\xff' : ('y', "yuml"),
          }

def _makeRE(limit):
	"""Returns a regular expression object that will match special characters
	up to the given limit."""
	return re.compile("([\x80-\\x%s])" % limit, re.M)
ALL = _makeRE('ff')
ONLY_WINDOWS = _makeRE('9f')

def _replHTML(match):
	"Replace the matched character with its HTML equivalent."
	return _repl(match, 1)

def _repl(match, html=0):
	"Replace the matched character with its HTML or ASCII equivalent."
	g = match.group(0)
	a = CHARS.get(g,g)
	if type(a) == types.TupleType:
		a = a[html]
		if html:
			a = '&' + a + ';'
	return a

def _dammit(t, html=0, fixWindowsOnly=0):
	"Turns ISO-Latin-1 into an ASCII representation, dammit."

	r = ALL
	if fixWindowsOnly:
		r = ONLY_WINDOWS
	m = _repl
	if html:
		m = _replHTML

	return re.sub(r, m, t)

def asciiDammit(t, fixWindowsOnly=0):
	"Turns ISO-Latin-1 into a plain ASCII approximation, dammit."
	return _dammit(t, 0, fixWindowsOnly)

def htmlDammit(t, fixWindowsOnly=0):
	"Turns ISO-Latin-1 into plain ASCII with HTML codes, dammit."
	return _dammit(t, 1, fixWindowsOnly=fixWindowsOnly)

def demoronise(t):
	"""Helper method named in honor of the original smart quotes
	remover, The Demoroniser:

	http://www.fourmilab.ch/webtools/demoroniser/"""
	return asciiDammit(t, 1)

# ==========================
# non-AsciiDammit functions (by Nick Matzke) follow...

def make_NaN_array(xdim, ydim):
	"""
	Make an empty floating-point array with the specified dimensions
	"""
	temparray = array( [NaN] * (xdim * ydim), dtype=float)
	temparray = reshape(temparray, (xdim, ydim))
	return temparray


def make_None_list_array(xdim, ydim):
	"""
	Make a list of lists ("array") with the specified dimensions
	"""

	# Create empty array
	temp_row_array = [None] * xdim
	temp_cols_array = []
	for i in range(0, ydim):
		temp_cols_array.append(temp_row_array)

	return temp_cols_array

def set_diags_to_none(list_array):
	"""
	Given a square list of lists (rows), set the diagonal to None.
	"""
	ydim = len(list_array)
	xdim = len(list_array[0])

	if xdim != ydim:
		print ''
		print 'ERROR: This list_array is not a square!'

	for index1 in range(0, ydim-1):
		for index2 in range(0, xdim-1):
			list_array[index1, index2] = None

	return list_array


def list1_items_in_list2(list1, list2):
	"""
	Returns the list of list1 items which are found in list2
	http://vermeulen.ca/python-techniques.html
	"""
	intersection = filter(lambda x:x in list1, list2)
	return intersection

def list1_items_not_in_list2(list1, list2):
	"""
	Returns the list of list1 items which are NOT found in list2
	http://vermeulen.ca/python-techniques.html
	"""
	difference=filter(lambda x:x not in list2,list1)
	return difference


def remove_nonxml_brackets(line):
	"""
	When a GBIF DarwinCore-formatted XML record is converted to ASCII, links in some of the metadata (e.g. "<rap.conservation.org>") can get interpreted as XML tags, which are then mismatched, causing a crash.  This converts those to e.g. "[rap.conservation.org]".
	"""



	match_strings = re.findall(r'<.*?>', line)

	list_of_valid_xml = ['<gbif', '</gbif', '<?xml', '<to:', '</to:', '<tc:', '</tc:', '<tn:', '</tn:', '<tcom:', '</tcom:', '<tax', '</tax', '<data', '</data', '<name', '</name', '<occurrence', '</occurrence' ]


	# If no match, continue to next loop
	if match_strings != []:		
		# If there is a match:
		for tempmatch in match_strings:
			# Skip legit tags
			valid_tag_found = False
			for valid_tag in list_of_valid_xml:
				if tempmatch.startswith(valid_tag):
					# if it matches, go to next tempmatch
					valid_tag_found = True
			if valid_tag_found == False:
				# if not, replace
				replacement = tempmatch.replace('<', '[')
				replacement = replacement.replace('>', ']')
				print 'Replacing "' + tempmatch + '" with "' + replacement + '"'
				line = line.replace(tempmatch, replacement)
				continue


	return line



def unescape(text):
	"""
	##
	# Removes HTML or XML character references and entities from a text string.
	#
	# @param text The HTML (or XML) source text.
	# @return The plain text, as a Unicode string, if necessary.
	# source: http://effbot.org/zone/re-sub.htm#unescape-html
	"""
	def fixup(m):
		text = m.group(0)
		if text[:2] == "&#":
			# character reference
			try:
				if text[:3] == "&#x":
					return unichr(int(text[3:-1], 16))
				else:
					return unichr(int(text[2:-1]))
			except ValueError:
				pass
		else:
			# named entity
			try:
				text = unichr(htmlentitydefs.name2codepoint[text[1:-1]])
			except KeyError:
				pass
		return text # leave as is
	return re.sub("&#?\w+;", fixup, text)


def fix_ampersand(line):
	"""
	Replaces "&" with "&amp;" in a string; this is otherwise 
	not caught by the unescape and unicodedata.normalize functions.
	"""
	return line.replace('&', '&amp;')



def fix_ASCII_file(fn_unfixed):
	"""
	# Search-replace to fix annoying 
	# non-ASCII characters in search results
	# 
	# inspiration:
	# http://www.amk.ca/python/howto/unicode
	# http://www.peterbe.com/plog/unicode-to-ascii
	"""
	if fn_unfixed.count('.') < 1:
		return None

	# Make a new filename, insert '_fixed' before the rightmost period
	fn_fixed_tuple = fn_unfixed.rpartition('.')
	fn_fixed = ''.join( [fn_fixed_tuple[0], '_fixed.', fn_fixed_tuple[2]] )

	# Open unfixed file and fix all non-ASCII characters
	fh = open(fn_unfixed, 'r')
	fh_fixed = open(fn_fixed, 'w')

	lines = fh.readlines()
	newlines = fix_ASCII_lines(lines)

	fh_fixed.write(newlines)
	fh_fixed.close()
	fh.close()
	return fn_fixed


def fix_ASCII_lines(lines):
	"""
	Convert each line in an input string into pure ASCII
	(This avoids crashes when printing to screen, etc.)
	"""
	# If type is string:
	if lines.__class__ == 'string'.__class__:
		lines2 = lines.split('\n')
	# If type is list:
	elif lines.__class__ == [1,2].__class__:
		lines2 = lines
	else:
		print ''
		print 'Error: fix_ASCII_lines takes only strings with endlines, or lists.'

	newstr_list = []
	for line in lines2:
		ascii_content5 = fix_ASCII_line(line)

		newstr_list.append(ascii_content5 + '\n')
	return ''.join(newstr_list)


def fix_ASCII_line(line):
	"""
	Run several functions that fix common problems found when printing e.g. GBIF's XML results to screen: non-ASCII characters, ampersands, and <links in brackets> which can get misinterpreted as unmatched XML tags.
	"""

	# Catch error of an empty string.
	if line == None:
		line = "None"
	if line == '':
		line = ''

	# library from here: http://effbot.org/zone/re-sub.htm#unescape-html
	ascii_content1 = unescape(line)

	# Public domain module to convert diverse characters to ASCII
	# http://newsbruiser.tigris.org/source/browse/~checkout~/newsbruiser/nb/lib/AsciiDammit.py
	ascii_content2 = asciiDammit(ascii_content1)

	# inspiration: http://www.amk.ca/python/howto/unicode
	ascii_content3 = unicodedata.normalize('NFKC', unicode(ascii_content2)).encode('ascii','ignore')

	# Fix the ampersand
	ascii_content4 = fix_ampersand(ascii_content3)
	ascii_content5 = remove_nonxml_brackets(ascii_content4)

	return ascii_content5

def element_text_to_string(txt):
	if txt == None:
		txt = ""
	return str(txt).strip()


def element_items_to_string(items):
	"""
	Input a list of items, get string back
	"""
	s = ""
	for item in items:
		s = s + " " + str(item)
		s = s.strip()

	return s


def get_str_subset(start, end, seq):

	index1 = start-1
	index2 = end

	newstring = str()

	for i in range(index1, index2):
		newstring = newstring + seq[i]

	#print 'len(newstring)=', str(len(newstring))
	#print 'start - stop', str(end - start + 1) 

	return newstring


def find_1st_match(list1, list2):
	"""
	Find the first match in two ordered lists.
	"""
	anc_list1 = list1
	anc_list2 = list2

	for anc1 in anc_list1:
		for anc2 in anc_list2:
			if anc1==anc2:
				return(anc1)

	return None



def kml_to_poly(filename):
	'''
	Takes a kml file generated by Google earth or similar and extract the polygon
	'''
	poly=[]
	fn = open(filename)
	f  = fn.read()
	x  = f.split('coordinates>')
	co = x[1].replace('\t','').replace('</','').strip()[:-2]
	co = co.split(',0 ')
	for e in co:
		e=e.split(',')
		t=(float(e[0]),float(e[1]))
		poly.append(t)

	return poly



class individual:
	'''
	Get individual records into biogeographer data structure
	'''
	def __init__(self,gbifrecord):
		self.record   = gbifrecord
		self.gbif_key = None
		self.country  = None
		self.species  = None
		self.lat      = None
		self.lon      = None
		self.area     = None
		self.InferCountry()
		self.InferKey()
		self.InferSpecies()
		self.InferLat()
		self.InferLong()
		self.InferArea()

	def InferKey(self):
		self.gbif_key = self.record.gbifkey
		return self.gbif_key

	def InferCountry(self):
		self.country = self.record.country
		return self.country

	def InferSpecies(self):
		self.species = self.record.namecomplete
		#T = SPcheck(self.species,'Species')
		#self.species = T.validname
		return self.species

	def InferLat(self):
		if self.record.lat:
			if type(self.record.lat) is float:
				self.lat = self.record.lat
			elif self.record.lat.startswith('Expected'):
				bline = self.record.lat.split()
				try:
					deg  = float(bline[7])
				except:
					deg  = float(bline[7][:-1])
				try:
					minu = float(bline[8][:-1])/60.0
				except:
					minu = 0.0
				try:
					seg  = float(bline[9][:-1])/3600
				except:
					seg = 0.0
				hem  = bline[-1]
				if hem.upper() == 'N':
					lat = deg+minu+seg
				elif hem.upper() == 'S':
					lat = (deg+minu+seg)*-1
				self.lat = lat

		return self.lat

	def InferLong(self):
		if self.record.lon:
			if type(self.record.lon) is float:
				self.lon = self.record.lon
			elif self.record.lon.startswith('Expected'):
				bline = self.record.lon.split()
				try:
					deg  = float(bline[7])
				except:
					deg  = float(bline[7][:-1])
				try:
					minu = float(bline[8][:-1])/60.0
				except:
					minu = 0.0
				try:
					seg  = float(bline[9][:-1])/3600
				except:
					seg = 0.0
				hem  = bline[-1]
				if hem.upper() == 'E':
					lon = deg+minu+seg
				elif hem.upper() == 'W':
					lon = (deg+minu+seg)*-1
				self.lon = lon
		return self.lon

	def InferArea(self):
		self.area = self.record.area

class Species:
	'''
	Summarizes individuals data per species
	'''
	def __init__(self, species='', recs='', polyname=''):
		self.species = species
		self.recs    = recs
		self.countries   = []
		self.individuals = []
		self.lats        = []
		self.longs       = []
		self.latr        = None
		self.lonr        = None
		self.GetIndiv()
		self.GetCountries()
		self.GetLats()
		self.GetLongs()
		self.GetLatRange()
		self.GetLongRange()

	def GetIndiv(self):
		for el in self.recs.obs_recs_list:
			if el.area == polyname:
				if el.namecomplete == self.species:
					i = individual(el)
					if i.lat and i.lon:
						self.individuals.append(i)
			else:
				continue
		return self.individuals

	def GetCountries(self):
		for i in self.individuals:
			self.countries.append(i.country)
		return self.countries

	def GetLats(self):
		for j in self.individuals:
			if j.lat != None:
				self.lats.append(j.lat)
		return self.lats

	def GetLongs(self):
		for k in self.individuals:
			if k.lon != None:
				self.longs.append(k.lon)
		return self.longs

	def GetLatRange(self):
		if self.lats:
			ma  = max(self.lats)
			mi  = min(self.lats)
			med = median(self.lats)
			self.latr = ( mi , med , ma )
		return self.latr

	def GetLongRange(self):
		if self.longs:
			ma  = max(self.longs)
			mi  = min(self.longs)
			med = median(self.longs)
			self.lonr = ( mi , med , ma )
		return self.latr

class SPcheck:
	'''
	Conversion of mike's access to CAS
	'''
	def __init__(self, SEARCH='Urotrygon aspidurus',TYPE='Species'):
		if len(SEARCH.split()) == 2:
			self.SEARCH = SEARCH
		else:
			newname = SEARCH.split()[0] + ' ' + SEARCH.split()[1]
			self.SEARCH = newname
		self.TYPE   = TYPE
		self.data   = None
		self.URL = None
		self.CGI = None
		self.HEADERS = None
		self.connection = None
		self.f = None
		self.validname = ''
		self.SetData()
		self.SetGlobals()
		self.SetConnection()
		self.GetConnectionResult()
		self.GetValidName()
		print 'Probando la validez de los nombres especificos'

	def SetData(self):
		# The POST data includes the search parameters and the submit button (which probably isn't even used)
		self.data = {'tbl':self.TYPE,'contains':self.SEARCH,'Submit':'Search'}
		# Encode POST data
		self.data = urllib.urlencode(self.data)
		return self.data

	def SetGlobals(self):
		self.URL    = "researcharchive.calacademy.org"
		self.CGI = "/research/Ichthyology/catalog/fishcatmain.asp"
		self.HEADERS = {
			"Accept"          : "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8",
			"Accept-Encoding" : "gzip, deflate",
			"Accept-Language" : "en-us,en;q=0.5",
			"Cache-Control"   : "max-age=0",
			"Content-Length"  : str(len(self.data)),
			"Content-Type"    : "application/x-www-form-urlencoded",
			"Host"            : "researcharchive.calacademy.org",
			"Referer"         : "http://researcharchive.calacademy.org/research/Ichthyology/catalog/fishcatmain.asp",
			"User-Agent"      : "Mozilla/5.0 (Macintosh; Intel Mac OS X 10.6; rv:11.0) Gecko/20100101 Firefox/11.0"
		}
		return

	def SetConnection(self):
		# Send POST request and get resulting data
		self.connection = httplib.HTTP(self.URL,80)
		self.connection.putrequest("POST", self.CGI)
		for header,value in self.HEADERS.iteritems(): 
			self.connection.putheader(header,value)
		self.connection.endheaders()
		self.connection.send(self.data)
		reply, msg, hdrs = self.connection.getreply()
		if reply != 200:
			print "ERROR: Website returned error",reply,msg
			sys.exit(0)

		return self.connection

	def GetConnectionResult(self):
		self.f = self.connection.getfile()
		return self.f

	def CurrentStatus(self, HTML, query):
		#  when not available
		#if HTML.find('Not available') != -1:
		#	pass
		#else:
		current = HTML.split('<b>Current status:</b>')
		if current[1].find('Valid as <i>') == -1:
			name = current[1][current[1].find('Synonym of <i>')+12:current[1].find('</i>')]
		else:
			name = current[1][current[1].find('Valid as <i>')+12:current[1].find('</i>')]
		return name	

	def cleanDuplicates(self,names):
		b = False
		count={}
		s = set(names)
		for i in s:
			count[i]=names.count(i)
		if self.SEARCH in s:
			b=True
		return b, count, s

	def GetValidName(self):
		names=[]
		for l in self.f:
			try:
				#if '<p class="result"><b>' in l:
				name = self.CurrentStatus(l,self.SEARCH)
				names.append(name)
			except:
				continue
		valid, count, s = self.cleanDuplicates(names)
		if valid == True:
			self.validname = self.SEARCH
		elif  len(s) == 0:
			pass		

		elif len(s) != 1:
			self.validname = max(count)

		else:
			self.validname = names[0]

def SpList(recs,polyname):
	splist=[]
	checked={}
	for r in recs.obs_recs_list:
		if r.area == polyname:
			#check if only species name reported, if not jut get the friggin first two entries
			if r.namecomplete in checked:
				species = checked[r.namecomplete]
			else:
				if len(r.namecomplete.split()) == 2:
					T = SPcheck(r.namecomplete,'Species')
				else:
					newname=r.namecomplete.split()[0] + ' ' + r.namecomplete.split()[-1]
					T = SPcheck(r.namecomplete,'Species')
				if T.validname == '':
					print 'ATENCION!!!! %s no se encuentra en la base de datos de CAS'\
						  ', se usara este nombre pero debe ser evaluado. Es posible '\
						  'que sea una mala anotacion.'%(r.namecomplete)
					species = r.namecomplete
				else:
					species = T.validname
				checked[r.namecomplete]=species

			if not species in splist:
				splist.append(species)
			for sp in splist:
				if sp == '':
					continue
				else:
					S = Species(species,recs,polyname)
	return splist, S

def distribution(sps,poly, polyname, step=5):
	'''
	giving the dictionary created with instances of the class species and the polygon
	will create a distribution histogram
	'''
	#get the max and mins giving the polygon
	if poly:
		max_lat , min_lat = max([x[1] for x in poly]) , min([x[1] for x in poly])
		max_lon , min_lon = max([x[0] for x in poly]) , min([x[0] for x in poly])
	else:
		max_lat , min_lat = 90 , -90
		max_lon , min_lon = 180 , -180	
	lat_range = range(int(round(min_lat)), int(round(max_lat)),step)
	lon_range = range(int(round(min_lon)), int(round(max_lon)),step)
	# for now lets just code latitud
	labels = ()# [item.get_text() for item in ax.get_yticklabels()]
	count_lat={}
	for i in range(len(lat_range)):
		b=0
		labels+=('%d to %d'%(lat_range[i],lat_range[i]+5),)
		#ax.set_yticklabels('%d to %d'%(lat_range[i],lat_range[i]+5), y=lat_range[i])
		for S in sps.itervalues():
			if S.latr[0] <= lat_range[i] <= S.latr[2]:
				b+=1
		count_lat[lat_range[i]]=(b)
	#plt.hist(count_lat,lat_range,rwidth=0.5)#,orientation='horizontal')
	C = (sorted(count_lat.items(), key=lambda x: x[0]))
	lat = [x[0] for x in C]
	count = tuple([y[1] for y in C])
	fig, ax = plt.subplots()
	index = array(lat)
	bar_width = step
	rects1 = plt.barh(index,count,bar_width)
	plt.xlabel('Frequency')
	plt.ylabel('Latitunidal ranges (degrees)')
	plt.title('Distribution of species in %s'%polyname)
	plt.yticks(index + bar_width/2, labels)
	plt.tight_layout()
	plt.show()
	ax.spines['top'].set_color('none')
	ax.spines['right'].set_color('none')
	ax.tick_params(top='off', right='off')
	plt.show()    

def write_csv(polyname,sps):
	outf = open(polyname+'.csv','w')
	line='Nombre cientifico\tLat desde\tPunto medio\tLat hasta'
	for sp, S in sps.iteritems():
		line+='\n'+sp+'\t'+S.latr[0]+'\t'+S.latr[1]+'\t'+S.latr[2]+'\t'
		for c in sorted(S.countries):
			line+=c+'\t'
	line = line[:-1]
	outf.write(line)

# ########################################################################################################### #
# End definitions ########################################################################################### #
# ########################################################################################################### #

# Aplication of the code#######################################################################################

#Say Hi!!
print 'Bienvenido a BIOGEOGRAPHER. Un script en Python para acceder y manipular las entradas del GBIF'
print

#Get the input from user
kmls=G('*.kml')
if not kmls:
	print 'No hay un archivo del area a buscar. Se buscara en todo el globo.'
	polyname = 'ALL'
	poly = False
# means that a kml file is going to be provided
else:
	if len(kmls) >= 2: 
		print 'more than one file provided working only with %s'%kmls[0]
		poly = kml_to_poly(kmls[0])
	elif len(kmls) == 1:
		poly = kml_to_poly(kmls[0])

	polyname = raw_input('Nombre la region escogida(preferiblemente sin espacios ni acentos): ')

while True:
	try:
		f = raw_input('Nombre del archivo con la lista de nombres especificos a buscar: ')
		f = open(f)
		break
	except:
		print '\t No hay nungun archivo con ese nombre, intentalo de nuevo.'

sps={}
time_queried=dt.datetime.now()
for sn in f:
	sn = sn.strip()
	print 'Procesando %s'%sn
	fn = '_'.join(sn.split())+'.pckl'
	#records to download per server request
	params = {'format': 'darwin', 'scientificname': sn}
	while 1:
		if os.path.isfile(fn):
			last , recs = P.load(open(fn))
			elapsed = time_queried - last
			days = elapsed.days
			if days <= 60:
				print 'La base de datos fue actualizada hace menos de dos meses. Se usara la base de datos local.'
				break
		else:    
			recs = GbifSearchResults()
			gbif_xmltree_list = recs.get_all_records_by_increment(params, 1000)
			P.dump((time_queried,recs),open(fn,'wb'))
			break
	#Set the location of search
	recs.classify_records_into_area(poly,polyname)
	# get the records of the area into species structure
	splist, S = SpList(recs,polyname)
	sps[sn]=S

distribution(sps , poly, polyname, step=5)
write_csv(polyname,sps)
print 'DONE'
#I = individual(recs)