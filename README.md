Biogeography
============

Scripts for Biogeography

In this repo you will find some of the scripts I have made to work with biogeography-like data'

biogeographer.py:
----------------
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

Requires: 
- Biopython
- xml.etree
- matplotlib and dependencies

GBIFer.py:
----------
This a program that given a GBIF zip dataset(s) will extract the taxonomic and geographic info
to be analyzed. This program supercedes biogeographer, since the GBIF API change. This script is in development and will contain more analytical tools.

Requires: 
- numpy
- matplotlib and dependencies


-------------
Disclaimer: This code has probably been superceded by GenGIS software plugins (https://github.com/beiko-lab/gengis), however the development will continue since is possible to include pythoncode into GenGIS
