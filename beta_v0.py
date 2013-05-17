#!/usr/bin/env python

import sys
import os
import pyfits as pf
import numpy as np
from optparse import OptionParser

def writeFITS(image):
	hdu = pf.PrimaryHDU(image)

	#write header parameters
	hdu.header.update('BSCALE',   1.000000000000E+00, 'PHYSICAL = PIXEL*BSCALE + BZERO')
        hdu.header.update('BZERO',    0.000000000000E+00)
        hdu.header.update('BTYPE',    'Brightness/Specific Intensity')
        hdu.header.update('SOURCE',   'Simulated skies')
        hdu.header.update('BUNIT',    'JY/BEAM', 'Brightness unit')
	hdu.header.update('EPOCH',    2.000000000000E+03)

	hdu.header.update('CTYPE1',   'RA---SIN')
	hdu.header.update('CRVAL1',   ra * _HRS2DEG)
	hdu.header.update('CDELT1',   cdelt)
	hdu.header.update('CROTA1',   0.000000000000E+00)
	hdu.header.update('CRPIX1',   npix/2.)
	hdu.header.update('CUNIT1',   'deg')

	hdu.header.update('CTYPE2',   'DEC--SIN')
	hdu.header.update('CRVAL2',   dec)
	hdu.header.update('CDELT2',   cdelt)
	hdu.header.update('CROTA2',   0.000000000000E+00)
	hdu.header.update('CRPIX2',   npix/2.)
	hdu.header.update('CUNIT2',   'deg')
	
	hdu.header.update('CTYPE3',   'STOKES  ')
	hdu.header.update('CRVAL3',   1.000000000000E+00)
	hdu.header.update('CDELT3',   1.000000000000E+00)
	hdu.header.update('CROTA3',   0.000000000000E+00)
	hdu.header.update('CRPIX3',   1.000000000000E+00)
	hdu.header.update('CUNIT3',   '        ')
	
	hdu.header.update('CTYPE4',   'FREQ    ')
	hdu.header.update('CRVAL4',   freq)
	hdu.header.update('CDELT4',   chanwid)
	hdu.header.update('CROTA4',   0.000000000000E+00)
	hdu.header.update('CRPIX4',   1.000000000000E+00)
	hdu.header.update('CUNIT4',   'HZ      ')

	if overwrite == True:
		hdu.writeto(filename,clobber=True)
	else:
		try:
			hdu.writeto(filename);
		except IOError as e:
			print e, '\nExiting without creating/overwriting FITS file.'
	return;

if __name__ == '__main__':

	parser = OptionParser(usage="""%prog [options]""",
        description="Create a FITS file with gas distributed according to the 2D-beta model");
        parser.add_option("-f","--fits",type="string",metavar="FITSFILE",default="beta.fits",
                          dest="file",help="Name of output FITS file (def: beta.fits)");
        parser.add_option("-n","--npix",type="int",metavar="NPIX",default=4096,
	                  dest="npix",help="Image size assumed to be npix X npix (def: 4096)");
        parser.add_option("-i","--imsize",type="float",metavar="IMSIZE",default=10.0,
	                  dest="imsize",help="Image size in arcminutes (def: 10.0)");
        parser.add_option("--ra",type="float",metavar="RA",default=0.0,
	                  dest="ra",help="Right Ascension of the phase centre in hours (def: 0.0)");
        parser.add_option("--dec",type="float",metavar="DEC",default=0.0,
	                  dest="dec",help="Declination of the phase centre in degrees (def: 0.0)");
        parser.add_option("-c","--rcore",type="int",metavar="CORE-R",default=None,
	                  dest="rcore",help="Radius of the cluster core (def: 1/5th of image size)");
        parser.add_option("-e","--ellip",type="float",metavar="ELLIPTICITY",default=0.0,
	                  dest="ellip",help="Ellipticity parameter (def: 0.0)");
        parser.add_option("-t","--theta",type="float",metavar="THETA",default=0.0,
	                  dest="theta",help="Ellipse orientation angle in degrees (def: 0.0)");
        parser.add_option("-s","--sjy",type="float",metavar="FLUXDENS",default=1.0,
	                  dest="sjy",help="Source flux density in Jansky (def: 1.0)");
        parser.add_option("-b","--beta",type="float",metavar="BETA",default=0.2,
	                  dest="beta",help="beta value in exponent (def: 0.2)");
        parser.add_option("-v","--freq",type="float",metavar="FREQ",default=1.4E9,
	                  dest="freq",help="Channel frequency in Hz (def: 1.4E9)");
        parser.add_option("--chanwid",type="float",metavar="CHANWIDTH",default=3.90625E5,
	                  dest="chanwid",help="Channel width in Hz (def: 3.90625E5 - KAT-7 wideband mode)");
        parser.add_option("-o","--overwrite",action="store_true",dest="overwrite",
	                  help="Overwrite FITS file if it exists");
	parser.set_defaults(overwrite=False);

	(options,args) = parser.parse_args()
	
	if(len(args)) != 0:
		parser.print_help();
		sys.exit(0);

	#Define necessary constants
	_DEG2RAD = np.pi /180.0
	_ARCMIN2DEG = 1.0/60.0
	_HRS2DEG = 15.0

	#Assign values to variables for ease of use.
	filename = options.file
	npix = options.npix
	imsize = options.imsize
	ra = options.ra
	dec = options.dec
	max_r = np.sqrt(2*np.power((npix/2.0),2));
	if options.rcore != None and options.rcore <= max_r:
		rcore = options.rcore
	else:
		rcore = max_r/5.0 # Cluster core size - rule of thumb 1/5ths of the cluster halo extent.
	ellip = options.ellip
	theta = options.theta * _DEG2RAD
	sjy = options.sjy
	beta = options.beta
	freq = options.freq
	chanwid = options.chanwid
	overwrite = options.overwrite

	cdelt = (imsize/npix)*_ARCMIN2DEG

	#Create gas distribution that follows the 2-D beta model (Lorentz model)
	l_arr,m_arr = np.indices((npix,npix));
	l_arr = l_arr * (- 1) + npix / 2
	m_arr = m_arr - npix / 2

	l_new = np.abs(l_arr)*np.cos(theta) + np.abs(m_arr)*np.sin(theta)
	m_new = np.abs(m_arr)*np.cos(theta) - np.abs(l_arr)*np.sin(theta)
	r = (np.sqrt(np.power(l_arr,2)*np.power((1.0-ellip),2) + np.power(m_arr,2)))/(1.0 - ellip)

	image = np.zeros(shape=(1,1,npix,npix))
	betahalo = sjy*np.power((1 + np.power((r/rcore),2)),-beta)

	#Write to FITS file
	image[0,0,:,:] = betahalo
	writeFITS(image);
