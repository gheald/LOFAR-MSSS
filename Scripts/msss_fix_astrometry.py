#!/usr/bin/env python

import argparse
import os
from astropy.table import Table
from numpy import *
from pylab import *
import glob
import pyfits
import datetime

def main(args):
	files = []
	if args.corrfiles != '__none__':
		for pattern in args.corrfiles.split(','):
			files += glob.glob(pattern)
	#print files
	if args.outlist == '__none__':
		wolf = False
	else:
		wolf = True
		olf = open(args.outlist,'a')
	cmd = '/home/heald/stilts tskymatch2 dec1="DEC(2000)" dec2=DEC error=25 ifmt2=ascii in1=%s in2=%s ofmt=votable-tabledata out=%s ra1="RA(2000)" ra2=RA'%(args.reftable,args.intable,args.outtable)
	os.system(cmd)
	t = Table.read(args.outtable)
	ra1 = t['RA_2000_']
	ra2 = t['RA']
	dec1 = t['DEC_2000_']
	dec2 = t['DEC']
	dra=ra1-ra2
	dra[dra>359.]-=360.
	dra[dra<-359.]+=360.
	ddec=dec1-dec2
	totnum = len(dra)
	raoffset=mean(dra)
	decoffset=mean(ddec)
	raoffset_std=std(dra)
	decoffset_std=std(ddec)
	meanoffset=sqrt((dra-raoffset)**2+(ddec-decoffset)**2)
	offseterr=mean(sqrt(raoffset_std**2+decoffset_std**2))
	cliprad = 1.5*offseterr
	print 'Clipping sources that matched beyond %.2f arcsec of mean'%(cliprad*3600.)
	remove=where(meanoffset>cliprad)
	print 'Removing',len(remove[0]),'sources with large offset'
	t.remove_rows(remove)
	totflux = t['Total_flux']
	pkflux = t['Peak_flux']
	fluxratio_mean = mean(totflux/pkflux)
	fluxratio_std = std(totflux/pkflux)
	print 'Mean flux ratio is',fluxratio_mean,'with std',fluxratio_std
	remove = where(totflux/pkflux > fluxratio_mean+1.*fluxratio_std)
	print 'Removing',len(remove[0]),'fluffy sources'
	t.remove_rows(remove)
	ra1 = t['RA_2000_']
	ra2 = t['RA']
	dec1 = t['DEC_2000_']
	dec2 = t['DEC']
	dra=ra1-ra2
	dra[dra>359.]-=360.
	dra[dra<-359.]+=360.
	ddec=dec1-dec2
	num = len(dra)
	print num,'of',totnum,'crossmatches used'
	raoffset=mean(dra)
	decoffset=mean(ddec)
	raoffset_error=std(dra)*3600./sqrt(float(num)) # for print only
	decoffset_error=std(ddec)*3600./sqrt(float(num)) # for print only
	print 'RA offset  (arcsec):',raoffset*3600.,'+/-',raoffset_error
	print 'Dec offset (arcsec):',decoffset*3600.,'+/-',decoffset_error
	# Difference is NVSS - MSSS
	# So, correction is additive
	if len(files) > 0:
		print 'Correcting %d files'%len(files)
		for f in files:
			print f
			if args.inplace:
				h = pyfits.open(f,mode='update')
			else:
				h = pyfits.open(f) # not in update mode
				newf = f.split('/')[-1].split('.fits')[0]+'.ast.fits'
				print 'Will write to',newf
				if wolf: print >>olf, newf
			hdr = h[0].header
			hdr['CRVAL1'] += raoffset
			hdr['CRVAL2'] += decoffset
			hdr.set('ASTR_FIX',datetime.datetime.now().strftime("%Y-%m-%d_%H:%M:%S"))
			hdr.set('RA_SHIFT',raoffset)
			hdr.set('DECSHIFT',decoffset)
			if args.inplace:
				h.close()
			else:
				h[0].header=hdr
				h.writeto(newf,clobber=True)
	else:
		print 'No FITS files provided, so no corrections are being applied.'
	if args.showplot or (args.plotfile != '__none__'):
		axhline(0.,ls='-',c='k',alpha=0.5)
		axvline(0.,ls='-',c='k',alpha=0.5)
		plot(dra*3600.,ddec*3600.,'ko',alpha=0.3,ms=10)
		axhline(decoffset*3600.,ls='--',c='r')
		axvline(raoffset*3600.,ls='--',c='r')
		xlabel('RA offset (arcsec)')
		ylabel('Dec offset (arcsec)')
		if args.plotfile != '__none__': savefig(args.plotfile,bbox_inches='tight')
		if args.showplot: show()

ap = argparse.ArgumentParser()
ap.add_argument('intable',help='Table to match for astrometry')
ap.add_argument('--reftable','-r',help='Path to table to use as reference',default='/home/heald/NVSS/CATALOG41.FIT')
ap.add_argument('--outtable','-o',help='Output table filename',default='./xmatch_astrometry.vot')
ap.add_argument('--showplot','-s',help='Show offset plot? [default False]',default=False,action='store_true')
ap.add_argument('--plotfile','-p',help='Plot figure file, can be used together with --showplot [default none]',default='__none__')
ap.add_argument('--corrfiles','-c',help='FITS files to correct, which can be one or more filename pattern(s) separated by commas [default none]',default='__none__')
ap.add_argument('--inplace','-i',help='Correct FITS files in place? If not, new FITS files will be written in the cwd that have ".ast.fits" as the new suffix [default False]',default=False,action='store_true')
ap.add_argument('--outlist','-l',help='If inplace is false, this specifies the file name that will collect the output fits file names [default none]',default='__none__')
args = ap.parse_args()
main(args)

