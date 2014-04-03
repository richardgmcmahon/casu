#!/usr/bin/python
"""
Query the DQC and writes a FITS file with the result:
  python querydqc.py output.fits



 VISTA Program IDs
 VHS:    179.A-2010 
 VIKING: 179.A-2004

 VST
 ATLAS: 177.A-3011

Requirements:
  It needs pyfits + numarray installed as well as the DateTime module from mx.
  It uses psycopg to comunicate with the postgres db.

Port forwarding Access 
dbi_vista="host=apm45.ast.cam.ac.uk dbname=vista user=vuser 
 password=dqc_user port=5432"

On a muon e.g. muon2 type: ssh -L 5432:apm45:5432 calx042

On local machine: ssh -L 5432:127.0.0.1:5432 muon2

On local machine:
dbi_vista="host=127.0.0.1 dbname=vista user=vuser password=dqc_user port=5432"


Caveats:
  The column order in the output file is random.


History:
  20060413 - Initial write. EAGS

  20060415 - Fix the case where columns in database are NULL. EAGS.
  20060414 - Updated to work with psycopg2. Check length of strings 
             from database information. Use database information to 
             get column types. EAGS

  20090316 - EAGS: Removed dependency on mx DateTime module
  20090731 - RGM: converted from numarray to numpy

  20101009 - RGM: added naxis1, naxis2 to query
  20100107 - EGS: bugfix and creation of vistaqc view
                  SELECT * FROM vistaqc

  20130501 - RGM: added verbose info and python library versions

  20130519 - RGM: added support for config file to store connection info

"""

import os
import sys
import types
import datetime

from time import strftime, gmtime

# FITS module
import pyfits

from UserDict import UserDict

#import numpy as np

from numpy import *

#except:
#  from numarray import *
#  from numarray.records import chararray

# Import Postgres module
try:
    import psycopg2 as psycopg
except:
    import psycopg

import ConfigParser
config = ConfigParser.RawConfigParser()
config.read('querydqc.cfg')

db='wfcam'
db='vista'
db='vst'


host = config.get(db,'host')
dbname = config.get(db,'dbname')
user = config.get(db,'user')
password = config.get(db,'password')
port = config.get(db,'port')

verbose=True
if verbose:
   print 'Python version:     ', sys.version
   print 'pyfits.__version__: ', pyfits.__version__
   print 'numpy.__version__:  ', __version__
   print 'psycopg version:    ', psycopg.__version__

class odict(UserDict):
    def __init__(self, dict = None):
        self._keys = []
        UserDict.__init__(self, dict)

    def __delitem__(self, key):
        UserDict.__delitem__(self, key)
        self._keys.remove(key)

    def __setitem__(self, key, item):
        UserDict.__setitem__(self, key, item)
        if key not in self._keys: self._keys.append(key)

    def clear(self):
        UserDict.clear(self)
        self._keys = []

    def copy(self):
        dict = UserDict.copy(self)
        dict._keys = self._keys[:]
        return dict

    def items(self):
        return zip(self._keys, self.values())

    def keys(self):
        return self._keys

    def popitem(self):
        try:
            key = self._keys[-1]
        except IndexError:
            raise KeyError('dictionary is empty')
        val = self[key]
        del self[key]
        return (key, val)

    def setdefault(self, key, failobj = None):
        UserDict.setdefault(self, key, failobj)
        if key not in self._keys: self._keys.append(key)

    def update(self, dict):
        UserDict.update(self, dict)
        for key in dict.keys():
            if key not in self._keys: self._keys.append(key)

    def values(self):
        return map(self.get, self._keys)

# Output file
output = sys.argv[1]

logdata='Current working directory: %s' % (os.getcwd())
print logdata

logdata='Executing: %s' % sys.argv[0]
print logdata

logdata = 'Owner: '+ os.getenv('USER') + '@' + os.getenv('HOST')
print logdata

# DQC connection string
#dbi_wfcam="host=apm25.ast.cam.ac.uk dbname=main-wfcam-jdqc user=wfcam password=wfcam_dqc port=5433"

#dbi_vista="host=apm25.ast.cam.ac.uk dbname=dqc user=vuser password=dqc_user port=5437"

#dbi_vista="host=apm25.ast.cam.ac.uk dbname=dqc user=vuser password=dqc_user port=5438"

#dbi_vista="host=apm45.ast.cam.ac.uk dbname=vista user=vuser password=dqc_user port=5432"

#host = config.get(db,'host')
#dbname = config.get(db,'dbname')
#user = config.get(db,'user')
#password = config.get(db,'password')
#port = config.get(db,'port')

dbi="host=" + host + " dbname=" + dbname + " user=" + user + " password=" + password +  " port=" + port

print 'dbi: ', dbi

#dbi_vista="host=" + host + " dbname=" + dbname + " user=" + user + " password=" + password +  " port=" + port

#print 'dbi_vista: ', dbi_vista

#dbi_vst="host=apm45.ast.cam.ac.uk dbname=vst user=vuser password=dqc_user port=5432"



#dbi=dbi_wfcam
#dbi=dbi_vista
#dbi=dbi_vst

"""
# Query. We could read this from a file. These tests could be added
as a dictionary and could be used as regression tests
"""

query_ukidss="""
SELECT
  * 
FROM
  Main.Mv_chips_las
LIMIT 10 
"""

query_ukidss_las="""
SELECT 
  chip_id,field_id,survey_id,stacked_id,
  tcenra, tcendec,
  fcenra,fcendec,
  filtname, chip,
  filename, catalogue,
  mjdobs, dateobs,
  telra, teldec, 
  cen_ra, cen_dec,
  airmass,
  ra_lo, ra_hi, dec_lo, dec_hi,
  moonfli, moonra, moondec,
  ndetect,
  magzpt,magzerr,extinct,
  magzptcor,apcor3,percorr,
  cnt_skylev, cnt_skyrms, mag_skylev, mag_skyrms,
  pix_seeing, sec_seeing, ellipt, maglimit,
  project, object, tile, field,
  naxis1, naxis2,
  ctype1, ctype2,
  crval1, crval2,
  crpix1, crpix2,
  cd11, cd12, cd21, cd22, pv21, pv22, pv23,
  exptime, nexp, nint, njitter, nmstep, readmode,
  wcsrms, numwcsfit,
  redpath
FROM
  mv.v_chips_las
"""


query_ukidss_gsc="""
SELECT 
  chip_id,field_id,survey_id,stacked_id,
  tcenra, tcendec,
  fcenra,fcendec,
  filtname, chip,
  filename, catalogue,
  mjdobs, dateobs,
  telra, teldec, 
  cen_ra, cen_dec,
  airmass,
  ra_lo, ra_hi, dec_lo, dec_hi,
  moonfli, moonra, moondec,
  ndetect,
  magzpt,magzerr,extinct,
  magzptcor,apcor3,percorr,
  cnt_skylev, cnt_skyrms, mag_skylev, mag_skyrms,
  pix_seeing, sec_seeing, ellipt, maglimit,
  project, object, tile, field,
  naxis1, naxis2,
  ctype1, ctype2,
  crval1, crval2,
  crpix1, crpix2,
  cd11, cd12, cd21, cd22, pv21, pv22, pv23,
  exptime, nexp, nint, njitter, nmstep, readmode,
  wcsrms, numwcsfit,
  redpath
FROM
  mv.v_chips_gsc
"""

query1="""
SELECT
  table_name
FROM
  information_schema.tables
WHERE
  table_schema='public'
"""

query2="""
SELECT
  *
FROM
  qc
"""

query3="""
SELECT
  *
FROM
  extension
"""

query4="""
SELECT
  *
FROM
  image
"""

query5="""
SELECT 
  Image.id AS "Image.id", Extension.id AS "Extension.id", 
  Extension.image_id AS "Extension.image_id", 
  Qc.id AS "QC.id" , 
  Qc.extension_id AS "Qc.extension_id", 
  extension.chipno AS "extension.chipno"
FROM
  Qc,
  Image,
  Extension
WHERE
 Image.id = Extension.image_id 
 AND Extension.id = Qc.extension_id
ORDER BY
  Image.id, Extension_id
"""

query6="""
SELECT 
  Image.*, Extension.*, Qc.*
FROM
  Qc,
  Image,
  Extension
WHERE
 Image.id = Extension.image_id AND
 Extension.id = Qc.extension_id
ORDER BY Image.id, Extension.id
"""

query7="""
SELECT 
  COUNT(*) 
FROM
  Qc,
  Image,
  Extension
WHERE
 Image.id = Extension.image_id AND
 Extension.id = Qc.extension_id
"""

query8="""
SELECT
  image.*, filter.*, programme.*, extension.*, qc.*
FROM 
  image, extension, qc, programme, filter
WHERE 
  image.id = extension.image_id
  and extension.id = qc.extension_id
  and image.programme_id = programme.id
  and image.insfilter_id = filter.id
"""


# view created 20100107 by EGS
query9="""
SELECT 
  * 
FROM 
  vistaqc 
WHERE 
  is_tile = 'True' 
LIMIT 1000
"""


query_all_tiles_thin="""
SELECT 
 image_id, version, 
 ra, dec, prog, filtname,
 mjd, dateobs, nightobs,
 obstatus, qcstatus,
 is_tile, is_stack
FROM 
  vistaqc 
WHERE 
  is_tile = 'True'
"""

query_all_pawprints_thin="""
SELECT 
 image_id, version, 
 ra, dec, prog, filtname,
 mjd, dateobs, nightobs,
 obstatus, qcstatus,
 is_tile, is_stack
FROM 
  vistaqc 
WHERE 
  is_tile = 'False' AND
  is_stack = 'True' 
"""




# VHS Pawprints
query10="""
SELECT 
  * 
FROM 
  vistaqc 
WHERE 
  is_tile <> 'True' 
  AND prog = '179.A-2010'
"""

# VHS Pawprints or Tiles
query10a="""
SELECT 
 *
FROM 
  vistaqc 
WHERE 
  is_tile = 'True' 
  /* is_tile <> 'True' 
  AND is_stack = 'True' 
  */
  /* VHS */
  AND prog = '179.A-2010'
"""

# VHS Pawprints/Tiles with selected columns
query10b="""
SELECT 
 image_id, version, 
 ra, dec, prog, filtname,
 filename, filepath,
 exptime, image_exptime, totexptime, 
 reqtime, ocsreqtime,
 mjd, dateobs, nightobs, obsstart, obsid, obstplno,  
 obsname, irunno,
 obstatus, qcstatus,
 amstart, amend, is_tile, is_stack, 
 naxis1, naxis2, chipno, 
 cenra, cendec,
 nobjects,
 pixsize, 
 magzpt, magzerr, magznpt, 
 ellipticity, ellpt_min, ellpt_max,
 seeing, 
 skynoise, skylevel,
 maglim,  
 stdcrms, numbrms,
 apcor1,apcor2, apcor3, apcor4, apcor5,
 offsetx, offsety,
 jitterx, jittery
FROM 
  vistaqc 
WHERE 
  is_tile = 'True' 
  /* is_tile <> 'True' */
  /* AND is_stack = 'True' */
  /* VHS */
  AND prog = '179.A-2010'
"""


query10c="""
/* VHS Tiles */
SELECT 
 *
FROM 
  vistaqc 
WHERE 
  is_tile = 'True' 
  /* is_tile <> 'True' */
  /* AND is_stack = 'True' */
  /* VHS */
  AND prog = '179.A-2010'
"""

query10d="""
/* All Tiles */
SELECT 
 *
FROM 
  vistaqc 
WHERE 
  is_tile = 'True' 
"""

# VIKING Pawprints
query11a="""
SELECT 
  * 
FROM 
  vistaqc 
WHERE 
  is_tile <> 'True' 
  AND is_stack = 'True' 
  /* VIKING */
  AND prog = '179.A-2004'
LIMIT 1000
"""

# VIKING Pawprints with selected columns
query11b="""
SELECT 
 image_id, version,
 ra, dec, prog, filtname,
 filename, filepath,
 exptime, image_exptime, totexptime, 
 reqtime, ocsreqtime,
 mjd, dateobs, nightobs, obsstart, obsid, obstplno,  
 obsname, irunno,
 obstatus, qcstatus,
 amstart, amend, is_tile, is_stack, 
 naxis1, naxis2, chipno, 
 cenra, cendec,
 nobjects,
 pixsize, 
 magzpt, magzerr, magznpt, 
 ellipticity, 
 seeing, 
 skynoise, skylevel,
 maglim, 
 stdcrms, numbrms,
 apcor1,apcor2, apcor3, apcor4, apcor5,
 offsetx, offsety,
 jitterx, jittery
FROM 
  vistaqc 
WHERE 
  is_tile <> 'True' 
  AND is_stack = 'True' 
  /* VIKING */
  AND prog = '179.A-2004'
"""

query12a="""
SELECT 
  * 
FROM 
  vstqc
WHERE
  prog = '177.A-3011'
  /* AND is_stack = 'True' */
  AND chipno = 1
  /* LIMIT 100000 */
"""

query12a="""
SELECT 
  is_stack, is_tile, naxis1, naxis2, version, chipno, dateobs, nightobs
FROM 
  vstqc
WHERE
  prog = '177.A-3011'
  /* AND is_stack = 'True' */
  /* AND chipno = 1  */
  /* LIMIT 100000 */
"""



query12b="""
SELECT 
 image_id, version,
 ra, dec, prog, filtname,
 filename, filepath,
 exptime, 
 mjd, dateobs, nightobs, obsstart, obsid, obstplno,  
 obsname, irunno,
 obstatus, qcstatus,
 amstart, amend, is_tile, is_stack, 
 naxis1, naxis2, chipno, 
 cenra, cendec,
 nobjects,
 pixsize, 
 magzpt, magzerr, magznpt, 
 ellipticity, 
 seeing, 
 skynoise, skylevel,
 maglim, 
 stdcrms, numbrms,
 apcor1,apcor2, apcor3, apcor4, apcor5
 FROM 
  vstqc 
"""

month = "201208"
progid  =  "179.A-2010"

#print "Program ID: " + progid + " ;" + month
print "Program ID: %s ;%s" % (progid, month)
#print """Program ID: %s ;%s""" % (progid, month)



query201208="""
SELECT 
  * 
FROM 
  vistaqc 
WHERE
  filename like '%201208%'
  /* VHS */
  AND prog = '179.A-2010'
""" 
#% (month, progid)

query=query201208
query=query10a
query=query_all_tiles_thin
query=query_all_pawprints_thin
#query=query11b

# vst
query=query12a

print 'sql: ', query

pause=raw_input("Type return to continue >> ")

#query=query_ukidss_las
#query=query_ukidss_gsc

# Connect to database
con=psycopg.connect(dbi)

# Get a cursor and execute query
cur=con.cursor()
cur.execute(query)

#help(cur)
print 'cur.rowcount(): ', cur.rowcount


# Return all results in a dictionary and information about columns.
# This is a bit obfuscated.
values = cur.fetchall()
#help(values)
print 'len(values): ', len(values)

keys = [r[0] for r in cur.description]     # Column name
ct = [r[1] for r in cur.description]       # Column type
ml = [r[2] for r in cur.description]       # Column length (used for strings)
dd=odict()
ctype={}
maxlength={}

for k,v in zip(keys,ct): ctype[k]=v

for k,v in zip(keys,ml):
    if v==None:
    	v=80
    maxlength[k]=v

#for k in keys: dd[k]=[]

for value in values:
    kdone=[]
    for k,v in zip(keys, value):
        if k not in kdone:
            if k not in dd.keys():
                dd[k]=[]
            dd[k].append(v)
            kdone.append(k)
        else:
            for i in range(10): # Need this for duplicate keys in the query. We append _1, _2, etc.
                newk = '%s_%d' % (k, i+1)
                if not (newk in dd.keys()):
                    dd[newk]=[]

                if not (newk in kdone):
                    dd[newk].append(v)
                    ctype[newk]=ctype[k]
                    maxlength[newk]=maxlength[k]
                    kdone.append(newk)
                    break

# Get column types
fmt=[]
keys = dd.keys()

for k in keys:
    if ctype[k] in [700, 701]:                         # Float/Double

        fmt.append('D')

    elif ctype[k] in [19, 18, 25, 1042, 1043]:         # String

        fmt.append('%sA' % maxlength[k])

    elif ctype[k] == 16:                               # Boolean

        fmt.append('L')

    elif ctype[k] in [20, 21, 23, 1700]:                     # Integer

        fmt.append('J')

    elif ctype[k] in [1114, 1184, 704, 1186, 1082]:    # Date

        fmt.append('J')

        for j in range(len(dd[k])):

            try:

                dd[k][j] = int(dd[k][j].strftime('%Y%m%d'))

            except:

                dd[k][j]=0

        #dd[k] = [int(r.strftime('%Y%m%d')) for r in dd[k]]

    else:

        print 'Column %s: Type not recognized (%s)' % (k, ctype[k])

        fmt.append('')



# We get the column types from the first result

cols=[]

icol=0

for k,f in zip(keys, fmt):

    icol=icol+1

    if f[-1]=='A': # If column type is string then we need a special case

        print 'Parsing format: A'
        try:

            # Some columns have None (null) as value. Convert them to zero.

            indx=[i for i in range(len(dd[k])) if dd[k][i]==None]

            for i in indx: dd[k][i]=''

            # Now create the column

            #c = pyfits.Column(name=k, format=f, array=pyfits.chararray.array(dd[k]))
            c = pyfits.Column(name=k, format=f, array=pyfits.np.array(dd[k]))

            #print 'test1',icol, len(dd[k]), k, f, '0 0'

            #print 'test2',icol, len(dd[k]), k, f, min(c), max(c)

            print icol, len(dd[k]), k, f, min(dd[k]), max(dd[k])

            cols.append(c)

        except:

            print '%s Error in column %s (%s)' % (icol, k, f)

    elif f<>'':                 # Otherwise the general case

        try:

            # Some columns have None (null) as value. Convert them to zero.

            indx=[i for i in range(len(dd[k])) if dd[k][i]==None]

            for i in indx: dd[k][i]=0

            # Now create the column

            data = array(dd[k])

            c = pyfits.Column(name=k, format=f, array=data)

            print icol, len(dd[k]), k, f, data.min(), data.max()

            print icol, len(dd[k]), k, f, min(data), max(data)

            cols.append(c)

        except:

            print '%s Error in column %s (%s)' % (icol, k, f)



# Create table

coldefs = pyfits.ColDefs(cols)

tabhdu = pyfits.new_table(coldefs)



# Write table. Fails if already exists.

tabhdu.writeto(output)





