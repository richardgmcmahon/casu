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
password=???? port=5432"

On a muon e.g. muon2 type: ssh -L 5432:apm45:5432 calx042

On local machine: ssh -L 5432:127.0.0.1:5432 muon2

On local machine:
dbi_vista="host=127.0.0.1 dbname=vista user=vuser password=??? port=5432"


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


from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
from builtins import input
try:
   input = raw_input
except NameError:
   pass

import os
import sys
import socket
import types
import datetime

import astropy
import pandas as pd

from time import strftime, gmtime

# FITS module
#import pyfits



from astropy.table import Table
from astropy.io import fits as pyfits


import numpy as np

from tablestats import *

#except:
#  from numarray import *
#  from numarray.records import chararray

# Import Postgres module
try:
    import psycopg2 as psycopg
except:
    import psycopg

import sqlutilpy
print('sqlutilpy.__version__', sqlutilpy.__version__)


def getargs(verbose=False):
    """

    Template getargs function

    Usage

    python getargs.py --help


    def getargs():

    ....

    if __name__=='__main__':

        args = getargs()
        debug = args.debug()



    parse command line arguements

    not all args are active

    """
    import sys
    import pprint
    import argparse

    # there is probably a version function out there
    __version__ = '0.1'

    description = 'Search CASU DQC database'
    epilog = """WARNING: Not all options may be supported
             """
    parser = argparse.ArgumentParser(
        description=description, epilog=epilog,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)


    default_db='vista'
    default_db='vst'
    default_outfile='/tmp/tmp.fits'
    default_sqlfile='querydqc.sql'

    parser.add_argument("-d", "--db", dest="db",
                        default=default_db,
                        help="name of CASU DQC database [wfcam, vista, vst]")

    parser.add_argument("-o", "--outfile", dest="outfile",
                        default=default_outfile,
                        help="output filename")

    parser.add_argument("-q", "--query", dest="sqlfile",
                        default=default_sqlfile,
                        help="SQL query filename")

    # default type is string
    parser.add_argument("--string",
                        help="string input")

    parser.add_argument("--float", type=float,
                        help="float input")

    parser.add_argument("--configfile",
                        default=None,
                        help="configuration file")

    parser.add_argument("--infile",
                        help="Input file name")

    parser.add_argument("--debug",
                        action='store_true',
                        help="debug option")

    parser.add_argument("--pause",
                        action='store_true',
                        help="Pausing option")

    parser.add_argument("--verbose", default=verbose,
                        action='store_true',
                        help="Verbose option")

    parser.add_argument("--version",
                        action='store_true',
                        help="verbose option")

    args = parser.parse_args()

    if args.debug or args.verbose:
        print()
        print('Number of arguments:', len(sys.argv),
              'arguments: ', sys.argv[0])

    if args.debug or args.verbose:
        print()
        for arg in vars(args):
            print(arg, getattr(args, arg))
        print()

    if args.debug or args.verbose:
        print()
        pprint.pprint(args)
        print()

    if args.version:
        print('version:', __version__)
        sys.exit(0)

    return args


def getconfig(configfile=None, debug=False, silent=False):
    """
    read config file

    config = getconfig()
    parameter_value = config.get('INPUTS', 'parameter')

    Note the Python 2 ConfigParser module has been renamed to configparser
    in Python 3 so it better to use import configparser in Python 2 for
    future proofing

    see also getconfig.cfg

    TODO: catch exceptions

    Support for lists:

    see:

    https://stackoverflow.com/questions/335695/lists-in-configparser

    https://github.com/cacois/python-configparser-examples

    look in cwd, home and home/.config

    home/.config not implemented yet


    """
    import os
    import configparser

    # read the configuration file
    # config = configparser.RawConfigParser()
    config = configparser.SafeConfigParser()

    print('__file__', __file__)
    configfile_default = os.path.splitext(__file__)[0] + '.cfg'
    print('configfile_default:', configfile_default)
    if configfile is None:
        configfile_default = os.path.splitext(__file__)[0] + '.cfg'
        configfile = configfile_default
        if debug:
           print('__file__', __file__)
           print('configfile_default:', configfile_default)

    print('Open configfile:', configfile)
    if debug:
        print('Open configfile:', configfile)

    try:
        if not silent:
            print('Reading config file', configfile)

        try:
            config.read(configfile)
        except IOError:
            print('config file', configfile, "does not exist")
            configfile = os.path.join(os.environ["HOME"], configfile)
            print('trying ', configfile)
            config.read(configfile)

    except Exception as e:
        print('Problem reading config file: ', configfile)
        print(e)

    if debug:
        print('configfile:', configfile)
        print('sections:', config.sections())
        for section_name in config.sections():
            print('Section:', section_name)
            print('Options:', config.options(section_name))
            for name, value in config.items(section_name):
                print('  %s = %s' % (name, value))
        print()

        for section_name in config.sections():
            print()
            print('Section:', section_name)
            for name, value in config.items(section_name):
                print('  %s = %s' % (name, value))
                print(section_name, ':',
                      name, config.get(section_name, name))
        print()

    return config




def main(query=None, outfile=None):

    import time

    t0 = time.time()

    logdata = 'Current working directory: %s' % (os.getcwd())
    print(logdata)

    logdata='Executing: %s' % sys.argv[0]
    print(logdata)

    print('Owner: '+ os.getenv('USER'))
    print('Host: '+ socket.gethostname())
    logdata = 'Owner: ' + os.getenv('USER') + '@' + socket.gethostname()
    print(logdata)

    # DQC connection string
    #dbi_wfcam="host=apm25.ast.cam.ac.uk dbname=main-wfcam-jdqc user=wfcam password=wfcam_dqc port=5433"

    #dbi_vista="host=apm25.ast.cam.ac.uk dbname=dqc user=vuser password=dqc_user port=5437"

    #dbi_vista="host=apm25.ast.cam.ac.uk dbname=dqc user=vuser password=dqc_user port=5438"

    #dbi_vista="host=apm45.ast.cam.ac.uk dbname=vista user=vuser password=dqc_user port=5432"

    dbi = "host=" + host + " dbname=" + dbname + " user=" + user + " password=" + password +  " port=" + port

    print('dbi: ', dbi)

    # get the query
    sql_query = sqlquery()
    print('sql: ', sql_query)
    pause = input("Type return to continue >> ")

    print('astropy.__version__', astropy.__version__)
    print('pandas.__version__', pd.__version__)
    print('numpy.__version__', np.__version__)
    print('psycopg.__version__', psycopg.__version__)


    print('Elapsed time(secs):', time.time() - t0)

    UseSqlutil = True
    UsePandas = False

    if UseSqlutil:
        print('sqlutilpy.__version__', sqlutilpy.__version__)
        result = sqlutilpy.sqlutil.get(sql_query,
                                       db=db,
                                       host=host,
                                       user=user,
                                       password=password,
                                       asDict=True)

        table = Table(result)
        table.info()
        print()

        for column in table.itercols():
            try:
                if column.dtype == 'datetime64[us]':
                    print(column.name, column.dtype,
                      table[column.name][0],
                      table[column.name][-1])
                    print('using str32')
                    NewColumn = Table.Column(column.data,
                                             dtype='S32')
                    print('using str32')
                    table.replace_column(column.name, NewColumn)
                    print(column.name, column.dtype,
                      table[column.name][0],
                      table[column.name][-1])
                    print()
            except:
                pass

        table.info()
        print()

    if UsePandas:
        # Connect to database
        conn = psycopg.connect(dbi)

        # Get a cursor and execute query
        cursor = conn.cursor()
        cursor.execute(sql_query)
        print('Elapsed time(secs):', time.time() - t0)

        print('sql to pandas df')
        df = pd.read_sql(sql_query, conn)
        print('pandas:', df.shape)
        print(df.dtypes)

        print('Elapsed time(secs):', time.time() - t0)

        print(df.info(verbose=True))
        print(df.describe(include='all'))

        print('Elapsed time(secs):', time.time() - t0)

        table = Table.from_pandas(df)

        table.info()

        table.info('stats')
        print('Elapsed time(secs):', time.time() - t0)


    table.info('stats')
    print('astropy.__version__', astropy.__version__)
    print('pandas.__version__', pd.__version__)
    print('numpy.__version__', np.__version__)

    tablestats(table)
    print('Output file:', outfile)
    table.write(outfile, overwrite=True)
    print('Elapsed time(secs):', time.time() - t0)

    sys.exit()


    for column in table.itercols():
        try:
            if column.dtype == 'object':
                print(column.name, column.dtype,
                      table[column.name][0],
                      table[column.name][-1])
                NewColumn = Table.Column(column.data, dtype='bool')
                table.replace_column(column.name, NewColumn)
        except:
            pass

    table.info()
    table.info('stats')
    print('Elapsed time(secs):', time.time() - t0)

    print('astropy.__version__', astropy.__version__)
    print('pandas.__version__', pd.__version__)
    print('numpy.__version__', np.__version__)

    print('Output file:', outfile)
    table.write(outfile)

    print('Elapsed time(secs):', time.time() - t0)

    return table


def sqlquery():

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
      /* prog = '177.A-3011' */
      /* AND is_stack = 'True' */
      AND chipno = 1
      /* LIMIT 100000 */
    """

    query12a="""
    SELECT
      is_stack, is_tile, naxis1, naxis2, version, chipno, dateobs, nightobs, prog
    FROM
      vstqc
    WHERE
      /* prog = '177.A-3011' */
      is_stack = 'True'
      AND chipno = 1
      /* LIMIT 100000 */
    """

    query12aa="""
    SELECT
      *
    FROM
      vstqc
    WHERE
      /* surveyname = 'ATLAS' */
      prog = '177.A-3011'
      /* is_stack = 'True'   */
      /* AND chipno = 1 */
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
    #print "Program ID: %s ;%s" % (progid, month)
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

    query = sql

    #% (month, progid)

    #query=query201208
    #query=query10a
    #query=query_all_tiles_thin
    #query=query_all_pawprints_thin
    #query=query11b

    # vst
    #query=query12a
    #query=sql

    return sql



if __name__ == '__main__':
    """

    """
    args = getargs()
    db = args.db
    print('db:', db)

    configfile = 'querydqc.cfg'
    config = getconfig(configfile=configfile, debug=True)

    host = config.get(db,'host')
    dbname = config.get(db,'dbname')
    user = config.get(db,'user')
    password = config.get(db,'password')
    port = config.get(db,'port')

    # Output file
    outfile = args.outfile
    print('Output file:', outfile)

    # Open and read the file as a single buffer
    sqlfile = args.sqlfile
    print('Reading:', sqlfile)

    fh = open(sqlfile, 'r')
    sql = fh.read()
    fh.close()

    print('sql:', sql)

    verbose=True
    if verbose:
        print('Python version:', sys.version)
        #print 'pyfits.__version__:', pyfits.__version__
        print('numpy.__version__:', np.__version__)
        print('psycopg version:', psycopg.__version__)

    main(outfile=outfile)
