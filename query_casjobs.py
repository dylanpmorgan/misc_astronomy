__author__ = 'dpmorg'

import finding_cwdm.code.fileprep as fp
import casjobs
import pdb
import numpy as np
from astropy.table import Table, vstack
import time
import sys
import os.path

def main(query, casID=None, casPass=None, context='DR12'):
    """
    Submits the input query to the SDSS database and downloads the result
    locally.

    -- Sends query
    -- Monitors query until Finished
    -- Downloads the output file
    -- Forms the filenames from returned data and adds it to the input file.
    """

    # Try to connect with casjobs
    try:
        jobs = casjobs.CasJobs(casID, casPass)
    except:
        sys.exit("\n Connecting to casjobs failed, did you enter the ""
                 "correct casID and password? \n")

    counter = 0
    for ra, dec in zip(data['ra'], data['dec']):
        the_time =
        table_name =
        temp_filename =

        if os.path.isfile(output_path+file_name):
            counter += 1
            print 'Temp file already found -- Moving on: '+str(file_name)
            continue


        # Grab the text for the query.
        query = get_query(query_name(ra, dec, table_name))
        # Submit jobs to casjobs
        job_id = jobs.submit(query, context='DR12')

        counter = 0
        for ra, dec in zip(data['ra'], data['dec']):
            # Table name in casjobs
            table_name = 'extinct_'+str(counter)
            # Temporary file name, will be deleted.
            file_name = 'temp_extinct_'+"{:.5f}".format(ra)+'_'+"{:.5f}".format(dec)+'.fits'

            # test if file already exists, then continue
            if os.path.isfile(output_path+file_name):
                counter += 1
                print 'File already found -- Moving on: '+str(file_name)
                continue

            # Form query
            query = sdss_extinct_query(ra, dec, table_name)

            # Submit job to casjobs
            job_id = jobs.submit(query, context='DR12')

            # Monitor the status of the job.
            monitor_status = 0
            while monitor_status != 5:
                monitor_status = jobs.monitor(job_id)[0]
                print monitor_status
                if monitor_status == 4:
                    sys.exit('ERROR: Query failed.')
                #time.sleep(5)

            # Request and download the output from the query.
            jobs.request_and_get_output(table_name, 'FITS', output_path+file_name)

            # Read in the datafile
            t = Table.read(output_path+file_name)

            if len(t) == 1:
                # Save exinction information
                data[counter]['EXTINCTION_U'] = t['extinction_u'][0]
                data[counter]['EXTINCTION_G'] = t['extinction_g'][0]
                data[counter]['EXTINCTION_R'] = t['extinction_r'][0]
                data[counter]['EXTINCTION_I'] = t['extinction_i'][0]
                data[counter]['EXTINCTION_Z'] = t['extinction_z'][0]

            jobs.drop_table(table_name)

            print 'Finished: '+str(file_name)
            data.write('/Users/dpmorg/gdrive/research/cwdm_flares/new/data/'+'s82_cwdm_phot_ext.fits', overwrite=True)
            counter += 1

        pdb.set_trace()
        data.write('/Users/dpmorg/gdrive/research/cwdm_flares/new/data/'+'s82_cwdm_phot_ext.fits', overwrite=True)

def submit_query(query, casID=None, casPass=None):
    """
    Submits the input query to the SDSS database and downloads the result
    locally.

    -- Sends query
    -- Monitors query until Finished
    -- Downloads the output file
    -- Forms the filenames from returned data and adds it to the input file.
    """




def submit_individual_query(query, casID=None, casPass=None, input_file=None):
    '''
    Queries for single objects in SDSS spectroscopic database
    and returns all spectra for the given object.

    -- Sends query
    -- Monitors query until Finished
    -- Downloads the outputfile
    -- Forms the filenames from returned data and adds it to the input
    file.
    '''

    try:
        data = pyfits.getdata(input_file)
        ndata = len(data)
        print(" Input file found, %f lines." %(ndata)

        if len(data) > 1000:
            print("WARNING: Your file has more than 1000 entries,"
                  "when given an input file, this script runs queries for each"
                  "object individually. This query may take a long time and"
                  "make some folks mad. Try uploading the table into casjobs first"
                  "and then scripting your query to use that table as an input.")
    else:
        sys.exit("\n No input file, run query script using submit_query")

    try:
        ra_in = data['RA']
    else:
        sys.exit("No RA found in file")

    try:
        dec_in = data['DEC']
    else:
        sys.exit("No DEC found in file")

    output = pd.DataFrame()
    for ra, dec in zip(ra_in, dec_in):
        # Table name in casjobs
        the_time = str(time.time())
        table_name = 'query_individual_'+the_time
        # Temporary file name, will be deleted.
        file_name = 'temp_allspec_'+"{:.5f}".format(ra)+'_'+"{:.5f}".format(dec)+'.fits'

        # Form query
        query = sdss_allspec_query(ra, dec, table_name)

        # Submit job to casjobs
        job_id = jobs.submit(query, context='DR12')

        # Monitor the status of the job.
        monitor_status = 0
        while monitor_status != 5:
            monitor_status = jobs.monitor(job_id)[0]
            print monitor_status

        # Request and download the output from the query.
        jobs.request_and_get_output(table_name, 'FITS', output_path+file_name)

        # Read in the datafile
        if counter == 0:
            t = Table.read(output_path+file_name)
        elif counter > 0:
            t0 = Table.read(output_path+file_name)
            t = vstack([t, t0])

        jobs.drop_table(table_name)

        print 'Finished: '+str(file_name)
        counter += 1

    pdb.set_trace()
    t.write('/Users/dpmorg/gdrive/research/cwdm_flares/new/data/'+'s82_cwdm_allspec.fits', overwrite=True)




def s82_lc_query(raOut, decOut, queryname="casjobs", arcSearch=2.0):

    arcSearch = arcSearch/3600.

    query = """ SELECT ALL
    p.objid, p.ra AS ra, p.dec AS dec,
    p.skyVersion, p.run, p.rerun, p.camcol, p.field, p.flags,
    f.mjd_u, f.mjd_g, f.mjd_r, f.mjd_i, f.mjd_z,
    p.psfMag_u, p.psfMagErr_u,
    p.psfMag_g, p.psfMagErr_g,
    p.psfMag_r, p.psfMagErr_r,
    p.psfMag_i, p.psfMagErr_i,
    p.psfmag_z, p.psfMagErr_z,
    p.cx, p.cy,
    p.dered_u, p.dered_g, p.dered_r, p.dered_i, p.dered_z,
    p.extinction_u, p.extinction_g, p.extinction_r, p.extinction_i, p.extinction_z,
    p.flags_u, p.flags_g, p.flags_r, p.flags_i, p.flags_z into mydb.%s from PhotoObjAll AS p
    JOIN FIELD as f on p.skyVersion = f.skyVersion and p.run = f.run and p.camcol = f.camcol and p.field = f.field
    WHERE p.ra BETWEEN %s AND %s AND
          p.dec BETWEEN  %s AND %s
          AND ((flags_u & 0x10000000) != 0)
          AND ((flags_g & 0x10000000) != 0)
          --AND ((flags_r & 0x10000000) != 0)
          --AND ((flags_i & 0x10000000) != 0)
          --AND ((flags_z & 0x10000000) != 0)
          -- detected in BINNED1
          AND ((flags_u & 0x8100000c00a4) = 0)
          AND ((flags_g & 0x8100000c00a4) = 0)
          --AND ((flags_r & 0x8100000c00a4) = 0)
          --AND ((flags_i & 0x8100000c00a4) = 0)
          --AND ((flags_z & 0x8100000c00a4) = 0)
          -- not EDGE, NOPROFILE, PEAKCENTER, NOTCHECKED, PSF_FLUX_INTERP,
          -- SATURATED, or BAD_COUNTS_ERROR
          AND (((flags_u & 0x400000000000) = 0) or (psfmagerr_u <= 0.5))
          AND (((flags_g & 0x400000000000) = 0) or (psfmagerr_g <= 0.5))
          --AND (((flags_r & 0x400000000000) = 0) or (psfmagerr_r <= 0.5))
          --AND (((flags_i & 0x400000000000) = 0) or (psfmagerr_i <= 0.5))
          --AND (((flags_z & 0x400000000000) = 0) or (psfmagerr_z <= 0.5))
          -- not DEBLEND_NOPEAK or small PSF error
          -- (substitute psfmagerr in other band as appropriate)
          AND (((flags_u & 0x100000000000) = 0) or (flags_u & 0x1000) = 0)
          AND (((flags_g & 0x100000000000) = 0) or (flags_g & 0x1000) = 0)
          --AND (((flags_r & 0x100000000000) = 0) or (flags_r & 0x1000) = 0)
          --AND (((flags_i & 0x100000000000) = 0) or (flags_i & 0x1000) = 0)
          --AND (((flags_z & 0x100000000000) = 0) or (flags_z & 0x1000) = 0)
          -- not INTERP_CENTER or not COSMIC_RAY
    """ %(queryname, raOut-arcSearch, raOut+arcSearch, decOut-arcSearch, decOut+arcSearch)

    return query

def sdss_allspec_query(raOut, decOut, queryname="casjobs", arcSearch=2.0):

    arcSearch = arcSearch/3600.

    query = """ SELECT ALL
    s.specobjid, s.objid,
    s.ra as ra, s.dec as dec,
    s.plate, s.mjd, s.fiberid,
    sp.snMedian, sp.snMedian_u, sp.snMedian_g, sp.snMedian_r, sp.snMedian_i, sp.snMedian_z,
    s.run, sp.instrument, sp.run2d INTO mydb.%s
    FROM SpecPhotoAll as s
        LEFT JOIN SpecObjAll as sp on s.specobjid = sp.specobjid
    WHERE
      sp.snMedian > 2.0 AND
      s.ra BETWEEN %s AND %s AND
      s.dec BETWEEN  %s AND %s
    ORDER BY s.ra
    """ %(queryname, raOut-arcSearch, raOut+arcSearch, decOut-arcSearch, decOut+arcSearch)

    return query

def sdss_extinct_query(raOut, decOut, queryname="casjobs", arcSearch=2.0):

    arcSearch = arcSearch/3600.

    query = """ SELECT ALL
    p.objid, p.ra as ra, p.dec as dec,
    p.extinction_u, p.extinction_g, p.extinction_r, p.extinction_i, extinction_z
    INTO mydb.%s
    FROM PhotoPrimary as p
    WHERE
      p.ra BETWEEN %s AND %s AND
      p.dec BETWEEN  %s AND %s
    ORDER BY p.ra
    """ %(queryname, raOut-arcSearch, raOut+arcSearch, decOut-arcSearch, decOut+arcSearch)

    return query

def sdss_allspec(filename):
    '''
    Queries for single objects in SDSS spectroscopic database
    and returns all spectra for the given object.

    -- Sends query
    -- Monitors query until Finished
    -- Downloads the outputfile
    -- Forms the filenames from returned data and adds it to the input
       file.
    '''
    output_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/data/temp/'
    data = Table.read(filename)

    casID = '191359142'
    casPass = 'parker'

    jobs = casjobs.CasJobs(casID, casPass)

    counter = 0
    t = Table()
    for ra, dec in zip(data['RA'], data['DEC']):
        # Table name in casjobs
        table_name = 'allspec_'+str(counter)
        # Temporary file name, will be deleted.
        file_name = 'temp_allspec_'+"{:.5f}".format(ra)+'_'+"{:.5f}".format(dec)+'.fits'

        # Form query
        query = sdss_allspec_query(ra, dec, table_name)

        # Submit job to casjobs
        job_id = jobs.submit(query, context='DR12')

        # Monitor the status of the job.
        monitor_status = 0
        while monitor_status != 5:
            monitor_status = jobs.monitor(job_id)[0]
            print monitor_status

        # Request and download the output from the query.
        jobs.request_and_get_output(table_name, 'FITS', output_path+file_name)

        # Read in the datafile
        if counter == 0:
            t = Table.read(output_path+file_name)
        elif counter > 0:
            t0 = Table.read(output_path+file_name)
            t = vstack([t, t0])

        jobs.drop_table(table_name)

        print 'Finished: '+str(file_name)
        counter += 1

    pdb.set_trace()
    t.write('/Users/dpmorg/gdrive/research/cwdm_flares/new/data/'+'s82_cwdm_allspec.fits', overwrite=True)


def s82_lc(filename):
    output_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/data/s82/fits/'
    data = Table.read(filename)

    casID = '191359142'
    casPass = 'parker'

    jobs = casjobs.CasJobs(casID, casPass)

    s82_file_lis = []
    counter = 0
    for ra, dec in zip(data['RA'], data['DEC']):
        table_name = 'cwdms82_'+str(counter)
        file_name = '{0:03d}'.format(counter)+'_specs82cwdm_'+"{:.5f}".format(ra)+'_'+"{:.5f}".format(dec)+'.fits'

        query = s82_lc_query(ra, dec, table_name)

        pdb.set_trace()

        job_id = jobs.submit(query, context='Stripe82')

        while monitor_status != 5:
            monitor_status = jobs.monitor(job_id)[0]
            print monitor_status

        jobs.request_and_get_output(table_name, 'FITS', output_path+file_name)

        jobs.drop_table(table_name)

        s82_file_lis.append(str(file_name))
        print 'Downloaded '+str(file_name)
        counter += 1

    pdb.set_trace()
    data['S82_FILE'] = s82_file_lis
    data.write(output_path+file_name, overwrite=True)

def sdss_extinct(filename):
    '''
    Queries for single objects in SDSS spectroscopic database
    and returns all spectra for the given object.

    -- Sends query
    -- Monitors query until Finished
    -- Downloads the outputfile
    -- Forms the filenames from returned data and adds it to the input
       file.
    '''
    output_path = '/Users/dpmorg/gdrive/research/cwdm_flares/new/data/temp/'
    data = Table.read(filename)

    casID = '191359142'
    casPass = 'parker'

    jobs = casjobs.CasJobs(casID, casPass)

    # Add extinction values
    if not 'EXTINCTION_U' in data.colnames: data['EXTINCTION_U'] = 0.
    if not 'EXTINCTION_G' in data.colnames: data['EXTINCTION_G'] = 0.
    if not 'EXTINCTION_R' in data.colnames: data['EXTINCTION_R'] = 0.
    if not 'EXTINCTION_I' in data.colnames: data['EXTINCTION_I'] = 0.
    if not 'EXTINCTION_Z' in data.colnames: data['EXTINCTION_Z'] = 0.

    counter = 0
    for ra, dec in zip(data['ra'], data['dec']):
        # Table name in casjobs
        table_name = 'extinct_'+str(counter)
        # Temporary file name, will be deleted.
        file_name = 'temp_extinct_'+"{:.5f}".format(ra)+'_'+"{:.5f}".format(dec)+'.fits'

        # test if file already exists, then continue
        if os.path.isfile(output_path+file_name):
            counter += 1
            print 'File already found -- Moving on: '+str(file_name)
            continue

        # Form query
        query = sdss_extinct_query(ra, dec, table_name)

        # Submit job to casjobs
        job_id = jobs.submit(query, context='DR12')

        # Monitor the status of the job.
        monitor_status = 0
        while monitor_status != 5:
            monitor_status = jobs.monitor(job_id)[0]
            print monitor_status
            if monitor_status == 4:
                sys.exit('ERROR: Query failed.')
            #time.sleep(5)

        # Request and download the output from the query.
        jobs.request_and_get_output(table_name, 'FITS', output_path+file_name)

        # Read in the datafile
        t = Table.read(output_path+file_name)

        if len(t) == 1:
            # Save exinction information
            data[counter]['EXTINCTION_U'] = t['extinction_u'][0]
            data[counter]['EXTINCTION_G'] = t['extinction_g'][0]
            data[counter]['EXTINCTION_R'] = t['extinction_r'][0]
            data[counter]['EXTINCTION_I'] = t['extinction_i'][0]
            data[counter]['EXTINCTION_Z'] = t['extinction_z'][0]

        jobs.drop_table(table_name)

        print 'Finished: '+str(file_name)
        data.write('/Users/dpmorg/gdrive/research/cwdm_flares/new/data/'+'s82_cwdm_phot_ext.fits', overwrite=True)
        counter += 1

    pdb.set_trace()
    data.write('/Users/dpmorg/gdrive/research/cwdm_flares/new/data/'+'s82_cwdm_phot_ext.fits', overwrite=True)
