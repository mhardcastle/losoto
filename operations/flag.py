#!/usr/bin/env python
# -*- coding: utf-8 -*-

# This operation for LoSoTo implement a flagging procedure

import logging
from operations_lib import *

logging.debug('Loading FLAG module.')

def run( step, parset, H ):

    import numpy as np
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    axesToMed = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    clipLevel = parset.getFloat('.'.join(["LoSoTo.Steps", step, "SigLevel"]), 0. )
    
    if len(axesToMed) < 1:
        logging.error("Please specify axes.")
        return 1
    if clipLevel == 0.:
        logging.error("Please specify significance level above which to clip.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("Flagging soltab: "+soltab._v_name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        # some checks
        if len(axesToMed) < 1:
            logging.error("Please specify axes.")
            return 1

        if clipLevel == 0.:
            logging.error("Please specify significance level above which to clip.")
            return 1

        for i, axis in enumerate(axesToMed[:]):
            if axis not in sf.getAxesNames():
                del axesToMed[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')

        before_count=0
        after_count=0
        total=0
        for vals, coord in sf.getValuesIter(returnAxes=axesToMed):

            total+=len(vals)
            before_count+=np.count_nonzero(np.isnan(vals))

            # clipping
            # first find the median and standard deviation
            valmedian = np.median(vals)
            clipvalue = valmedian * clipLevel
            np.putmask(vals, vals > clipvalue, np.nan)
            clipvalue = valmedian / clipLevel
            np.putmask(vals, vals < clipvalue, np.nan)
        
            after_count+=np.count_nonzero(np.isnan(vals))

            # writing back the solutions
            coord = removeKeys(coord, axesToMed)
            sw.setSelection(**coord)
            sw.setValues(vals)

        sw.addHistory('FLAG (over %s with %s sigma cut)' % (axesToMed, clipLevel))
        logging.info('Clip: %i points initially bad, %i after flagging (%f %%)' % (before_count,after_count,float(after_count)/total))
        
    return 0


