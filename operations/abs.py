#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Take absolute value. Needed if amplitudes are negative!

import logging
from operations_lib import *

logging.debug('Loading ABS module.')

def run( step, parset, H ):

    import numpy as np
    from h5parm import solFetcher, solWriter

    soltabs = getParSoltabs( step, parset, H )
    ants = getParAxis( step, parset, H, 'ant' )
    pols = getParAxis( step, parset, H, 'pol' )
    dirs = getParAxis( step, parset, H, 'dir' )

    axesToAbs = parset.getStringVector('.'.join(["LoSoTo.Steps", step, "Axes"]), [] )
    
    if len(axesToAbs) < 1:
        logging.error("Please specify axes.")
        return 1

    for soltab in openSoltabs( H, soltabs ):

        sf = solFetcher(soltab)
        sw = solWriter(soltab)

        logging.info("ABSing soltab: "+soltab._v_name)

        sf.setSelection(ant=ants, pol=pols, dir=dirs)

        # some checks

        for i, axis in enumerate(axesToAbs[:]):
            if axis not in sf.getAxesNames():
                del axesToAbs[i]
                logging.warning('Axis \"'+axis+'\" not found. Ignoring.')


        for vals, coord in sf.getValuesIter(returnAxes=axesToAbs):

            valsnew=np.abs(vals)

            # writing back the solutions
            coord = removeKeys(coord, axesToAbs)
            sw.setSelection(**coord)
            sw.setValues(valsnew)

        sw.addHistory('ABS (over %s)' % axesToAbs)

        
    return 0


