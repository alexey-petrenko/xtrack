import sys
import pathlib
import numpy as np

import pymask as pm
import xtrack as xt
import xpart as xp
import xobjects as xo

from cpymad.madx import Madx

test_data_folder = pathlib.Path(
        __file__).parent.joinpath('../test_data').absolute()

def test_twiss():

    path = test_data_folder.joinpath('hllhc14_input_mad/')

    mad = Madx(command_log="mad_final.log")
    mad.call(str(path.joinpath("final_seq.madx")))
    mad.use(sequence="lhcb1")
    mad.twiss()
    mad.readtable(file=str(path.joinpath("final_errors.tfs")),
                  table="errtab")
    mad.seterr(table="errtab")
    mad.set(format=".15g")
    twmad = mad.twiss(rmatrix=True, chrom=True)

    line = xt.Line.from_madx_sequence(
            mad.sequence['lhcb1'], apply_madx_errors=True)
    part_ref = xp.Particles(mass0=xp.PROTON_MASS_EV, q0=1,
                            gamma0=mad.sequence.lhcb1.beam.gamma)

    for context in xo.context.get_test_contexts():
        print(f"Test {context.__class__}")
        tracker = xt.Tracker(_context=context, line=line)

        twxt = tracker.twiss(particle_ref=part_ref)

        for name in ['mb.b19r5.b1', 'mb.b19r1.b1',
                     'ip1', 'ip2', 'ip5', 'ip8',
                     'mbxf.4l1', 'mbxf.4l5']:

            imad = list(twmad['name']).index(name+':1')
            ixt = list(twxt['name']).index(name)

            assert np.isclose(twxt['betx'][ixt], twmad['betx'][imad],
                              atol=0, rtol=3e-4)
            assert np.isclose(twxt['bety'][ixt], twmad['bety'][imad],
                              atol=0, rtol=3e-4)
            assert np.isclose(twxt['dx'][ixt], twmad['dx'][imad], atol=1e-2)
            assert np.isclose(twxt['dy'][ixt], twmad['dy'][imad], atol=1e-2)
            assert np.isclose(twxt['dpx'][ixt], twmad['dpx'][imad], atol=3e-4)
            assert np.isclose(twxt['dpy'][ixt], twmad['dpy'][imad], atol=3e-4)

            assert np.isclose(twxt['s'][ixt], twmad['s'][imad], atol=5e-6)
            assert np.isclose(twxt['x'][ixt], twmad['x'][imad], atol=5e-6)
            assert np.isclose(twxt['y'][ixt], twmad['y'][imad], atol=5e-6)
            assert np.isclose(twxt['px'][ixt], twmad['px'][imad], atol=1e-7)
            assert np.isclose(twxt['py'][ixt], twmad['py'][imad], atol=1e-7)
