import pathlib
import json
import numpy as np

import xobjects as xo
import xline as xl
import xtrack as xt

fname_sequence = '../../test_data/lhc_no_bb/line_and_particle.json'

num_turns = int(2000)
n_part = 50

####################
# Choose a context #
####################

context = xo.ContextCpu()
#context = xo.ContextCupy()
#context = xo.ContextPyopencl('0.0')

buf = context.new_buffer()

##################
# Get a sequence #
##################

with open(fname_sequence, 'r') as fid:
     input_data = json.load(fid)
sequence = xl.Line.from_dict(input_data['line'])

# Add Gamma Factory IP
theta_l = 2.6*np.pi/180 # rad
nx = 0; ny = -np.sin(theta_l); nz = -np.cos(theta_l)

GF_IP = xt.IonLaserIP(_buffer=buf,
                      laser_direction_nx = 0,
                      laser_direction_ny = ny,
                      laser_direction_nz = nz,
                      laser_energy=5e-3, # J
                     )
sequence.append_element(GF_IP, 'GammaFactory_IP')

##################
# Build TrackJob #
##################
tracker = xt.Tracker(_context=context, _buffer=buf, sequence=sequence)

######################
# Get some particles #
######################
particles = xt.Particles(_context=context,
                         p0c=6500e9,
                         x=np.random.uniform(-1e-3, 1e-3, n_part),
                         px=np.random.uniform(-1e-5, 1e-5, n_part),
                         y=np.random.uniform(-2e-3, 2e-3, n_part),
                         py=np.random.uniform(-3e-5, 3e-5, n_part),
                         zeta=np.random.uniform(-1e-2, 1e-2, n_part),
                         delta=np.random.uniform(-1e-4, 1e-4, n_part),
                         )
#########
# Track #
#########

print('Tracking...')
tracker.track(particles, num_turns=num_turns, turn_by_turn_monitor=True)

print('Saving data...')
zeta  = tracker.record_last_track.zeta
delta = tracker.record_last_track.delta
p0c   = tracker.record_last_track.p0c
pID   = tracker.record_last_track.particle_id
parent_pID = tracker.record_last_track.parent_particle_id
state = tracker.record_last_track.state

#print(np.shape(delta))
#print(delta)

with open('track.json', "w") as f:
    json.dump({
        'zeta': zeta.tolist(),
        'delta': delta.tolist(),
        'p0c': p0c.tolist(),
        'pID': pID.tolist(),
        'parent_pID': parent_pID.tolist(),
        'state': state.tolist()
              }, f, indent=1)

print('Done.')