import json
import numpy as np

import time
import xobjects as xo
import xline as xl
import xtrack as xt

# MPI:
from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


num_turns = int(1e6)
n_part = 300
N_saved_states = 11

# ion beam dimensions:
sigma_z = 0.063 # m
sigma_dp = 2e-4 # relative ion momentum spread
sigma_x = 1.047e-3 # m
sigma_y = 0.83e-3  # m

####################
# Choose a context #
####################

context = xo.ContextCpu()
#context = xo.ContextCpu(omp_num_threads=5)
#context = xo.ContextCupy()
#context = xo.ContextPyopencl('0.0')

buf = context.new_buffer()

##################
# Get a sequence #
##################

#sequence based on the LinearTransferMatrix element:
M = xt.LinearTransferMatrix(_buffer=buf, _context=context,
            Q_s=0.0131,
            beta_s=sigma_z/sigma_dp)

sequence = xl.Line(
    elements=[M],
    element_names=['LinearTransferMatrix']
)

#fname_sequence = '../../test_data/sps_w_spacecharge/line_no_spacecharge_and_particle.json'
#with open(fname_sequence, 'r') as fid:
#     input_data = json.load(fid)
#sequence = xl.Line.from_dict(input_data['line'])
#
#items, names = sequence.get_elements_of_type(xl.Cavity)
#for itm, name in zip(items, names):
#    print(name + ':\t' + str(itm))
#    
## Modify voltage:
#name = 'acta.31637'
#idx = sequence.find_element_ids(name)[0]
#cav = sequence.elements[idx]
#cav.voltage = 7e6 # V

# Gamma Factory Proof-of-Principle Experiment at the SPS:
# http://cds.cern.ch/record/2690736
# The parameters are from Table 3 (page 20)

# Ion mass
m_u = 931.49410242e6 # eV/c^2 -- atomic mass unit
A = 207.98 # Lead-208
Z = 82  # Number of protons in the ion (Lead)
Ne = 3 # Number of remaining electrons (Lithium-like)
m_e = 0.511e6 # eV/c^2 -- electron mass
m_p = 938.272088e6 # eV/c^2 -- proton mass
c = 299792458.0 # m/s

m_ion = A*m_u + Ne*m_e # eV/c^2

equiv_proton_momentum = 236e9 # eV/c = gamma_p*m_p*v

gamma_p = np.sqrt( 1 + (equiv_proton_momentum/m_p)**2 ) # equvalent gamma for protons in the ring

# R = pc/qB = > R*B = pc/q => p*c/(Z-Ne) = p_proton*c
# => p*c = (Z-Ne)*p_proton*c

p0c = equiv_proton_momentum*(Z-Ne) # eV/c
gamma = np.sqrt( 1 + (p0c/m_ion)**2 ) # ion relativistic factor
print("Ion gamma = %e" % gamma)

beta = np.sqrt(1-1/(gamma*gamma)) # ion beta

# laser-ion beam collision angle
theta_l = 2.6*np.pi/180 # rad
nx = 0; ny = -np.sin(theta_l); nz = -np.cos(theta_l)

# Ion excitation energy:
hw0 = 230.823 # eV
hc = 0.19732697e-6 # eV*m (ħc)
lambda_0 = 2*np.pi*hc/hw0 # m -- ion excitation wavelength

lambda_l = lambda_0*gamma*(1 + beta*np.cos(theta_l)) # m -- laser wavelength
# Shift laser wavelength for fast longitudinal cooling:
lambda_l = lambda_l*(1+sigma_dp) # m

laser_frequency = c/lambda_l # Hz
sigma_w = 2*np.pi*laser_frequency*sigma_dp
#sigma_w = 2*np.pi*laser_frequency*sigma_dp/2 # for fast longitudinal cooling

sigma_t = 1/sigma_w # sec -- Fourier-limited laser pulse
print('Laser pulse dration sigma_t = %.2f ps' % (sigma_t/1e-12))

print('Laser wavelength = %.2f nm' % (lambda_l/1e-9))

# Add Gamma Factory IP:
GF_IP = xt.IonLaserIP(_buffer=buf,
                      laser_direction_nx = 0,
                      laser_direction_ny = ny,
                      laser_direction_nz = nz,
                      laser_energy         = 5e-3, # J
                      laser_duration_sigma = sigma_t, # sec
                      laser_wavelength = lambda_l, # m
                      laser_waist_radius = 1.3e-3, # m
                      ion_excitation_energy = hw0, # eV
                      ion_excited_lifetime  = 76.6e-12, # sec
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
                         mass0 = m_ion, # eV/c^2
                         q0    = Z-Ne,
                         p0c   = p0c, # eV
                         x     = np.random.normal(scale=sigma_x, size=n_part),
                         px    = np.random.normal(scale=1e-5, size=n_part),
                         y     = np.random.normal(scale=sigma_y, size=n_part),
                         py    = np.random.normal(scale=1e-5, size=n_part),
                         zeta  = np.random.normal(scale=sigma_z, size=n_part),
                         delta = np.random.normal(scale=sigma_dp, size=n_part),
                         )

#########
# Track #
#########

print('Tracking...')
t0 = time.time()

turns_between_saves = int(num_turns/(N_saved_states-1))

saved_states = []
saved_turns  = []

turn = 0
saved_states.append(particles.copy())
saved_turns.append(turn)
for isave in range(N_saved_states-1):
    tracker.track(particles, num_turns=turns_between_saves)
    turn = turn + turns_between_saves
    saved_states.append(particles.copy())
    saved_turns.append(turn)
    if rank == 0:
        s_left = (num_turns-turn)*(time.time()-t0)/turn
        h_left = np.floor(s_left/3600)
        m_left = np.floor((s_left - h_left*3600)/60)
        s_left = s_left - h_left*3600 - m_left*60
        print(f"turn = {turn}, {h_left:.0f} h {m_left:.0f} min {s_left:.0f} sec left...")

#print(f"time = {time.time()-t0 : .1f} sec")

# MPI data management:
# from https://rabernat.github.io/research_computing/parallel-programming-with-mpi-for-python.html

zeta_all  = None
x_all     = None
y_all     = None
delta_all = None
p0c_all   = None
pID_all   = None
state_all = None

if rank == 0:
    zeta_all  = np.empty(n_part*size*N_saved_states, dtype='float64')
    x_all     = np.empty(n_part*size*N_saved_states, dtype='float64')
    y_all     = np.empty(n_part*size*N_saved_states, dtype='float64')
    delta_all = np.empty(n_part*size*N_saved_states, dtype='float64')
    p0c_all   = np.empty(n_part*size*N_saved_states, dtype='float64')
    pID_all   = np.empty(n_part*size*N_saved_states, dtype='int64')
    state_all = np.empty(n_part*size*N_saved_states, dtype='int64')

zeta   = np.array([part.zeta  for part in saved_states]).transpose()
x      = np.array([part.x     for part in saved_states]).transpose()
y      = np.array([part.y     for part in saved_states]).transpose()
delta  = np.array([part.delta for part in saved_states]).transpose()
p0c    = np.array([part.p0c   for part in saved_states]).transpose()
pID    = np.array([part.particle_id for part in saved_states]).transpose()
state  = np.array([part.state for part in saved_states]).transpose()

comm.Gather(zeta.flatten(),  zeta_all,  root=0)
comm.Gather(x.flatten(),     x_all,     root=0)
comm.Gather(y.flatten(),     y_all,     root=0)
comm.Gather(delta.flatten(), delta_all, root=0)
comm.Gather(p0c.flatten(),   p0c_all,   root=0)
comm.Gather(pID.flatten(),   pID_all,   root=0)
comm.Gather(state.flatten(), state_all, root=0)

#print(f'MPI: size={size}, rank={rank}: zeta = \n{zeta}')

if rank == 0:
    print('Saving data...')
    #print(f'zeta_all = \n{zeta_all}')
    #print(f'zeta to save = \n{zeta_all.reshape(-1,N_saved_states).transpose()}')
    
    with open(f'track.json', "w") as f:
        json.dump({
            'turns': saved_turns,
            'zeta':   zeta_all.reshape(-1,N_saved_states).transpose().tolist(),
            'x':         x_all.reshape(-1,N_saved_states).transpose().tolist(),
            'y':         y_all.reshape(-1,N_saved_states).transpose().tolist(),
            'delta': delta_all.reshape(-1,N_saved_states).transpose().tolist(),
            'p0c':     p0c_all.reshape(-1,N_saved_states).transpose().tolist(),
            'pID':     pID_all.reshape(-1,N_saved_states).transpose().tolist(),
    #        'parent_pID': parent_pID.tolist(),
            'state': state_all.reshape(-1,N_saved_states).transpose().tolist()
                  }, f, indent=1)

    print('Done.')
