# copyright ############################### #
# This file is part of the Xtrack Package.  #
# Copyright (c) CERN, 2021.                 #
# ######################################### #

import time                          #!skip-doc
import numpy as np

import xtrack as xt
import xpart as xp

######################################
# Create a dummy collimation process #
######################################

class DummyInteractionProcess:
    '''
    I kill some particles. I kick some others by an given angle
    and I generate some secondaries with the opposite angles.
    '''
    def __init__(self, fraction_lost, fraction_secondary, length, kick_x):

        self.fraction_lost = fraction_lost
        self.fraction_secondary = fraction_secondary
        self.kick_x = kick_x
        self.length = length

        self.drift = xt.Drift(length=self.length)


    def interact(self, particles):

        self.drift.track(particles)

        n_part = particles._num_active_particles

        # Kill some particles
        mask_kill = np.random.uniform(size=n_part) < self.fraction_lost
        particles.state[:n_part][mask_kill] = -1 # special flag`


        # Generate some more particles
        mask_secondary = np.random.uniform(size=n_part) < self.fraction_secondary
        n_products = np.sum(mask_secondary)
        if n_products>0:
            products = {
                's': particles.s[:n_part][mask_secondary],
                'x': particles.x[:n_part][mask_secondary] + self.length * (
                    particles.px[:n_part][mask_secondary] + self.kick_x),
                'px': particles.px[:n_part][mask_secondary] + self.kick_x,
                'y': particles.y[:n_part][mask_secondary],
                'py': particles.py[:n_part][mask_secondary],
                'zeta': particles.zeta[:n_part][mask_secondary],
                'delta': particles.delta[:n_part][mask_secondary],

                'mass_ratio': particles.x[:n_part][mask_secondary] *0 + 1.,
                'charge_ratio': particles.x[:n_part][mask_secondary] *0 + 1.,

                'parent_particle_id': particles.particle_id[:n_part][mask_secondary],
                'at_element': particles.at_element[:n_part][mask_secondary],
                'at_turn': particles.at_turn[:n_part][mask_secondary],
                }
        else:
            products = None

        return products

############################################################
# Create a beam interaction from the process defined above #
############################################################

# kick_x is large so that the particles are lost on the first aperture element
interaction_process=DummyInteractionProcess(length=1., kick_x=4e-2,
                                            fraction_lost=0.0,
                                            fraction_secondary=0.2)
collimator = xt.BeamInteraction(length=interaction_process.length,
                                      interaction_process=interaction_process)

########################################################
# Create a line including the collimator defined above #
########################################################

line = xt.Line(elements=[
    xt.Drift(length=2.),
    xt.LimitEllipse(a=2e-2, b=2e-2),
    collimator,
    xt.LimitEllipse(a=2e-2, b=2e-2),
    xt.Drift(length=10.),
    ])

#################
# Build tracker #
#################

tracker = xt.Tracker(line=line, global_xy_limit=1e3)

##########################
# Build particles object #
##########################

# We prepare empty slots to store the product particles that will be
# generated during the tracking.
particles = xp.Particles(
        _capacity=200000,
        x=np.zeros(100000))

#########
# Track #
#########

t1 = time.time()                                          #!skip-doc
tracker.track(particles)
t2 = time.time()                                          #!skip-doc

                                                          #!skip-doc
print(f'{t2-t1=:.2f}')                                    #!skip-doc

############################
# Loss location refinement #
############################

part_before = particles.copy()

loss_loc_refinement = xt.LossLocationRefinement(tracker,
                                            n_theta = 360,
                                            r_max = 0.5, # m
                                            dr = 50e-6,
                                            ds = 0.05,
                                            save_refine_trackers=True)

loss_loc_refinement.refine_loss_location(particles)

assert np.all(part_before.s == particles.s)