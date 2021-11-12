#ifndef XTRACK_IONLASERIP_H
#define XTRACK_IONLASERIP_H

/*gpufun*/
void IonLaserIP_track_local_particle(IonLaserIPData el, LocalParticle* part0){

    //The algorithm is partially from https://anaconda.org/petrenko/psi_beam_vs_laser

    double nx  = IonLaserIPData_get_laser_direction_nx(el);
    double ny  = IonLaserIPData_get_laser_direction_ny(el);
    double nz  = IonLaserIPData_get_laser_direction_nz(el);
    
    double laser_x = IonLaserIPData_get_laser_x(el);
    double laser_y = IonLaserIPData_get_laser_y(el);
    double laser_z = IonLaserIPData_get_laser_z(el);
    double w0 = IonLaserIPData_get_laser_waist_radius(el);

    double laser_energy = IonLaserIPData_get_laser_energy(el);
    double laser_waist_shift = IonLaserIPData_get_laser_waist_shift(el);

    double laser_sigma_t = IonLaserIPData_get_laser_duration_sigma(el);
    //double laser_sigma_w = 1/laser_sigma_t; // rad/sec -- assuming Fourier-limited pulse

    double laser_wavelength = IonLaserIPData_get_laser_wavelength(el); // Hz
    double c = 299792458.0; // m/s

    double ion_excited_lifetime = IonLaserIPData_get_ion_excited_lifetime(el); // sec
    double ion_excitation_energy = IonLaserIPData_get_ion_excitation_energy(el); // eV
    double eV = 1.602176634e-19; // J
    
    // constants for the Map_of_Excitation vs OmegaRabi and Detuning:
    int64_t N_OmegaRabiTau_values = IonLaserIPData_get_N_OmegaRabiTau_values(el);
    int64_t N_DeltaDetuningTau_values = IonLaserIPData_get_N_DeltaDetuningTau_values(el);
    double  OmegaRabiTau_max = IonLaserIPData_get_OmegaRabiTau_max(el);
    double  DeltaDetuningTau_max = IonLaserIPData_get_DeltaDetuningTau_max(el);
    double  dOmegaRabiTau = OmegaRabiTau_max/(N_OmegaRabiTau_values-1.0);
    double  dDeltaDetuningTau = DeltaDetuningTau_max/(N_DeltaDetuningTau_values-1.0);
    
    double laser_Rayleigh_length = PI*w0*w0/laser_wavelength;
    //printf("\nlaser_Rayleigh_length=%e m\n",laser_Rayleigh_length); exit(1);

    double p0c = LocalParticle_get_p0c(part0); // eV
    double m0  = LocalParticle_get_mass0(part0); // eV/c^2
    //printf("m0=%e,p0c=%e,p0c/m0=%f\n",m0,p0c,p0c/m0);
    double hbar = 1.054571817e-34; // J*sec
    
    
    //double gamma0 = sqrt(1.0 + p0c*p0c/(m0*m0));
    //double beta0  = sqrt(1.0 - 1.0/(gamma0*gamma0));
    double OmegaTransition = ion_excitation_energy*eV/hbar; // rad/sec
    //printf("2.0*PI*c/OmegaTransition = %e m, laser_wavelength=%e m\n",
    //        2.0*PI*c/OmegaTransition, laser_wavelength/((1.0+beta0)*gamma0)); exit(0);

    // Maximum laser intensity (at the focal point)
    double I0 = sqrt(2/PI)*(laser_energy/laser_sigma_t)/(PI*w0*w0); // W/m^2
    
    double state,delta,z,x,y,px,py,pc,gamma,beta,beta_x,beta_y,beta_z,vx,vy,vz,tcol;
    double r2, Z_to_laser_focus, I, OmegaRabi, OmegaRabiTau, DeltaDetuningTau, w;
    double laser_omega_ion_frame, cos_theta, excitation_probability, rnd;
    
    int64_t row, col, idx;

    //start_per_particle_block (part0->part)
    
        state = LocalParticle_get_state(part);
        delta = LocalParticle_get_delta(part);
        z     = LocalParticle_get_zeta(part);
        x     = LocalParticle_get_x(part);
        y     = LocalParticle_get_y(part);
        px    = LocalParticle_get_px(part);
        py    = LocalParticle_get_py(part);
    
        pc = p0c*(1.0+delta); // eV
        gamma = sqrt(1.0 + pc*pc/(m0*m0));
        beta  = sqrt(1.0 - 1.0/(gamma*gamma));
        beta_x  = px*p0c/m0/gamma;
        beta_y  = py*p0c/m0/gamma;
        beta_z  = sqrt(beta*beta - beta_x*beta_x -beta_y*beta_y);

        vx  = c*beta_x; // m/sec
        vy  = c*beta_y;
        vz  = c*beta_z;

        //printf("\n\n beta_x=%e, beta_y=%e, beta_z=%e \n\n",beta_x,beta_y,beta_z);
        //exit(1);
    
        // Collision of ion with the laser pulse:
        // The position of the laser beam center is rl=rl0+ct*n. We can find the moment
        // when a particle with a position r=r0+vt collides with the laser as the moment
        // when r−rl is perpendicular to n. Then (r−rl,n)=0, which yields the equation
        // (r0,n)+(v,n)t−(rl0,n)−ct(n,n)=0. Hence
        // tcol=(r0−rl0,n)/[c−(v,n)]
    
        tcol = ( (x-laser_x)*nx + (y-laser_y)*ny + (z-laser_z)*nz ) / (c - (vx*nx+vy*ny+vz*nz)); // sec
        // r^2 to the laser center = |r-rl| at the moment tcol:
        r2 = (\
                pow(x+vx*tcol - (laser_x+c*nx*tcol), 2) + \
                pow(y+vy*tcol - (laser_y+c*ny*tcol), 2) + \
                pow(z+vz*tcol - (laser_z+c*nz*tcol), 2) \
             ); // m

        Z_to_laser_focus = laser_waist_shift - tcol*c; // m
        
        // Laser beam size at the point of collision:
        w = w0*sqrt(1.0+pow(Z_to_laser_focus/laser_Rayleigh_length, 2));
        // Max. laser intensity experienced by the ion (in the ion rest frame):
        I = 4.0*gamma*gamma * I0*(w0/w)*(w0/w)*exp(-2.0*r2/(w*w)); // W/m^2
        OmegaRabi = \
            (hbar*c/(ion_excitation_energy*eV)) * \
            sqrt(I*2*PI/(ion_excitation_energy*eV*ion_excited_lifetime)); // rad/sec
    
        OmegaRabiTau = OmegaRabi*laser_sigma_t/(2.0*gamma); // in the ion rest frame
    
        // Detuning from the ion transition resonance in the ion rest frame:
        
        cos_theta = -(nx*vx + ny*vy + nz*vz)/(beta*c);
        laser_omega_ion_frame = (2.0*PI*c/laser_wavelength)*(1.0+beta*cos_theta)*gamma;
        DeltaDetuningTau = fabs(
            (OmegaTransition - laser_omega_ion_frame)*laser_sigma_t/(2.0*gamma)
        );
        //printf("DeltaDetuningTau = %e\n", DeltaDetuningTau); if (ii>15) {exit(0);}
            
        //Test:
        if (state > 0)
        {
            // Map_of_Excitation vs OmegaRabi and Detuning:
            //double v0 = IonLaserIPData_get_Map_of_Excitation(el, 0);
            //double v1 = IonLaserIPData_get_Map_of_Excitation(el, 1);
            //double v2 = IonLaserIPData_get_Map_of_Excitation(el, 999);
            //printf("\n\nTest v0=%e, v1=%e, v2=%e\n\n",v0,v1,v2); exit(1);
            if (DeltaDetuningTau < DeltaDetuningTau_max && OmegaRabiTau > OmegaRabiTau_max)
            {
                // In case of a very high laser field:
                LocalParticle_set_state(part, 2); // Excited particle
            }
            else if (DeltaDetuningTau < DeltaDetuningTau_max &&
                OmegaRabiTau > dOmegaRabiTau/10.0)
            {
                // N_OmegaRabiTau_values  N_DeltaDetuningTau_values
                //   OmegaRabiTau_max       DeltaDetuningTau_max
                 
                row = (int)floor(OmegaRabiTau/dOmegaRabiTau);
                col = (int)floor(DeltaDetuningTau/dDeltaDetuningTau);
                idx = row*N_DeltaDetuningTau_values + col;
                
                excitation_probability = IonLaserIPData_get_Map_of_Excitation(el, idx);
                //printf("%.3f, %.3f, %.3f\n",OmegaRabiTau,DeltaDetuningTau,
                //      excitation_probability);
                //if (ii>19998){exit(0);}
                
                rnd = (float)rand()/(float)(RAND_MAX);
                if ( rnd < excitation_probability )
                {
                    LocalParticle_set_state(part, 2); // Excited particle
                    
                    // photon recoil (from emitted photon!):
                    rnd = (float)rand()/(float)(RAND_MAX);
                    LocalParticle_add_to_energy(part,
                                                -ion_excitation_energy*rnd*2.0*gamma); // eV
                } else {
                    LocalParticle_set_state(part, 1); // Still particle
                }
            }
            else {
                LocalParticle_set_state(part, 1); // Still particle
            }
        }
        
        // ...
        // laser cooling demostration -- temporary !!!!
        
        //if (delta > 0)
        //{
        //    //LocalParticle_add_to_energy(part, -5e6); // eV
        //}    

    //end_per_particle_block
    
}

#endif
