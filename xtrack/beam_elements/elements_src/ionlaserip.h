#ifndef XTRACK_IONLASERIP_H
#define XTRACK_IONLASERIP_H

/*gpufun*/
void IonLaserIP_track_local_particle(IonLaserIPData el, LocalParticle* part0){

    //The algorithm is from https://anaconda.org/petrenko/psi_beam_vs_laser

    double const nx  = IonLaserIPData_get_laser_direction_nx(el);
    double const ny  = IonLaserIPData_get_laser_direction_ny(el);
    double const nz  = IonLaserIPData_get_laser_direction_nz(el);
    double const laser_x = IonLaserIPData_get_laser_x(el);
    double const laser_y = IonLaserIPData_get_laser_y(el);
    double const laser_z = IonLaserIPData_get_laser_z(el);
    double const w0 = IonLaserIPData_get_laser_waist_radius(el);

    double const laser_energy = IonLaserIPData_get_laser_energy(el);
    double const laser_waist_shift = IonLaserIPData_get_laser_waist_shift(el);

    double const laser_sigma_t = IonLaserIPData_get_laser_duration_sigma(el);
    double const laser_sigma_w = 1/laser_sigma_t; // rad/sec -- assuming Fourier-limited pulse

    double const laser_wavelength = IonLaserIPData_get_laser_wavelength(el); // Hz
    double const c = 299792458.0; // m/s

    double const ion_excited_lifetime = IonLaserIPData_get_ion_excited_lifetime(el); // sec
    double const ion_excitation_energy = IonLaserIPData_get_ion_excitation_energy(el); // eV
    double const eV = 1.602176634e-19; // J
    
    double const laser_Rayleigh_length = PI*w0*w0/laser_wavelength;
    //printf("\nlaser_Rayleigh_length=%e m\n",laser_Rayleigh_length); exit(1);

    double const p0c = LocalParticle_get_p0c(part0); // eV
    double const m0  = LocalParticle_get_mass0(part0); // eV/c^2
    //printf("m0=%e,p0c=%e,p0c/m0=%f\n",m0,p0c,p0c/m0);
    double const hbar = 1.054571817e-34; // J*sec

    // Maximum laser intensity (at the focal point)
    double const I0 = sqrt(2/PI)*(laser_energy/laser_sigma_t)/(PI*w0*w0); // W/m^2
    
    double state,delta,z,x,y,px,py,pc,gamma,beta,beta_x,beta_y,beta_z,vx,vy,vz,tcol;
    double r2, Z_to_laser_focus, I, OmegaRabi, w;

    //start_per_particle_block (part0->part)
            
        // Map_of_Excitation_vs_Intensity_and_Detuning
        double v0 = \
            IonLaserIPData_get_Map_of_Excitation_vs_Intensity_and_Detuning(el, 0);
        double v1 = \
            IonLaserIPData_get_Map_of_Excitation_vs_Intensity_and_Detuning(el, 1);
        double v2 = \
            IonLaserIPData_get_Map_of_Excitation_vs_Intensity_and_Detuning(el, 999);
        //printf("\n\nTest v0=%e, v1=%e, v2=%e\n\n",v0,v1,v2); exit(1);

    
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
    
        //printf("OmegaRabi*laser_sigma_t/(2*gamma) = %e\n",OmegaRabi*laser_sigma_t/(2*gamma)); exit(0);
    
        //Test:
        if (state > 0)
        {
            if (OmegaRabi*laser_sigma_t/(2*gamma) > 1.5)
            {
                LocalParticle_set_state(part, 2);
            }
            else {
                LocalParticle_set_state(part, 1);
            }
        }
        
        // ...
        // laser cooling demostration -- temporary !!!!
        
        if (delta > 0)
        {
            //LocalParticle_add_to_energy(part, -5e6); // eV
        }    

    //end_per_particle_block

}

#endif
