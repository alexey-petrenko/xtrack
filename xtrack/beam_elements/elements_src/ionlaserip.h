#ifndef XTRACK_IONLASERIP_H
#define XTRACK_IONLASERIP_H

/*gpufun*/
void IonLaserIP_track_local_particle(IonLaserIPData el, LocalParticle* part0){
    //first test with code form SRotation:

    //start_per_particle_block (part0->part)
    //The algorithm is from https://anaconda.org/petrenko/psi_beam_vs_laser
        
        double const nx = IonLaserIPData_get_laser_direction_nx(el);
        double const ny = IonLaserIPData_get_laser_direction_ny(el);
        double const nz = IonLaserIPData_get_laser_direction_nz(el);
    
        double const p0c = LocalParticle_get_p0c(part);

        double delta = LocalParticle_get_delta(part);
        double x     = LocalParticle_get_x(part);
        double y     = LocalParticle_get_y(part);
        double px    = LocalParticle_get_px(part);
        double py    = LocalParticle_get_py(part);
    
        // Collision of ion with the laser pulse
        
        // The position of the laser beam center is rl=a+c(t−tl)n. We can find the moment
        // when a particle with a position r=r0+vt collides with the laser as the moment
        // when r−rl is perpendicular to n. Then (r−rl,n)=0, which yields the equation
        // (r0,n)+(v,n)t−(a,n)−c(t−tl)(n,n)=0. Hence
        // tcol=[(r0−a,n)+ctl]/[c−(v,n)]
    
        // ...
        // laser cooling demostration -- temporary !!!!
        
        if (delta > 0)
        {
            LocalParticle_add_to_energy(part, -1e6); // eV
        }
    

    //end_per_particle_block

}

#endif
