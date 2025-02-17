// copyright ############################### //
// This file is part of the Xtrack Package.  //
// Copyright (c) CERN, 2021.                 //
// ######################################### //

#ifndef XTRACK_LIMITRECT_H
#define XTRACK_LIMITRECT_H

/*gpufun*/
void LimitRect_track_local_particle(LimitRectData el, LocalParticle* part0){

    double const min_x = LimitRectData_get_min_x(el);
    double const max_x = LimitRectData_get_max_x(el);
    double const min_y = LimitRectData_get_min_y(el);
    double const max_y = LimitRectData_get_max_y(el);

    //start_per_particle_block (part0->part)

        double const x = LocalParticle_get_x(part);
        double const y = LocalParticle_get_y(part);

	int64_t const is_alive = (int64_t)(
                      (x >= min_x) &&
		      (x <= max_x) &&
		      (y >= min_y) &&
		      (y <= max_y) );

	// I assume that if I am in the function is because
    	if (!is_alive){
           LocalParticle_set_state(part, 0);
	}

    //end_per_particle_block

}

#endif
