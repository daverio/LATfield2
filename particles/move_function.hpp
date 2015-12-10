#ifndef LATFIELD2_MOVE_FUNCTION_HPP
#define LATFIELD2_MOVE_FUNCTION_HPP




/**
 * \addtogroup prartClass
 * @{
 */

using namespace LATfield2;


void move_particles_simple(double dtau,
                               double lat_resolution,
                               part_simple * part,
                               double * ref_dist,
                               part_simple_info partInfo,
                               Field<Real> ** fields,
                               Site * sites,
                               int nfield,
                               double * params,
                               double * outputs,
                               int noutputs){


    //double a;
    //a = 1 + 23;
    for (int l=0;l<3;l++) (*part).pos[l] += dtau*(*part).vel[l];
   
}

/**@}*/

#endif
