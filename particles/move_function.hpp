#ifndef LATFIELD2_MOVE_FUNCTION_HPP
#define LATFIELD2_MOVE_FUNCTION_HPP


/*
 move particles method for gevolution, made to work with any particles type
 */

using namespace LATfield2;


void move_particles_gevolution(double dtau,
                               double lat_resolution,
                               part_simple * part,
                               double * ref_dist,
                               part_simple_info partInfo,
                               Field<Real> ** fields,
                               Site * sites,
                               int nfield,
                               double * params){

    
    if(nfield == 1)
    {
#define Bi (fields[0])
#define xB (sites[0])
#define rescaleB (params[0])
        
        Real b[3];
        
        b[0] = (*Bi)(xB, 0) * (1.-ref_dist[1]) * (1.-ref_dist[2]);
        b[1] = (*Bi)(xB, 1) * (1.-ref_dist[0]) * (1.-ref_dist[2]);
        b[2] = (*Bi)(xB, 2) * (1.-ref_dist[0]) * (1.-ref_dist[1]);
        b[1] += (*Bi)(xB+0, 1) * ref_dist[0] * (1.-ref_dist[2]);
        b[2] += (*Bi)(xB+0, 2) * ref_dist[0] * (1.-ref_dist[1]);
        b[0] += (*Bi)(xB+1, 0) * ref_dist[1] * (1.-ref_dist[2]);
        b[2] += (*Bi)(xB+1, 2) * (1.-ref_dist[0]) * ref_dist[1];
        b[0] += (*Bi)(xB+2, 0) * (1.-ref_dist[1]) * ref_dist[2];
        b[1] += (*Bi)(xB+2, 1) * (1.-ref_dist[0]) * ref_dist[2];
        b[1] += (*Bi)(xB+2+0, 1) * ref_dist[0] * ref_dist[2];
        b[0] += (*Bi)(xB+2+1, 0) * ref_dist[1] * ref_dist[2];
        b[2] += (*Bi)(xB+1+0, 2) * ref_dist[0] * ref_dist[1];

        for (int l=0;l<3;l++) (*part).pos[l] += dtau*((*part).vel[l] + b[l] / rescaleB);
#undef Bi
#undef xB
#undef rescaleB
        
    }
    else
    {
        for (int l=0;l<3;l++) (*part).pos[l] += dtau*(*part).vel[l];
    }
    
    
}



#endif
