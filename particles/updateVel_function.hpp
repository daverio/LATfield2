#ifndef LATFIELD2_UPDATEVEL_FUNCTION_HPP
#define LATFIELD2_UPDATEVEL_FUNCTION_HPP


/*
 move particles method for gevolution, made to work with any particles type
 */

using namespace LATfield2;


Real updateVel_gevolution(double dtau,
                          double lat_resolution,
                          part_simple * part,
                          double * ref_dist,
                          part_simple_info partInfo,
                          Field<Real> ** fields,
                          Site * sites,
                          int nfield,
                          double * params,
                          double * outputs,
                          int noutputs)
{


#define phi (fields[0])
#define xP (sites[0])
#define Bi (fields[1])
#define xB (sites[1])   
#define H_conformal params[0] 
#define rescaleB    params[1]
#define flag_init   params[2]

      
    double gradPhi[3]={0,0,0};
    double pgradB[3]={0,0,0};
    double v2;
    double dtauH_2 = dtau * H_conformal * 0.5l;
    
    gradPhi[0] = (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((*phi)(xP+0) - (*phi)(xP));
    gradPhi[1] = (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((*phi)(xP+1) - (*phi)(xP));
    gradPhi[2] = (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((*phi)(xP+2) - (*phi)(xP));
    gradPhi[0] += ref_dist[1] * (1.-ref_dist[2]) * ((*phi)(xP+1+0) - (*phi)(xP+1));
    gradPhi[1] += ref_dist[0] * (1.-ref_dist[2]) * ((*phi)(xP+1+0) - (*phi)(xP+0));
    gradPhi[2] += ref_dist[0] * (1.-ref_dist[1]) * ((*phi)(xP+2+0) - (*phi)(xP+0));
    gradPhi[0] += (1.-ref_dist[1]) * ref_dist[2] * ((*phi)(xP+2+0) - (*phi)(xP+2));
    gradPhi[1] += (1.-ref_dist[0]) * ref_dist[2] * ((*phi)(xP+2+1) - (*phi)(xP+2));
    gradPhi[2] += (1.-ref_dist[0]) * ref_dist[1] * ((*phi)(xP+2+1) - (*phi)(xP+1));
    gradPhi[0] += ref_dist[1] * ref_dist[2] * ((*phi)(xP+2+1+0) - (*phi)(xP+2+1));
    gradPhi[1] += ref_dist[0] * ref_dist[2] * ((*phi)(xP+2+1+0) - (*phi)(xP+2+0));
    gradPhi[2] += ref_dist[0] * ref_dist[1] * ((*phi)(xP+2+1+0) - (*phi)(xP+1+0));

    if(nfield==2 && Bi != NULL)
    {
        pgradB[0] = ((1.-ref_dist[2]) * ((*Bi)(xB+0,1) - (*Bi)(xB,1)) + ref_dist[2] * ((*Bi)(xB+2+0,1) - (*Bi)(xB+2,1))) * (*part).vel[1] / rescaleB;
        pgradB[0] += ((1.-ref_dist[1]) * ((*Bi)(xB+0,2) - (*Bi)(xB,2)) + ref_dist[1] * ((*Bi)(xB+1+0,2) - (*Bi)(xB+1,2))) * (*part).vel[2] / rescaleB;
        pgradB[0] += (1.-ref_dist[1]) * (1.-ref_dist[2]) * ((ref_dist[0]-1.) * (*Bi)(xB-0,0) + (1.-2.*ref_dist[0]) * (*Bi)(xB,0) + ref_dist[0] * (*Bi)(xB+0,0)) * (*part).vel[0] / rescaleB;
        pgradB[0] += ref_dist[1] * (1.-ref_dist[2]) * ((ref_dist[0]-1.) * (*Bi)(xB+1-0,0) + (1.-2.*ref_dist[0]) * (*Bi)(xB+1,0) + ref_dist[0] * (*Bi)(xB+1+0,0)) * (*part).vel[0] / rescaleB;
        pgradB[0] += (1.-ref_dist[1]) * ref_dist[2] * ((ref_dist[0]-1.) * (*Bi)(xB+2-0,0) + (1.-2.*ref_dist[0]) * (*Bi)(xB+2,0) + ref_dist[0] * (*Bi)(xB+2+0,0)) * (*part).vel[0] / rescaleB;
        pgradB[0] += ref_dist[1] * ref_dist[2] * ((ref_dist[0]-1.) * (*Bi)(xB+2+1-0,0) + (1.-2.*ref_dist[0]) * (*Bi)(xB+2+1,0) + ref_dist[0] * (*Bi)(xB+2+1+0,0)) * (*part).vel[0] / rescaleB;
        
        pgradB[1] = ((1.-ref_dist[0]) * ((*Bi)(xB+1,2) - (*Bi)(xB,2)) + ref_dist[0] * ((*Bi)(xB+1+0,2) - (*Bi)(xB+0,2))) * (*part).vel[2] / rescaleB;
        pgradB[1] += ((1.-ref_dist[2]) * ((*Bi)(xB+1,0) - (*Bi)(xB,0)) + ref_dist[2] * ((*Bi)(xB+1+2,0) - (*Bi)(xB+2,0))) * (*part).vel[0] / rescaleB;
        pgradB[1] += (1.-ref_dist[0]) * (1.-ref_dist[2]) * ((ref_dist[1]-1.) * (*Bi)(xB-1,1) + (1.-2.*ref_dist[1]) * (*Bi)(xB,1) + ref_dist[1] * (*Bi)(xB+1,1)) * (*part).vel[1] / rescaleB;
        pgradB[1] += ref_dist[0] * (1.-ref_dist[2]) * ((ref_dist[1]-1.) * (*Bi)(xB+0-1,1) + (1.-2.*ref_dist[1]) * (*Bi)(xB+0,1) + ref_dist[1] * (*Bi)(xB+0+1,1)) * (*part).vel[1] / rescaleB;
        pgradB[1] += (1.-ref_dist[0]) * ref_dist[2] * ((ref_dist[1]-1.) * (*Bi)(xB+2-1,1) + (1.-2.*ref_dist[1]) * (*Bi)(xB+2,1) + ref_dist[1] * (*Bi)(xB+2+1,1)) * (*part).vel[1] / rescaleB;
        pgradB[1] += ref_dist[0] * ref_dist[2] * ((ref_dist[1]-1.) * (*Bi)(xB+2+0-1,1) + (1.-2.*ref_dist[1]) * (*Bi)(xB+2+0,1) + ref_dist[1] * (*Bi)(xB+2+0+1,1)) * (*part).vel[1] / rescaleB;
        
        pgradB[2] = ((1.-ref_dist[1]) * ((*Bi)(xB+2,0) - (*Bi)(xB,0)) + ref_dist[1] * ((*Bi)(xB+2+1,0) - (*Bi)(xB+1,0))) * (*part).vel[0] / rescaleB;
        pgradB[2] += ((1.-ref_dist[0]) * ((*Bi)(xB+2,1) - (*Bi)(xB,1)) + ref_dist[0] * ((*Bi)(xB+2+0,1) - (*Bi)(xB+0,1))) * (*part).vel[1] / rescaleB;
        pgradB[2] += (1.-ref_dist[0]) * (1.-ref_dist[1]) * ((ref_dist[2]-1.) * (*Bi)(xB-2,2) + (1.-2.*ref_dist[2]) * (*Bi)(xB,2) + ref_dist[2] * (*Bi)(xB+2,2)) * (*part).vel[2] / rescaleB;
        pgradB[2] += ref_dist[0] * (1.-ref_dist[1]) * ((ref_dist[2]-1.) * (*Bi)(xB+0-2,2) + (1.-2.*ref_dist[2]) * (*Bi)(xB+0,2) + ref_dist[2] * (*Bi)(xB+0+2,2)) * (*part).vel[2] / rescaleB;
        pgradB[2] += (1.-ref_dist[0]) * ref_dist[1] * ((ref_dist[2]-1.) * (*Bi)(xB+1-2,2) + (1.-2.*ref_dist[2]) * (*Bi)(xB+1,2) + ref_dist[2] * (*Bi)(xB+2+1,2)) * (*part).vel[2] / rescaleB;
        pgradB[2] += ref_dist[0] * ref_dist[1] * ((ref_dist[2]-1.) * (*Bi)(xB+1+0-2,2) + (1.-2.*ref_dist[2]) * (*Bi)(xB+1+0,2) + ref_dist[2] * (*Bi)(xB+1+0+2,2)) * (*part).vel[2] / rescaleB;
    }
        
    gradPhi[0] += pgradB[0];
    gradPhi[1] += pgradB[1];
    gradPhi[2] += pgradB[2];
    
    gradPhi[0] /= lat_resolution;
    gradPhi[1] /= lat_resolution;
    gradPhi[2] /= lat_resolution;  
    
    //update vel
    v2 = 0.;
    if (flag_init == 0)
    {
        for(int i=0;i<3;i++)
        {
            (*part).vel[i] = ( ( (1.0 - dtauH_2)* (*part).vel[i] ) - (gradPhi[i]*dtau) ) / (1.0 + dtauH_2);
            v2 += (*part).vel[i] * (*part).vel[i];
        }
    }
    else
    {
        for(int i=0;i<3;i++)
        {
            (*part).vel[i] = -gradPhi[i] / (1.5 * H_conformal);
            //(*part).vel[i] -= 0.002;
            v2 += (*part).vel[i] * (*part).vel[i];
        }
    }
    
    return v2;
    
#undef phi
#undef xP
#undef Bi
#undef xB
#undef H_conformal
#undef rescaleB
#undef flag_init  
}



#endif
