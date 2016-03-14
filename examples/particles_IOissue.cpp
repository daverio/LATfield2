#include <stdlib.h>
#include "LATfield2.hpp"

using namespace LATfield2;

int main(int argc, char **argv)
{
    
    
    int n,m;
    int io_groupe_size,io_size;
    string str_filename;
    int npts = 64;
    int numparts = 64;
    Real  latresolution =0.1;
    
    
    for (int i=1 ; i < argc ; i++ ){
        if ( argv[i][0] != '-' )
            continue;
        switch(argv[i][1]) {
            case 'n':
                n = atoi(argv[++i]); //size of the dim 1 of the processor grid
                break;
            case 'm':
                m =  atoi(argv[++i]); //size of the dim 2 of the processor grid
                break;
        }
    }
    

    parallel.initialize(n,m);

        int dim=3;
        int halo=1;
        int khalo=1;
        
        
        Lattice lat_part(dim,npts,0);
        Lattice lat(dim,npts,halo);
        
        Field<Real> phi(lat,1);
        Field<Real> B(lat,3);
        Field<Real> Tij(lat,3,3,LATfield2::symmetric);
        
        Real boxSize[3];
        for(int i=0;i<3;i++)boxSize[i] = latresolution * lat_part.size(i);
        
        double timerRef;
        double timerWrite,timerLoad,timerWriteServer;
        
        double timerProjScalar,timerProjVector,timerProjTensor;
        double timerCommScalar,timerCommVector,timerCommTensor;
        
        double timerMove,timerVel;
        
        Site x(lat);
        
        for(x.first();x.test();x.next())
        {
            phi(x)=0.;
            for(int i=0;i<3;i++)B(x,i)=0.;
        }
        
        
        //cout<<"start part init done"<<endl;
        
        part_simple_info particles_global_info;
        part_simple_dataType particles_dataType;
        
        
        particles_global_info.mass=0.1;
        particles_global_info.relativistic=false;
        set_parts_typename(&particles_global_info,"part_simple");

        
        
        Particles<part_simple,part_simple_info,part_simple_dataType> parts;
        parts.initialize(particles_global_info,particles_dataType,&lat_part,boxSize);

	part_simple part;
        long index =0;
        int ratio = numparts/npts;
        //ratio;
        
        Site xp(lat_part);
        for(xp.first();xp.test();xp.next())
        {
            for(int i=0;i<ratio;i++)
                for(int j=0;j<ratio;j++)
                    for(int k=0;k<ratio;k++){
            
            part.ID=index;
            part.pos[0]= (Real)xp.coord(0) * (Real)boxSize[0] / (Real)npts;
            part.pos[1]= (Real)xp.coord(1) * (Real)boxSize[1] / (Real)npts;
            part.pos[2]= (Real)xp.coord(2) * (Real)boxSize[2] / (Real)npts;
            part.vel[0]=1.0;
            part.vel[1]=1.0;
            part.vel[2]=1.0;
            //part.mass=0.22;
            parts.addParticle_global(part);
            index++;
            }
        }

        cout<<"implementation done"<<endl;
        
        
        timerRef = MPI_Wtime();
        parts.saveHDF5("bench_part",2);
        timerWrite = MPI_Wtime() - timerRef;
        
        
        timerRef = MPI_Wtime();
        parts.loadHDF5("bench_part",2);
        timerLoad = MPI_Wtime() - timerRef;
        
}

