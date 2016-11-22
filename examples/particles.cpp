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
            case 'i':
                io_size =  atoi(argv[++i]);
                break;
            case 'g':
                io_groupe_size = atoi(argv[++i]);
                break;
            case 'o':
                str_filename = argv[++i];
                break;
            case 's':
                npts = atoi(argv[++i]);
                break;
            case 'p':
                numparts = atoi(argv[++i]);
                break;
            case 'r':
                latresolution = atof(argv[++i]);
                break;
        }
    }

#ifndef EXTERNAL_IO
    parallel.initialize(n,m);
#else

    parallel.initialize(n,m,io_size,io_groupe_size);
    if(parallel.isIO()) ioserver.start();
    else
    {

#endif
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


        //cout<<"init done"<<endl;

        part_simple part;

        long index =0;
        /*
        for(int i=0;i<numparts;i++)
            for(int j=0;j<numparts;j++)
                for(int k=0;k<numparts;k++){
                    part.ID=index;
                    part.pos[0]= (Real)i * (Real)boxSize[0] / (Real)numparts;
                    part.pos[1]= (Real)j * (Real)boxSize[1] / (Real)numparts;
                    part.pos[2]= (Real)k * (Real)boxSize[2] / (Real)numparts;
                    part.vel[0]=1.0;
                    part.vel[1]=1.0;
                    part.vel[2]=1.0;
                    //part.mass=0.22;
                    parts.addParticle_global(part);
                    index++;
        }
        */
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

#ifdef EXTERNAL_IO
        while( !ioserver.openOstream());
        parts.saveHDF5_server_open("data/bench_part_server");
#endif


        timerRef = MPI_Wtime();
        parts.saveHDF5("bench_part",1);
        timerWrite = MPI_Wtime() - timerRef;


        timerRef = MPI_Wtime();
        parts.loadHDF5("bench_part",1);
        timerLoad = MPI_Wtime() - timerRef;



#ifdef EXTERNAL_IO
        timerRef = MPI_Wtime();
        parts.saveHDF5_server_write();
        timerWriteServer = MPI_Wtime() - timerRef;

        ioserver.closeOstream();
#endif


        /*
        // cout<<"write done"<<endl;

        projection_init(&phi);
        projection_init(&B);
        projection_init(&Tij);

        //cout<<"init proj done"<<endl;


        timerRef = MPI_Wtime();
        scalarProjectionCIC_project(&parts,&phi);
        timerProjScalar = MPI_Wtime() - timerRef;

        //cout<<"scalar proj done"<<endl;


        timerRef = MPI_Wtime();
        scalarProjectionCIC_comm(&phi);
        timerCommScalar = MPI_Wtime() - timerRef;

        //cout<<"scalar comm done"<<endl;

        timerRef = MPI_Wtime();
        vectorProjectionCICNGP_project(&parts,&B);
        timerProjVector = MPI_Wtime() - timerRef;

        //cout<<"vector proj done"<<endl;

        timerRef = MPI_Wtime();
        vectorProjectionCICNGP_comm(&B);
        timerCommVector = MPI_Wtime() - timerRef;

        //cout<<"vector comm done"<<endl;

        timerRef = MPI_Wtime();
        symtensorProjectionCICNGP_project(&parts,&Tij);
        timerProjTensor = MPI_Wtime() - timerRef;


        //cout<<"tensor proj done"<<endl;

        timerRef = MPI_Wtime();
        symtensorProjectionCICNGP_comm(&Tij);
        timerCommTensor = MPI_Wtime() - timerRef;

        //cout<<"tensor comm done"<<endl;

        timerRef = MPI_Wtime();
        parts.updateVel(&updateVel_simple,1.0);
        timerVel = MPI_Wtime() - timerRef;


        //cout<<"update vel done"<<endl;

        timerRef = MPI_Wtime();
        parts.moveParticles(&move_particles_simple,1.0/ (1.0* (double)npts));
        timerMove = MPI_Wtime() - timerRef;

        //cout<<"move done"<<endl;

        ofstream textfile;
        if(parallel.isRoot())
        {
            textfile.open(str_filename.c_str(),ios::out | ios::app);
            textfile<< n<<","<<m<< ","<< io_size<< ","<<io_groupe_size<< ","<< npts << ","<< numparts << ",";
            textfile<< timerWrite << ","<<  timerLoad << ","<<  timerWriteServer << ","<<  timerProjScalar << ","<<  timerCommScalar << ",";
            textfile<< timerProjVector << ","<<  timerCommVector << ","<< timerProjTensor << ","<< timerCommTensor << ","<< timerVel << ","<< timerMove <<endl;



            textfile.close();
        }
        */
#ifdef EXTERNAL_IO
        ioserver.stop();
    }
#endif
}
