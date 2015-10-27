#include <stdlib.h>
#include "LATfield2.hpp"


using namespace LATfield2;



int main(int argc, char **argv)
{
	
	
	int n,m;
    int io_groupe_size,io_size;
	
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
		}
	}
	
	
    parallel.initialize(n,m,io_size,io_groupe_size);
    
    if(parallel.isIO()) IO_Server.start();
    else
    {
        int dim=3;
        int halo=2;
        int khalo=2;
        int npts[3]={64 ,64,64};
        int numparts=16;
        
        Real boxSize[3]={1.,1.,1.};
        
        Real  latresolution = get_lattice_resolution(npts,boxSize);
        
        
        Lattice lat_part(dim,npts,0);
        Lattice lat(dim,npts,halo);
        
        Field<Real> phi(lat,1);
        Field<Real> B(lat,3);
        
        Site x(lat);
        
        for(x.first();x.test();x.next())
        {
            phi(x)=0.;
            for(int i=0;i<3;i++)B(x,i)=0.;
        }
        
        part_simple_info particles_global_info;
        particles_global_info.mass=0.1;
        particles_global_info.relativistic=false;
        set_parts_typename(&particles_global_info,"part_simple");
        
        
        
        part_simple_dataType particles_dataType;
        
        Particles<part_simple,part_simple_info,part_simple_dataType> parts,parts_verif;
        
        parts.initialize(particles_global_info,particles_dataType,&lat_part,boxSize);
        parts_verif.initialize(particles_global_info,particles_dataType,&lat_part,boxSize);
        
 
        part_simple part;
        long index =0;
        for(int i=0;i<numparts;i++)
            for(int j=0;j<numparts;j++)
                for(int k=0;k<numparts;k++){
                    part.ID=index;
                    part.pos[0]= (Real)i * (Real)boxSize[0] / (Real)numparts;
                    part.pos[1]= (Real)j * (Real)boxSize[1] / (Real)numparts;
                    part.pos[2]= (Real)k * (Real)boxSize[2] / (Real)numparts;
                    part.vel[0]=0.1;
                    part.vel[1]=0.1;
                    part.vel[2]=0.1;
                    //part.mass=0.22;
                    parts.addParticle_global(part);
                    index++;
                }
            }
        }
    
        
        parts.saveHDF5("test",2);
        
         
         parts_verif.loadHDF5("testserver",2);
         
         Site xpart(parts.lattice());
         
         std::list<part_simple>::iterator it,it_verif;
         
        
         for(xpart.first();xpart.test();xpart.next())
         {
             if(parts.field()(xpart).size!=0)
             {
                 for(it=parts.field()(xpart).parts.begin(), it_verif=parts_verif.field()(xpart).parts.begin();it != parts.field()(xpart).parts.end();++it,++it_verif)
                 {
                     if((*it).ID != (*it_verif).ID)cout<< "argarg: "<< (*it).ID << " : " << (*it_verif).ID<<endl;
                 }
             }
         }
    
        
        IO_Server.openOstream();
        parts.saveHDF5_server_open("testserver");
        parts.saveHDF5_server_write();
        
        IO_Server.closeOstream();
        
        
        //sleep(1);
        
        //parts.saveHDF5_server_open("testserver");
        //parts.saveHDF5_server_write();
        
        //IO_Server.closeOstream();
        
        
        
        Field<Real> * listField_move[1];
        //listField_move = new Field<Real>*[1];
        
        listField_move[0] = &B;
        
        
        Field<Real> * listField_updateVel[2];
        //listField_move = new Field<Real>*[2];
        
        listField_updateVel[0] = &phi;
        listField_updateVel[1] = &B;
        
        
        
        
        double H_conformal=1;
        int flag_init=1;
        Real rescaleB = 1;
        
        double updateVel_gevolution_params[3];
        updateVel_gevolution_params[0]=H_conformal;
        updateVel_gevolution_params[1]=rescaleB;
        updateVel_gevolution_params[2]=flag_init;
        
        double output[6];
        int output_type[6];
        
        
        output_type[0]=MIN;
        output_type[1]=SUM;
        output_type[2]=MAX;
        output_type[3]=MAX_LOCAL;
        output_type[4]=MIN_LOCAL;
        output_type[5]=SUM_LOCAL;
        
        
         for(int i=0;i<1;i++)
         {
             parts.updateVel(&updateVel_gevolution,1,listField_updateVel,2,updateVel_gevolution_params,output,output_type,6);
             parts.moveParticles(&move_particles_gevolution,1,listField_move,0,&rescaleB,output,output_type,6);
         }
        
        
        /*
         
         
         has_maxi_mass<part_simple> hmm_part;
         has_maxi_mass<part_simple_info> hmm_info;
         
         
         
         //cout<<"mass part:" << hmm_part.gos() << endl;
         //cout<<"mass info:" << hmm_info.gos() << endl;
         
         projection_init(&phi);
         size_t offset_mass_simple = hmm_part.gos();
         
         
         scalarProjectionCIC_project(&parts,&phi);
         scalarProjectionCIC_comm(&phi);
         
         vectorProjectionCICNGP_project(&parts,&B);
         vectorProjectionCICNGP_comm(&B);
         
         B.saveHDF5("test.h5");
         
        projection_init(&phi);
        scalarProjectionCIC_project(&parts,&phi);
        scalarProjectionCIC_comm(&phi);
        
        projection_init(&B);
        vectorProjectionCICNGP_project(&parts,&B);
        vectorProjectionCICNGP_comm(&B);
        
        B.saveHDF5("test.h5");
        
         */

        IO_Server.stop();
    }
	
	
    
}

