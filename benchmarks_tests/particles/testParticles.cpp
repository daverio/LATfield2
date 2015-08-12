#include <stdlib.h>
#include "LATfield2.hpp"


using namespace LATfield2;



int main(int argc, char **argv)
{
	
	
	int n,m;

	
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
    int khalo=2;
	int npts[3]={64,128,64};    
    int numparts=16;

    Real boxSize[3]={1.,2.,1.};

    Real  latresolution;
    latresolution = get_lattice_resolution(npts,boxSize);
    
    
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
    particles_global_info.mass=0.21;
    particles_global_info.type_name="part_simple";
    particles_global_info.relativistic=true;
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
                part.vel[1]=0.2;
                part.vel[2]=0.3;
                //part.mass=0.22;
                parts.addParticle_global(part);
                index++;
            }
  
    
    /*
    parts.saveHDF5("test",3);
    
    MPI_Barrier(MPI_COMM_WORLD);
    parts_verif.loadHDF5("test",3);
    
    Site x(parts.lat());
    
    std::list<part_simple>::iterator it,it_verif;
    
    
    for(x.first();x.test();x.next())
    {
        if(parts.field_part()(x).size!=0)
        {
            for(it=parts.field_part()(x).parts.begin(), it_verif=parts_verif.field_part()(x).parts.begin();it != parts.field_part()(x).parts.end();++it,++it_verif)
            {
                if((*it).ID != (*it_verif).ID)cout<< "argarg: "<< (*it).ID << " : " << (*it_verif).ID<<endl;
            }
        }
    }
     */
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
        //cout<<i<<endl;
        parts.coutPart(0);
        parts.updateVel(&updateVel_gevolution,1,listField_updateVel,2,updateVel_gevolution_params,output,output_type,6);
        parts.moveParticles(&move_particles_gevolution,1,listField_move,0,&rescaleB,output,output_type,6);
    }
    
    
    
    
    
    has_maxi_mass<part_simple> hmm_part;
    has_maxi_mass<part_simple_info> hmm_info;
    
    
    
    //cout<<"mass part:" << hmm_part.gos() << endl;
    //cout<<"mass info:" << hmm_info.gos() << endl; 
    
    projection_init(&phi);
    size_t offset_mass_simple = hmm_part.gos();
    
    
    scalarProjectionCIC_project(&parts,&phi);
    scalarProjectionCIC_comm(&phi);
    
    vectorProjectionCICNGP_project(&parts,&B);
    //vectorProjectionCIC_comm(&phi);
}

