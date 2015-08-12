#ifndef LATFIELD2_PARTICLES_HPP
#define LATFIELD2_PARTICLES_HPP

#include "LATfield2_particle_description.hpp"
#include "particles_tools.hpp"
#include "move_function.hpp"
#include "updateVel_function.hpp"

#ifdef HDF5
#include "LATfield2_particlesIO.h"
#endif


CREATE_MEMBER_DETECTOR(mass)
CREATE_MEMBER_DETECTOR(ID)
CREATE_MEMBER_DETECTOR(vel)
CREATE_MEMBER_DETECTOR(pos)
CREATE_MEMBER_DETECTOR(type_name)

CREATE_MEMBER_DETECTOR_MAXI(mass)

#define GLOBAL_MASS     0
#define INDIVIDUAL_MASS 1
#define NO_MASS         2


#define SUM             1
#define MIN             2
#define MAX             4
#define SUM_LOCAL       8
#define MIN_LOCAL      16
#define MAX_LOCAL      32

using namespace LATfield2;

template <typename part, typename part_info, typename part_dataType>
class Particles;
#include "projections.hpp"

/*! \struct partList
 \breif structure which contains the particles lists, template for the particles field.
 */
template <typename  part>
struct  partList{
    int size;
    part parti;
    std::list<part>  parts;
    std::list<part>  partsTemp;
    partList() : size(0), parts(), partsTemp(){}
};

template <typename part, typename part_info, typename part_dataType>
class Particles
{

public:
  Particles(){;};
  ~Particles();

  void initialize(part_info part_global_info,
		  part_dataType part_datatype,
		  Lattice * lat_part,
		  Real boxSize[3]);
 
    
  void getPartProcess(part pcl,int * ranks);
  void getPartNewProcess(part pcl,int * ranks);
  void getPartCoord(part pcl,int * coord);
  void getPartCoordLocal(part pcl,int * coord);
  bool addParticle_global(part newPart);
    
    
    
    Real updateVel(Real (*updateVel_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,int,double*,double*,int),
                   double dtau,
                   Field<Real> ** fields,
                   int nfields,
                   double * params=NULL,
                   double * output=NULL,
                   int * reduce_type=NULL,
                   int noutput=0);
                   //Field<Real> * phi, Field<Real> * Bi, double rescaleB, double H_conformal, double dtau, int flag_init);
    
    void moveParticles( void (*move_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,int,double*,double*,int),
                       double dtau,
                       Field<Real> ** fields=NULL,
                       int nfields=0,
                       double * params=NULL,
                       double * output=NULL,
                       int * reduce_type=NULL,
                       int noutput=0);

  void saveHDF5(string filename_base, int fileNumber);
  void loadHDF5(string filename_base, int fileNumber);

    void coutPart(long ID);
    
    Lattice & lattice(){return lat_part_;};
    Field<partList<part> > & field(){return field_part_;};
    Real res(){return lat_resolution_;};
    part_info * parts_info();
    
    int mass_type(){return mass_type_;};
    size_t mass_offset(){return mass_offset_;}; 

private:

  part_info part_global_info_;
  part_dataType part_datatype_;

  Lattice lat_part_;
  Real  lat_resolution_;
  Real boxSize_[3];

    

  Field<partList<part> > field_part_;


  int mass_type_;
  size_t mass_offset_;

  long numParticles_;

};

template <typename part, typename part_info, typename part_dataType>
part_info * Particles<part,part_info,part_dataType>::parts_info()
{
    return &part_global_info_;
}

template <typename part, typename part_info, typename part_dataType>
Particles<part,part_info,part_dataType>::~Particles()
{
    field_part_.dealloc();
}
template <typename part, typename part_info, typename part_dataType> 
void Particles<part,part_info,part_dataType>::initialize(part_info part_global_info,
							 part_dataType part_datatype,
							 Lattice * lat_part,
							 Real  boxSize[3])
{
    
  part_global_info_ = part_global_info;
  COUT << "Initialization of the particles: "<< part_global_info_.type_name <<endl;

  part_datatype_=part_datatype;
  numParticles_ = 0;

  lat_part_.initialize(lat_part->dim(),lat_part->size(),0);

  field_part_.initialize(lat_part_);
  field_part_.alloc();

  //lat_resolution_ = lat_resolution;
  lat_resolution_ = get_lattice_resolution(lat_part->size(),boxSize);
  for(int i = 0;i < 3;i++)boxSize_[i]=boxSize[i];

  //cout<< "init done"<<endl;


  //particles structure verification

  if(has_ID<part>::value == 0){
    COUT<< "particles have to have the variable \"ID\" define as a long"<<endl;
    exit(-1);

  }
  if(has_pos<part>::value == 0){
    COUT<< "particles have to have the variable \"pos\" define as a Real[3]"<<endl;
    exit(-2);
  }
  if(has_vel<part>::value == 0){
    COUT<< "particles have to have the variable \"vel\" define as a Real[3]"<<endl;
    exit(-3);

  }


  //does the mass is global or individual
    
    has_maxi_mass<part> part_has_mass;
    has_maxi_mass<part_info> info_has_mass;

    
    
    if(part_has_mass.gos() != -1)
    {
        mass_offset_ = part_has_mass.gos();
        COUT<< "particles have individual mass" << endl;
        mass_type_= INDIVIDUAL_MASS;
    }
    else if(info_has_mass.gos() != -1)
    {
        mass_offset_ = info_has_mass.gos();
        COUT << "all particles have their mass set to: " << *(double*)((char*)&part_global_info_ + mass_offset_)<<endl;
        mass_type_= GLOBAL_MASS;
    }
    else
    {
        COUT<< "Particles have to have a mass!!! In part or in part_info."<<endl;
        COUT<< "no mass detected..."<<endl;
        mass_type_= GLOBAL_MASS;
    }
   
    
}


template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::getPartProcess(part pcl,int * ranks)
{
    int coord[3];
    getPartCoord(pcl,coord);
    ranks[0] = lat_part_.getRankDim0(coord[2]);
    ranks[1] = lat_part_.getRankDim1(coord[1]);
}
template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::getPartNewProcess(part pcl,int * ranks)
{
    int coord[3];
    for(int i=0;i<3;i++)coord[i] = (int)(floor(pcl.pos[i]/lat_resolution_)) %lat_part_.size(i); 
    
    if(pcl.pos[2] >=boxSize_[2]) ranks[0] = parallel.grid_size()[0];
    else if (pcl.pos[2] < 0.) ranks[0] = -1;
    else ranks[0] = lat_part_.getRankDim0(coord[2]);
    
    if(pcl.pos[1] >=boxSize_[1]) ranks[1] = parallel.grid_size()[1];
    else if (pcl.pos[1] < 0.) ranks[1] = -1;
    else ranks[1] = lat_part_.getRankDim1(coord[1]);
}


template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::getPartCoord(part pcl,int * coord)
{
    for(int i=0;i<3;i++)coord[i] = (int)(floor(pcl.pos[i]/lat_resolution_)) %lat_part_.size(i); 
}



template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::getPartCoordLocal(part pcl,int * coord)
{
    for(int i=0;i<3;i++)coord[i] = (int)(floor(pcl.pos[i]/lat_resolution_)) %lat_part_.size(i);
    coord[2]-=lat_part_.coordSkip()[0];
    coord[1]-=lat_part_.coordSkip()[1];
}

template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::coutPart(long ID)
{
    Site x(lat_part_);
    typename std::list<part>::iterator it;
    
    for(x.first();x.test();x.next())
    {
        if((field_part_)(x).size!=0)
        {
            for(it=(field_part_)(x).parts.begin(); it != (field_part_)(x).parts.end();++it)
            {
                if((*it).ID==ID)cout<< "Parallel ranks: ("<<parallel.grid_rank()[1]<<","<<parallel.grid_rank()[0]<<") ; "<< part_global_info_.type_name<<": "<<*it<<endl;
            }
        }
        
    }
    
}


template <typename part, typename part_info, typename part_dataType>
bool Particles<part,part_info,part_dataType>::addParticle_global(part newPart)
{
  Site x(lat_part_);
  int coord[3];

  this->getPartCoord(newPart,coord);
    
  if(x.setCoord(coord))
    {
      field_part_(x).size += 1;
      field_part_(x).parts.push_back(newPart);
      numParticles_ +=1;
      return true;
    }
  else
    {
      return false;
    }
}

template <typename part, typename part_info, typename part_dataType>
Real Particles<part,part_info,part_dataType>::updateVel(Real (*updateVel_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,int,double*,double*,int),
               double dtau,
               Field<Real> ** fields,
               int nfields,
               double * params,
               double * output,
               int * reduce_type,
               int noutput)
{
    if(nfields<1)
    {
        COUT<< "LATfield2::Particles::updateVel: updateVel needs at least 1 input field" <<endl;
        exit(-121);
    }
    
    Site xPart(lat_part_);
    
    Site  sites[nfields];
    for(int i=0;i<nfields;i++)sites[i].initialize(fields[i]->lattice());
    
    
    typename std::list<part>::iterator it;
    double frac[3];
    double x0;
    Real maxvel = 0.;
    Real v2;
    
    //cout<<"arg"<<endl;

    
    double * output_temp;
    output_temp =new double[noutput];
    
    if(noutput>0)for(int i=0;i<noutput;i++)
    {
        //COUT<<reduce_type[i]<<endl;
        if(reduce_type[i] & (SUM | SUM_LOCAL))
        {
            output[i]=0;
            
            //COUT<< "sum" <<endl;
        }
        else if(reduce_type[i] & (MIN | MIN_LOCAL))
        {
            output[i]=9223372036854775807;
            //COUT<<"min"<<endl;
        }
        else if(reduce_type[i] & (MAX | MAX_LOCAL))
        {
            output[i]=-9223372036854775807;
            //COUT<<"max"<<endl;
        }
    }
   
    for(int i=0;i<nfields;i++) sites[i].first();
    for(xPart.first() ; xPart.test(); xPart.next())
    {
        if(field_part_(xPart).size!=0)
        {
            for (it=(field_part_)(xPart).parts.begin(); it != (field_part_)(xPart).parts.end(); ++it)
            {
                for (int l=0; l<3; l++)
                    frac[l] = modf( (*it).pos[l] / lat_resolution_, &x0);
                
                
                v2 = updateVel_funct(dtau,
                           lat_resolution_,
                           &(*it),
                           frac,
                           part_global_info_,
                           fields,
                           sites,
                           nfields,
                           params,
                           output_temp,
                           noutput);
                 
                if(v2>maxvel)maxvel=v2;
                
                if(noutput>0)for(int i=0;i<noutput;i++)
                {
                    //COUT<<reduce_type[i]<<endl;
                    if(reduce_type[i] & (SUM | SUM_LOCAL))
                    {
                        output[i]+=output_temp[i];
                    }
                    else if(reduce_type[i] & (MIN | MIN_LOCAL))
                    {
                        if(output[i]>output_temp[i])output[i]=output_temp[i];
                    }
                    else if(reduce_type[i] & (MAX | MAX_LOCAL))
                    {
                        if(output[i]<output_temp[i])output[i]=output_temp[i];
                    }
                }

                
                
            }
        }
        for(int i=0;i<nfields;i++) sites[i].next();
    }

    if(noutput>0)for(int i=0;i<noutput;i++)
    {
        //COUT<<reduce_type[i]<<endl;
        if(reduce_type[i] & SUM)
        {
            parallel.sum(output[i]);
        }
        else if(reduce_type[i] & MIN)
        {
            parallel.min(output[i]);
        }
        else if(reduce_type[i] & MAX)
        {
            parallel.max(output[i]);
        }
    }

    delete[] output_temp;
    
    
    return sqrt(maxvel);


}


template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::moveParticles( void (*move_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,int,double*,double*,int),
                                                            double dtau,
                                                            Field<Real> ** fields,
                                                            int nfields,
                                                            double * params,
                                                            double * output,
                                                            int * reduce_type,
                                                            int noutput)
{
#ifdef DEBUG_MOVE
    cout<<parallel.rank()<<"; move start"<<endl;
#endif
    
    parallel.barrier();
    
    LATfield2::Site x(lat_part_);
    LATfield2::Site xNew(lat_part_);
    
    
    //LATfield2::Site xB;
    
    LATfield2::Site * sites;
    
    typename std::list<part>::iterator it,itTemp;
    //Real b[3];
    double frac[3];
    double x0;
    int localCoord[3];
    int newLocalCoord[3];
    
    bool flag[4];
    
        
    part **sendBuffer;
    part **recBuffer;
    
    
    sendBuffer = new part*[6];
    recBuffer = new part*[6];
    
    part * pTemp;
    
    long p;
    long bufferSize[6];
    long bufferSizeRec[6];
        
    std::list<part>  part_moveProc[8];
    
    
    if(nfields!=0)
    {
        sites = new LATfield2::Site[nfields];
        for(int i = 0;i<nfields;i++)
        {
            sites[i].initialize(fields[i]->lattice());
            sites[i].first();
        }
    }
    
    
    double * output_temp;
    output_temp =new double[noutput];
    
    if(noutput>0)for(int i=0;i<noutput;i++)
    {
        //COUT<<reduce_type[i]<<endl;
        if(reduce_type[i] & (SUM | SUM_LOCAL))
        {
            output[i]=0;
            
            //COUT<< "sum" <<endl;
        }
        else if(reduce_type[i] & (MIN | MIN_LOCAL))
        {
            output[i]=9223372036854775807;
            //COUT<<"min"<<endl;
        }
        else if(reduce_type[i] & (MAX | MAX_LOCAL))
        {
            output[i]=-9223372036854775807;
            //COUT<<"max"<<endl;
        }
    }

    
    
    for(x.first();x.test();x.next())
    {
        if(field_part_(x).size!=0)
        {
            for(int i=0;i<3;i++)localCoord[i] = x.coordLocal(i);
            
            for(it=field_part_(x).parts.begin(); it != field_part_(x).parts.end();)
            {
                itTemp = it;
                ++it;
                
#ifdef DEBUG_CONTROLSPEED
                double sizeTest[3]={boxSize_[0], boxSize_[1]/parallel.grid_size()[1], boxSize_[2]/parallel.grid_size()[0]};
                for(int l=0;l<3;l++)
                {
                    if(dtau*(*itTemp).vel[l] > sizeTest[l])
                    {
                        cout<<"velocity too big, abort request"<<endl;
                        cout<<"dtau * vel: "<< dtau*(*itTemp).vel[0]*lat_resolution_<<" , "<<dtau*(*itTemp).vel[1]*lat_resolution_<<" , "<<dtau*(*itTemp).vel[2]*lat_resolution_<< " , dtau: "<<dtau<<endl;
                        parallel.abortForce();
                    }
                }
#endif
                for (int l=0; l<3; l++)
                    frac[l] = modf( (*itTemp).pos[l] / lat_resolution_, &x0);
                
                move_funct(dtau,
                           lat_resolution_,
                           &(*itTemp),
                           frac,
                           part_global_info_,
                           fields,
                           sites,
                           nfields,
                           params,
                           output_temp,
                           noutput);
                
                
                if(noutput>0)for(int i=0;i<noutput;i++)
                {
                    //COUT<<reduce_type[i]<<endl;
                    if(reduce_type[i] & (SUM | SUM_LOCAL))
                    {
                        output[i]+=output_temp[i];
                    }
                    else if(reduce_type[i] & (MIN | MIN_LOCAL))
                    {
                        if(output[i]>output_temp[i])output[i]=output_temp[i];
                    }
                    else if(reduce_type[i] & (MAX | MAX_LOCAL))
                    {
                        if(output[i]<output_temp[i])output[i]=output_temp[i];
                    }
                }

            
                
                int partRanks[2];
                int thisRanks[2];
                getPartNewProcess((*itTemp),partRanks); 
                thisRanks[0] = parallel.grid_rank()[0];
                thisRanks[1] = parallel.grid_rank()[1];
                
                for(int i=0;i<3;i++)
                {
                    if((*itTemp).pos[i]<0)itTemp->pos[i] +=  boxSize_[i];
                    if((*itTemp).pos[i]>=boxSize_[i])itTemp->pos[i]-=boxSize_[i];
                }
                
               
                if(partRanks[0]==thisRanks[0] && partRanks[1]==thisRanks[1])
                {
                    int r[3];
                    getPartCoord(*itTemp, r);
                    if(!xNew.setCoord(r))
                    {
                        cout<<"arg"<< *itTemp <<" ; "<< partRanks[0]<< " , " << thisRanks[0] <<endl;
                    }
                    
                    getPartCoordLocal(*itTemp, newLocalCoord);
                    if(localCoord[0]!=newLocalCoord[0] || localCoord[1]!=newLocalCoord[1] || localCoord[2]!=newLocalCoord[2] )
                    {
                        xNew.setCoordLocal(newLocalCoord);
                        field_part_(xNew).partsTemp.splice(field_part_(xNew).partsTemp.end(),field_part_(x).parts,itTemp);
                        field_part_(xNew).size += 1;
                        field_part_(x).size -= 1;
                    }
                    
                }
                else if(partRanks[1]==thisRanks[1]-1)
                {
                    if(partRanks[0]==thisRanks[0])
                    {
                        part_moveProc[2].splice(part_moveProc[2].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]-1)
                    { 
                        part_moveProc[0].splice(part_moveProc[0].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]+1)
                    { 
                        part_moveProc[1].splice(part_moveProc[1].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else
                    {
                        cout<< "particle : "<<(*itTemp).ID<<" have move to far away (more than 1 proc)."<<endl;
                    }
                }
                else if(partRanks[1]==thisRanks[1]+1)
                {
                    if(partRanks[0]==thisRanks[0])
                    {
                        part_moveProc[5].splice(part_moveProc[5].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]-1)
                    {
                        part_moveProc[3].splice(part_moveProc[3].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]+1)
                    {
                        part_moveProc[4].splice(part_moveProc[4].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else
                    {
                        cout<< "particle : "<<(*itTemp).ID<<" have move to far away (more than 1 proc)."<<endl;
                    }
                }
                else if(partRanks[1]==thisRanks[1])
                {
                    if(partRanks[0]==thisRanks[0]-1)
                    {
                        part_moveProc[6].splice(part_moveProc[6].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]+1)
                    {
                        part_moveProc[7].splice(part_moveProc[7].end(),field_part_(x).parts,itTemp);
                        field_part_(x).size -= 1;
                    }
                    else
                    {
                        cout<< "particle : "<<(*itTemp).ID<<" have move to far away (more than 1 proc)."<<endl;
                    }
                }
                else
                {
                    cout<< "particle : "<<(*itTemp).ID<<" have move to far away (more than 1 proc)."<<endl;
                }
                
            }
        }
        
        if(nfields!=0) for(int i=0;i<nfields;i++) sites[i].next();
        
    }
    for(x.first();x.test();x.next())if((field_part_)(x).size!=0)(field_part_)(x).parts.splice((field_part_)(x).parts.end(), (field_part_)(x).partsTemp );
    
    if(noutput>0)for(int i=0;i<noutput;i++)
    {
        //COUT<<reduce_type[i]<<endl;
        if(reduce_type[i] & SUM)
        {
            parallel.sum(output[i]);
        }
        else if(reduce_type[i] & MIN)
        {
            parallel.min(output[i]);
        }
        else if(reduce_type[i] & MAX)
        {
            parallel.max(output[i]);
        }
    }
    
    delete[] output_temp;
    

    //remove the number of part...
    for(int i=0;i<8;i++)
    {
        numParticles_ -=  part_moveProc[i].size();
    }
    
    
   
    ////////////////////////////////////////////////////////////////
    //First send Y direction
    
    
    //pack data
    for(int i=0;i<6;i++)
    {
        bufferSize[i]=part_moveProc[i].size();
        if( bufferSize[i]!=0 )
        {
            sendBuffer[i] = new part[bufferSize[i]];
            for(it=part_moveProc[i].begin(),p=0; it != part_moveProc[i].end();++it,p++)sendBuffer[i][p] = (*it);
            part_moveProc[i].clear();
        }
    
    }
    
    if(parallel.grid_rank()[1]%2==0)
	{
		////////////////////////////////
		///send/rec to/from higher rank
		////////////////////////////////
		
        if(parallel.grid_rank()[1]!=parallel.grid_size()[1]-1)/// si pas le dernier alors envoie au +1
	    {
            //send
            parallel.send_dim1( &bufferSize[3], 3, parallel.grid_rank()[1]+1);
            for(int i=3;i<6;i++)
            {   
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[1]+1);
                }
            }
            
            //recieve
            parallel.receive_dim1( &bufferSizeRec[3], 3, parallel.grid_rank()[1]+1);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[1]+1);
                    
                }
            }
	    }
        
        //////////////////////////////
        ///send/rec to/from lower rank
        //////////////////////////////
		
        if(parallel.grid_rank()[1] != 0)     /// si pas le premier alors envoie au -1
	    {
            //send
            parallel.send_dim1( bufferSize, 3, parallel.grid_rank()[1]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[1]-1);
                }
            }
            //recieve
            parallel.receive_dim1( bufferSizeRec, 3, parallel.grid_rank()[1]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[1]-1);
                    
                }
            }
            
	    }
        else if(parallel.grid_size()[1]%2==0)  /// si pair et = 0 alors envoi au dernier
	    {
            //send
            parallel.send_dim1( bufferSize, 3, parallel.grid_size()[1]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , parallel.grid_size()[1]-1);
                }
            }
            //recieve
            parallel.receive_dim1( bufferSizeRec, 3, parallel.grid_size()[1]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , parallel.grid_size()[1]-1);
                    
                }
            }
	    }
        
        
	}
    else
    {
        //////////////////////////////
        ///rec/send from/to lower rank
        //////////////////////////////
        
        //tous recois du -1 puis envois au -1
        
        //recieve
        parallel.receive_dim1( bufferSizeRec, 3, parallel.grid_rank()[1]-1);
        for(int i=0;i<3;i++)
        {
            if(bufferSizeRec[i]!=0)
            {
                recBuffer[i]=new part[bufferSizeRec[i]];
                parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[1]-1);
            }
        }
        //send
        parallel.send_dim1( bufferSize, 3, parallel.grid_rank()[1]-1);
        for(int i=0;i<3;i++)
        {
            if(bufferSize[i]!=0)
            {
                parallel.send_dim1( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[1]-1);
            }
        }
        
        
        //////////////////////////////
        //rec/send from/to higher rank
        /////////////////////////////
        
        if(parallel.grid_rank()[1]!=parallel.grid_size()[1]-1)//si pas dernier alors recoi du +1
        {
            
            //recieve
            parallel.receive_dim1( &bufferSizeRec[3], 3, parallel.grid_rank()[1]+1);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[1]+1);
                }
            }
            
            //send
            parallel.send_dim1( &bufferSize[3], 3, parallel.grid_rank()[1]+1);
            for(int i=3;i<6;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[1]+1);
                }
            }
            
            
        }
        else // si dernier alors pair et recoi du premier
        {
            
            //recieve
            parallel.receive_dim1( &bufferSizeRec[3], 3, 0);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , 0);
                    
                }
            }
            
            //send
            parallel.send_dim1( &bufferSize[3], 3,0);
            for(int i=3;i<6;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , 0);
                }
            }
        }
    }
	
	//if unpair :
	/////////////
	
    if(parallel.grid_size()[1]%2!=0)
    {
        if(parallel.grid_rank()[1]==0)
        {
            //send
            parallel.send_dim1( bufferSize, 3, parallel.grid_size()[1]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , parallel.grid_size()[1]-1);
                }
            }
            //recieve
            parallel.receive_dim1( bufferSizeRec, 3, parallel.grid_size()[1]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , parallel.grid_size()[1]-1);
                }
            }
        }
        if(parallel.grid_rank()[1]==parallel.grid_size()[1]-1)
        {
            
            //recieve
            parallel.receive_dim1( &bufferSizeRec[3], 3, 0);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim1( recBuffer[i], bufferSizeRec[i] , 0);
                }
            }
            
            //send
            parallel.send_dim1( &bufferSize[3], 3,0);
            for(int i=3;i<6;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim1( sendBuffer[i], bufferSize[i] , 0);
                }
            }
        }
    }
    
    //unpack local list: rec[2] and rec[5]
    
    //add partnum
    numParticles_ += bufferSizeRec[2] + bufferSizeRec[5];
    for(int i=0;i<bufferSizeRec[2];i++)
    {
        this->getPartCoordLocal(recBuffer[2][i],newLocalCoord);
        x.setCoordLocal(newLocalCoord);
        
        field_part_(x).size += 1;
        field_part_(x).parts.push_back(recBuffer[2][i]);
#ifdef DEBUG_MOVE
        int verif;
        verif=addParticle_global(recBuffer[2][i]);
        
        if(verif != 1)
        {
            cout<<parallel.rank()<<"; MOVEBUF2 partID "<< recBuffer[2][i].ID<<" is not in the correct proc. " << recBuffer[2][i]<<endl;
        }
#endif
    }
    
    
    for(int i=0;i<bufferSizeRec[5];i++)
    {
        this->getPartCoordLocal(recBuffer[5][i],newLocalCoord);
        
        x.setCoordLocal(newLocalCoord);
        
        field_part_(x).size += 1;
        field_part_(x).parts.push_back(recBuffer[5][i]);        
#ifdef DEBUG_MOVE
        int verif;
        verif=addParticle_global(recBuffer[5][i]);
        
        if(verif != 1)
        {
            cout<<parallel.rank()<<"; MOVEBUF5 partID "<< recBuffer[5][i].ID<<" is not in the correct proc. "<<recBuffer[5][i] <<endl;
        }
#endif

    }
     
    
#ifdef DEBUG_MOVE
    cout<<parallel.rank()<<"; move : end of buffer 2 and 5 copy"<<endl;
#endif
    
    
    pTemp = sendBuffer[0];
    sendBuffer[0]=recBuffer[0];
    recBuffer[0]=pTemp;
    if(bufferSize[0]!=0)delete[] recBuffer[0];
    bufferSize[0]=bufferSizeRec[0];
        
    pTemp = sendBuffer[1];
    sendBuffer[1]=recBuffer[3];
    recBuffer[3]=pTemp;
    if(bufferSize[1]!=0)delete[] recBuffer[3];
    bufferSize[1]=bufferSizeRec[3];
    
    pTemp = sendBuffer[3];
    sendBuffer[3]=recBuffer[1];
    recBuffer[1]=pTemp;
    if(bufferSize[3]!=0)delete[] recBuffer[1];
    bufferSize[3]=bufferSizeRec[1];
    
    pTemp = sendBuffer[4];
    sendBuffer[4]=recBuffer[4];
    recBuffer[4]=pTemp;
    if(bufferSize[4]!=0)delete[] recBuffer[4];
    bufferSize[4]=bufferSizeRec[4];

    if(bufferSize[2]!=0)delete[] sendBuffer[2];
    if(bufferSizeRec[2]!=0)delete[] recBuffer[2];
    if(bufferSize[5]!=0)delete[] sendBuffer[5];
    if(bufferSizeRec[5]!=0)delete[] recBuffer[5];
    
    
    
    //pack list 6 & 7 into buffer 2 & 5
    
    bufferSize[2]=part_moveProc[6].size();
    if( bufferSize[2]!=0 )
    {
        sendBuffer[2] = new part[bufferSize[2]];
        for(it=part_moveProc[6].begin(),p=0; it != part_moveProc[6].end();++it,p++)sendBuffer[2][p]=(*it);
        part_moveProc[6].clear();
    }
    
    bufferSize[5]=part_moveProc[7].size();
    if( bufferSize[5]!=0 )
    {
        sendBuffer[5] = new part[bufferSize[5]];
        for(it=part_moveProc[7].begin(),p=0; it != part_moveProc[7].end();++it,p++)sendBuffer[5][p]=(*it);
        part_moveProc[7].clear();
    }
    
   //send z
    
    
    if(parallel.grid_rank()[0]%2==0)
	{
		////////////////////////////////
		///send/rec to/from higher rank
		////////////////////////////////
		
		if(parallel.grid_rank()[0]!=parallel.grid_size()[0]-1)/// si pas le dernier alors envoie au +1
		{
            //send
            parallel.send_dim0( &bufferSize[3], 3, parallel.grid_rank()[0]+1);
            for(int i=3;i<6;i++)
            {   
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[0]+1);
                }
                
            }
            
            //recieve
            parallel.receive_dim0( &bufferSizeRec[3], 3, parallel.grid_rank()[0]+1);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[0]+1);
                    
                }
            }
		}
		
		//////////////////////////////
		///send/rec to/from lower rank
		//////////////////////////////
		
		if(parallel.grid_rank()[0] != 0)     /// si pas le premier alors envoie au -1
		{
            //send
            parallel.send_dim0( bufferSize, 3, parallel.grid_rank()[0]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[0]-1);
                }
            }
            //recieve
            parallel.receive_dim0( bufferSizeRec, 3, parallel.grid_rank()[0]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                   recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[0]-1);
                    
                }
            }
            
		}
		else if(parallel.grid_size()[0]%2==0)  /// si pair et = 0 alors envoi au dernier
		{
            //send
            parallel.send_dim0( bufferSize, 3, parallel.grid_size()[0]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , parallel.grid_size()[0]-1);
                }
            }
            //recieve
            parallel.receive_dim0( bufferSizeRec, 3, parallel.grid_size()[0]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , parallel.grid_size()[0]-1);
                    
                }
            }
		}
        
	}
	else
	{
		//////////////////////////////
		///rec/send from/to lower rank
		//////////////////////////////
        
		//tous recois du -1 puis envois au -1
		
        //recieve
        parallel.receive_dim0( bufferSizeRec, 3, parallel.grid_rank()[0]-1);
        for(int i=0;i<3;i++)
        {
            if(bufferSizeRec[i]!=0)
            {
                recBuffer[i]=new part[bufferSizeRec[i]];
                parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[0]-1);
                
            }
        }
        //send
		parallel.send_dim0( bufferSize, 3, parallel.grid_rank()[0]-1);
        for(int i=0;i<3;i++)
        {
            if(bufferSize[i]!=0)
            {
                parallel.send_dim0( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[0]-1);
            }
		}
        
        
		//////////////////////////////
		//rec/send from/to higher rank
		/////////////////////////////
        
        if(parallel.grid_rank()[0]!=parallel.grid_size()[0]-1)//si pas dernier alors recoi du +1
		{
            
            //recieve
            parallel.receive_dim0( &bufferSizeRec[3], 3, parallel.grid_rank()[0]+1);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , parallel.grid_rank()[0]+1);
                    
                }
            }
            
            //send
            parallel.send_dim0( &bufferSize[3], 3, parallel.grid_rank()[0]+1);
            for(int i=3;i<6;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , parallel.grid_rank()[0]+1);
                }
            }
            
		}
		else // si dernier alors pair et recoi du premier
		{
            
            //recieve
            parallel.receive_dim0( &bufferSizeRec[3], 3, 0);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , 0);
                }
            }
            
            //send
            parallel.send_dim0( &bufferSize[3], 3,0);
            for(int i=3;i<6;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , 0);
                }
			}
		}
        
	}
	
	//if unpair :
	/////////////
	
	if(parallel.grid_size()[0]%2!=0)
	{
		if(parallel.grid_rank()[0]==0)
		{
            //send
            parallel.send_dim0( bufferSize, 3, parallel.grid_size()[0]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , parallel.grid_size()[0]-1);
                }
            }
            //recieve
            parallel.receive_dim0( bufferSizeRec, 3, parallel.grid_size()[0]-1);
            for(int i=0;i<3;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , parallel.grid_size()[0]-1);
                    
                }
            }
		}
		if(parallel.grid_rank()[0]==parallel.grid_size()[0]-1)
		{
            
            //recieve
            parallel.receive_dim0( &bufferSizeRec[3], 3, 0);
            for(int i=3;i<6;i++)
            {
                if(bufferSizeRec[i]!=0)
                {
                    recBuffer[i]=new part[bufferSizeRec[i]];
                    parallel.receive_dim0( recBuffer[i], bufferSizeRec[i] , 0);
                    
                }
            }
            
            //send
            parallel.send_dim0( &bufferSize[3], 3,0);
            for(int i=3;i<6;i++)
            {
                if(bufferSize[i]!=0)
                {
                    parallel.send_dim0( sendBuffer[i], bufferSize[i] , 0);
                }
            }
		}
	}
    
    //unpack data
	for(int i=0;i<6;i++)
    {
	    numParticles_ += bufferSizeRec[i] ;
    }
    
    for(int p=0;p<6;p++)
    {
        for(int i=0;i<bufferSizeRec[p];i++)
        {
            this->getPartCoordLocal(recBuffer[p][i],newLocalCoord);
            
            x.setCoordLocal(newLocalCoord);
            
            field_part_(x).size += 1;
            field_part_(x).parts.push_back(recBuffer[p][i]);
        }
    }
    
     
    for(int i=0;i<6;i++)
    {
        if(bufferSize[i]!=0)delete[] sendBuffer[i];
        if(bufferSizeRec[i]!=0)delete[] recBuffer[i];
    }
    delete[] sendBuffer;
    delete[] recBuffer;
    
    if(nfields!=0) delete[] sites;
    
    
}


template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::saveHDF5(string filename_base, int fileNumber)
{
  if(parallel.grid_size()[0] % fileNumber != 0)
    {
      COUT<<"fileNumber need to be a devider of parallel.grid_size()[0], aborting..."<<endl;
      exit(111);
    }
  string filename;
    
    
    
  int numProcPerFile = parallel.size()/fileNumber;
  int numProcPerFileDim0 = parallel.grid_size()[0]/fileNumber;
  int whichFile  = parallel.grid_rank()[0] * fileNumber / parallel.grid_size()[0];
  //int rankInFile;
  long numParts[numProcPerFile];
  int ranksList[numProcPerFile]; 
  MPI_Comm fileComm;
  MPI_Group fileGroup;
  part * partlist;
  partlist = new part[numParticles_];
  long index;

  LATfield2::Site x(lat_part_);
  typename std::list<part>::iterator it;
       
  for(int m=0;m<parallel.grid_size()[1];m++)
    {
      for(int n = 0;n<numProcPerFileDim0;n++)
        {
	  ranksList[m+(n*parallel.grid_size()[1])]= parallel.grid2world(n + whichFile * numProcPerFileDim0,m);
        }
    }
    
  MPI_Group_incl(parallel.lat_world_group(),numProcPerFile,ranksList,&fileGroup);
  MPI_Comm_create(parallel.lat_world_comm(),fileGroup, &fileComm);
  //MPI_Group_rank(fileGroup, &rankInFile);
    
  /*  
  numParts[rankInFile] = numParticles_;
  for(int i=0;i<numProcPerFile;i++)
    {
      MPI_Bcast(&numParts[i],1,MPI_LONG,i,fileComm);
    }
  */

  index=0;
  for(x.first();x.test();x.next())
    {
      if(field_part_(x).size!=0)
        {
	  for(it=field_part_(x).parts.begin(); it != field_part_(x).parts.end();++it)
            {

	      for(int i=0;i<3;i++)
                {
		  partlist[index]=(*it);
                }
	      index++;
            }
        }
    }
  
  fileDsc fd;
  fd.fileNumber=fileNumber;
  fd.numParts=numParticles_;
  fd.world_size=parallel.size();
  fd.grid_size[0]=parallel.grid_size()[0];
  fd.grid_size[1]=parallel.grid_size()[1];
  fd.numProcPerFile =numProcPerFile;
    for(int i=0;i<3;i++){
        fd.boxSize[i] = boxSize_[i];
        fd.localBoxSize[i] = lat_resolution_ * (Real)(lat_part_.sizeLocal(i));
        fd.latSize[i]=lat_part_.size(i);
    }
    fd.localBoxOffset[0] = 0;
    fd.localBoxOffset[1] = lat_resolution_ * (Real)(lat_part_.coordSkip()[1]);
    fd.localBoxOffset[2] = lat_resolution_ * (Real)(lat_part_.coordSkip()[0]);
    
    
    Real fileBoxSize[fileNumber];
    for(int i=0;i<fileNumber;i++)fileBoxSize[i]=0;
    fileBoxSize[whichFile]=fd.localBoxSize[2];
    parallel.sum_dim0(fileBoxSize,fileNumber);
    
    //cout<<fileBoxSize[whichFile] <<"  ///"<<endl;
    
    Real fileBoxOffset[fileNumber];
    for(int i=0;i<fileNumber;i++)fileBoxOffset[i]=boxSize_[2]+1.;
    fileBoxOffset[whichFile]=fd.localBoxOffset[2];
    parallel.min_dim0(fileBoxOffset,fileNumber);
    
    //for(int i=0;i<fileNumber;i++)cout<<fileBoxOffset[i]<<" ... "<<endl;
    
    fd.fileBoxSize = fileBoxSize[whichFile];
    fd.fileBoxOffset = fileBoxOffset[whichFile];
    
    
    
  if(fileNumber==1)  filename = filename_base +".h5";
  else filename = filename_base + "_" + LATfield2::int2string(whichFile,999)+".h5";


  save_hdf5_particles(filename,
		      partlist,
		      part_global_info_,
		      part_datatype_,
		      fd,
		      fileComm);



  MPI_Comm_free(&fileComm);
  MPI_Group_free(&fileGroup);

  delete[] partlist;

}


template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::loadHDF5(string filename_base, int fileNumber)
{
    //get_fd_global.
    struct fileDsc fd[fileNumber];
   
    part * partList;
    long partList_size,partList_offset;
    
    string filename;
    part_info part_info_file;
    
    std::list<int> file_list;
    std::list<int> block_list;
    std::list<int>::iterator it;
    
    std::list<fileDsc> fileDscList;
    
    long * numParts_file;
    Real * localBoxOffset_file;
    Real * localBoxSize_file;
   
        if(fileNumber ==1)
        {
            get_fileDsc_global(filename_base + ".h5",fd[0]);
            get_partInfo(filename_base + ".h5",part_info_file,part_datatype_);
        }
        else
        {
            for(int i=0;i<fileNumber;i++)get_fileDsc_global(filename_base + "_" + LATfield2::int2string(i,999)+".h5",fd[i]);
            get_partInfo(filename_base + "_" + LATfield2::int2string(0,999)+".h5",part_info_file,part_datatype_);
        }
    
    
    
    if(fd[0].boxSize[0]!=boxSize_[0] || fd[0].boxSize[1]!=boxSize_[1] || fd[0].boxSize[2]!=boxSize_[2]){
        cout<<"LATfield2::Particles::loadHDF5  :  wrong boxSize, exiting"<<endl;
        exit(-111);
    }
    
    //verif if it is correct parlicles.
   
    string type_name_file = part_info_file.type_name;
    string type_name_mem = part_global_info_.type_name;
     
    if(type_name_file.compare(type_name_file))
    {
        cout<<"LATfield2::Particles::loadHDF5  :  wrong particles type, expecting "<<type_name_mem <<", exiting"<<endl;
        exit(-111);
    }
    
    part_global_info_ = part_info_file;
    
    //compute which file to read.
    
    Real localBoxOffset[3];
    Real localBoxSize[3];
    
    localBoxOffset[0] = 0;
    localBoxOffset[1] = lat_resolution_ * (Real)(lat_part_.coordSkip()[1]);
    localBoxOffset[2] = lat_resolution_ * (Real)(lat_part_.coordSkip()[0]);
    
    for(int i=0;i<3;i++)
    {
        localBoxSize[i] = lat_resolution_ * (Real)(lat_part_.sizeLocal(i));
    }
    
    for(int i=0;i<fileNumber;i++)
    {
        //cout<<"offset: "<<fd[i].fileBoxOffset<<endl;
        //cout<<"box: "<<fd[i].fileBoxSize<<endl;
        
        if( !(localBoxOffset[2] >= fd[i].fileBoxOffset+fd[i].fileBoxSize) &&
           !(fd[i].fileBoxOffset >= localBoxOffset[2] + localBoxSize[2]) ) file_list.push_back(i);
        
    }
    
    for(it = file_list.begin();it != file_list.end(); it++)
    {

        //load the file block list.
        numParts_file = new long[fd[(*it)].numProcPerFile];
        localBoxOffset_file = new Real[3 * fd[(*it)].numProcPerFile];
        localBoxSize_file = new Real[3 * fd[(*it)].numProcPerFile];
        
        if(fileNumber ==1)get_fileDsc_local(filename_base + ".h5",numParts_file,
                                            localBoxOffset_file,localBoxSize_file);
        else get_fileDsc_local(filename_base + "_" + LATfield2::int2string((*it),999)+".h5",
                               numParts_file,localBoxOffset_file,localBoxSize_file);
            
        
        //look if need to reed to block, if yes read it and add particles...
        for(int i=0;i<fd[(*it)].numProcPerFile;i++)
        {
            if( !(localBoxOffset[2] >= localBoxOffset_file[3*i+2] + localBoxSize_file[3*i+2]) &&
                !(localBoxOffset_file[3*i+2] >= localBoxOffset[2] + localBoxSize[2])  && 
               !(localBoxOffset[1] >= localBoxOffset_file[3*i+1] + localBoxSize_file[3*i+1]) &&
               !(localBoxOffset_file[3*i+1] >= localBoxOffset[1] + localBoxSize[1])  ){
                
                cout<< parallel.grid_rank()[0]<<";"<< parallel.grid_rank()[1] <<"ok:" <<*it<<" , "<< i <<endl;   
                
                //load the particles list...
                partList_size = numParts_file[i];
                partList_offset=0;
                for(int l=0;l<i;l++)partList_offset += numParts_file[l];
                partList = new part[partList_size];
                
                cout<< "list size: "<<partList_size <<endl;
                
                if(fileNumber ==1)get_part_sublist(filename_base + ".h5",
                                                   partList_offset,partList_size,partList,part_datatype_);
                else get_part_sublist(filename_base + "_" + LATfield2::int2string((*it),999)+".h5",
                                      partList_offset,partList_size,partList,part_datatype_);
                
                
                for(int p=0;p<partList_size;p++)
                {
                    this->addParticle_global(partList[p]);
                    
                    //cout<< "adding part: "<<partList[p].ID<<endl;
                }
                
                
                delete[] partList;
            }
        }
        
        delete[] numParts_file;
        delete[] localBoxOffset_file;
        delete[] localBoxSize_file;
        cout<< parallel.grid_rank()[0]<<";"<< parallel.grid_rank()[1] <<" ,  nparts: "<< numParticles_ << endl;
        
        
    }
    
    
}


#endif
