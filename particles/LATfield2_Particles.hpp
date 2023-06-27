#ifndef LATFIELD2_PARTICLES_HPP
#define LATFIELD2_PARTICLES_HPP



/**
 * \defgroup partModule Particles Module
 * @{
 */

/*! \brief Projections for scalars vectors and tensors using hybrid cloud-in-cell and nearest-grid-point


 */

/**
 * \defgroup projections Projections
 * @{
 */
/**@}*/


/*! \brief The Particles class maps particles on a Lattice and manages particle displacements.


 */
/**
 * \defgroup prartClass Particles Class
 * @{
 */
/**@}*/


/*! \brief description of a particle type.


 */

/**
 * \defgroup partdesc Particles description
 * @{
 */
/**@}*/


/**@}*/



/*! \file LATfield2_Particles.hpp
 \brief LATfield2_Particles.hpp contain the class Particles definition.
 \author David Daverio
 */



#ifndef DOXYGEN_SHOULD_SKIP_THIS
#include "LATfield2_particle_simple.hpp"
#include "LATfield2_particle_rk4.hpp"
#include "particles_tools.hpp"
#endif
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

#define MAX_NUMBER 9223372036854775807

using namespace LATfield2;

template <typename part, typename part_info, typename part_dataType>
class Particles;
#include "projections.hpp"

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template <typename  part>
struct  partList{
//    int size;
    std::forward_list<part>  parts;
//    std::list<part>  partsTemp;
    partList() : parts(){} //size(0), parts(), partsTemp(){}
};
#endif
/**
 * \addtogroup prartClass
 * @{
 */



/*!
 \brief function to set the name of a particle type.

*/
template<typename pInfo>
void set_parts_typename(pInfo *info, string type_name)
{

    info->type_name_size = type_name.size();
    strcpy(info->type_name,type_name.c_str());

}

/*! \class Particles
 \brief the Particles class aims to map a list of particles of a given type and to manage the displacement of those particles.

 The Particles class is a template class. It take as template the 3 structure which describe a particle type. The class maps the particles to a Lattice object and manages the displacement of the particles.

 */
template <typename part, typename part_info, typename part_dataType>
class Particles
{

public:
    //! Constructor.
  Particles(){;};
    //! destructor.
  ~Particles();

    /*!
     Initialization.
     \param part_info part_global_info : structure containing the global properties of the particles.
     \param part_dataType part_datatype : structure containing the datatype of every particle property.
     \param Lattice * lat_part : Lattice on which the particles are maped, must be a 3d lattice.
     \param Real boxSize[3] : size of each dimension of the Lattice in unit used for the particle positions. The resolution boxSize[i]/lat_part.size(i) need to be the same for each i (for each dimension)).
     */
  void initialize(part_info part_global_info,
		  part_dataType part_datatype,
		  Lattice * lat_part,
		  Real boxSize[3]);

    /*!
     Method to get the process which store a given particle.
     \param part pcl: Input: structure containing the individual property of a particle linked to this instance.
     \param int * ranks:  Output: pointer to an array of two integer. The position of the process containing the particle "part" in the compute process grid.
     */
  void getPartProcess(part pcl,int * ranks);
    /*!
     Method to get the process which store a given particle, but this one does not take into account the fact that the process grid is a torus. (meaning it can return a number smaller than 1 or bigger than grid_size-1). Used for the displacement of the particles.
     \param part pcl: Input: structure containing the individual property of a particle linked to this instance.
     \param int * ranks:  Output: pointer to an array of two integer. The position of the process containing the particle "part" in the compute process grid.
     */
  void getPartNewProcess(part pcl,int * ranks);

    /*!
     Method to get the coordinate of the cell containing a given particle.
     \param part pcl: Input: structure containing the individual property of a particle linked to this instance.
     \param int * coord:  Output: pointer to an array of 3 integer. The (global) coordinate in the lattice of the cell which contains the particle.
     */
  void getPartCoord(part pcl,int * coord);
    /*!
     Method to get the local coordinate of the cell containing a given particle. (the 0,0,0 is the lowest cells of the local part of the lattice)
     \param part pcl: Input: structure containing the individual property of a particle linked to this instance.
     \param int * coord:  Output: pointer to an array of 3 integer. The (global) coordinate in the lattice of the cell which contains the particle.
     */
  void getPartCoordLocal(part pcl,int * coord);
    /*!
     Method to add a particle. Glabal function, each process will call it. The particle will be added only in one process.

     \param part newPart: Particle to add (structure containing the individual property).
     */
  bool addParticle_global(part newPart);

  void prepare_RK();


    template<typename mappingClass>
    Real updateVel(Real (*updateVel_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,mappingClass *,int,double*,double*,int),
                   double dtau,
                   Field<Real> ** fields=NULL,
                   int nfields=0,
                   mappingClass * mc = NULL,
                   double * params=NULL,
                   double * output=NULL,
                   int * reduce_type=NULL,
                   int noutput=0);

    /*!
    Method to modify the velocity of the particle. This method can be used to modify any individual property of a particles.


    \param *updateVel_funct
    \param double dtau: variation of time.
    \param Field<Real> ** fields=NULL: array of pointer to field class.
    \param int nfields: size of the array fields.
    \param double * params: pointer to an array of double, used to pass constants.
    \param double * output: pointer to an array of double. This array is used to return statistics over the particles properties, The outputs are constructed within the function updateVel_funct, and then reduced over every particles. The reduction can be the sum, the minimum or the maximum over all particles or over the particles stored in this given process.
    \param int * reduce_type: array with same size of the output array. This array is used to specify the reduction type which can be: SUM,MIN,MAX,SUM_LOCAL,MIN_LOCAL,MAX_LOCAL
    \param int noutput: size of the arrays output and reduce_type.

    */

    Real updateVel(Real (*updateVel_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,int,double*,double*,int),
                   double dtau,
                   Field<Real> ** fields=NULL,
                   int nfields=0,
                   double * params=NULL,
                   double * output=NULL,
                   int * reduce_type=NULL,
                   int noutput=0);

/*    template<typename mappingClass>
    void moveParticles( void (*move_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,mappingClass *,int,double*,double*,int),
                       double dtau,
                       Field<Real> ** fields=NULL,
                       int nfields=0,
                       mappingClass * mc = NULL,
                       double * params=NULL,
                       double * output=NULL,
                       int * reduce_type=NULL,
                       int noutput=0);*/

    void moveParticles( void (*move_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,int,double*,double*,int),
                       double dtau,
                       Field<Real> ** fields=NULL,
                       int nfields=0,
                       double * params=NULL,
                       double * output=NULL,
                       int * reduce_type=NULL,
                       int noutput=0);

#ifdef HDF5
    /*!
     Method to save all particles of this instance using HDF5 data format.

     \param string filename_base: base name of the file. Must contains the complete path to the file. It should not contains any extension (automatically added "_XXX.h5")
     \param int fileNumber: Number of file which will be written Files written : filename_base_000.h5 to filename_base_fileNumber.h5)
     */
  void saveHDF5(string filename_base, int fileNumber);
    /*!
     Method to load particles to this instance using HDF5 data format. Lattice size has not to be the same than the one used to write the files and the number of processes neither. But the size of the lattice in unit used for the particle positions has to be the same (no possibility to zoom).

     \param string filename_base: base name of the files. Must contains the complete path to the file. It should not contains any extension (automatically added ".h5") or file number (file : filename_base_XXX.h5).
     \param int fileNumber: Number of files to read.
     */
  void loadHDF5(string filename_base, int fileNumber);
#endif
#ifdef EXTERNAL_IO
    /*!
     Method to open a particle file using the output server. (output only).

     \param string filename_base: base name of the file. Must contains the complete path to the file. It should not contains any extension (automatically added "_XXX.h5")
     */
    void saveHDF5_server_open(string filename_base);
    /*!
     Method to write data in a particle file using the output server. (output only). The file should have been opened using saveHDF5_server_open(...). But if not already open, the method will open a file with default filename: "defaultfilename".

     \param string filename_base: base name of the file. Must contains the complete path to the file. It should not contains any extension (automatically added "_XXX.h5")
     */
    void saveHDF5_server_write(string filename_base  = "defaultfilename");
#endif
    /*!
     Method to "cout" the particle with a given ID, should never be used expect for debbuging as very ineficient!

     \param long ID: ID of the particle which will be "cout"
     */
    void coutPart(long ID);

    /*!
     Method to get the Lattice on which the particles are maped.
     \return lat_part_
     */
    Lattice & lattice(){return lat_part_;};
    /*!
     Method to get the Field in which the particles lists are strored.
     \return field_part_
     */
    Field<partList<part> > & field(){return field_part_;};
    /*!
     Method to get the resolution of a cell, in units used for the particle positions.
     \return field_part_
     */
    Real res(){return lat_resolution_;};

    /*!
     Method to get the global properties of the particles type of this instance of the class Particles.
     \return part_global_info_: structure containing the global properties of the particles.
     */
    part_info * parts_info();

    /*!
     Method to get the mass type. The mass can be a global property (GLOBAL_MASS), a individual property (INDIVIDUAL_MASS) or can be not defined (NO_MASS).
     \return mass_type_
     */
    int mass_type(){return mass_type_;};
    /*!
     Method to get the offset of the mass property within its respective structure.
     \return mass_offset_
     */
    size_t mass_offset(){return mass_offset_;};

    long numParticles(){
      long temp = numParticles_;
      parallel.sum(temp);
      return temp;
    };

    void cout_particle_velocity_stats(const string name);

protected:

  part_info part_global_info_;
  part_dataType part_datatype_;

  Lattice lat_part_;
  Real  lat_resolution_;
  Real boxSize_[3];

  Field<partList<part>> field_part_;

  int mass_type_;
  size_t mass_offset_;

  long numParticles_;

#ifdef EXTERNAL_IO
    ioserver_file io_file_;
#endif

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
        mass_type_= NO_MASS;
        mass_offset_=-1;
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
    for(int i=0;i<3;i++)coord[i] = (int)(floor(pcl.pos[i]*lat_part_.size(i))) %lat_part_.size(i);

    if(pcl.pos[2] >=boxSize_[2]) ranks[0] = parallel.grid_size()[0];
    else if (pcl.pos[2] < 0.)
    {
      if (pcl.pos[2] + boxSize_[2] >= boxSize_[2]) ranks[0] = 0;
      else ranks[0] = -1;
    }
    else ranks[0] = lat_part_.getRankDim0(coord[2]);

    if(pcl.pos[1] >=boxSize_[1]) ranks[1] = parallel.grid_size()[1];
    else if (pcl.pos[1] < 0.)
    {
      if (pcl.pos[1] + boxSize_[1] >= boxSize_[1]) ranks[1] = 0;
      else ranks[1] = -1;
    }
    else ranks[1] = lat_part_.getRankDim1(coord[1]);
}


template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::getPartCoord(part pcl,int * coord)
{
    for(int i=0;i<3;i++)coord[i] = (int)(floor(pcl.pos[i]*lat_part_.size(i))) %lat_part_.size(i);
}



template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::getPartCoordLocal(part pcl,int * coord)
{
    for(int i=0;i<3;i++)coord[i] = (int)(floor(pcl.pos[i]*lat_part_.size(i))) %lat_part_.size(i);
    coord[2]-=lat_part_.coordSkip()[0];
    coord[1]-=lat_part_.coordSkip()[1];
}

template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::coutPart(long ID)
{
    Site x(lat_part_);
    typename std::forward_list<part>::iterator it;

    for(x.first();x.test();x.next())
    {
//        if((field_part_)(x).size!=0)
//        {
          for(it=(field_part_)(x).parts.begin(); it != (field_part_)(x).parts.end();++it)
          {
              if((*it).ID==ID)cout<< "Parallel ranks: ("<<parallel.grid_rank()[1]<<","<<parallel.grid_rank()[0]<<") ; "<< part_global_info_.type_name<<": "<<*it<< " , lattice position:"<<x <<endl;
          }
//        }

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
      //field_part_(x).size += 1;
      field_part_(x).parts.push_front(newPart);
      numParticles_ +=1;
      return true;
    }
  else
    {
      return false;
    }
}

template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::cout_particle_velocity_stats(const string text)
{
  Site  xPart(lat_part_);
  typename std::forward_list<part>::iterator it;

  double min[3] = {1000000000,1000000000,1000000000};
  double max[3] = {-1000000000,-1000000000,-1000000000};
  long count = 0;
  double mean[3] = {0.0,0.0,0.0};

  for(xPart.first() ; xPart.test(); xPart.next())
  {
    //if(field_part_(xPart).size!=0)
    //{
      for (it=(field_part_)(xPart).parts.begin(); it != (field_part_)(xPart).parts.end(); ++it)
      {
        for(int i=0;i<3;i++)
        {
          mean[i] += (*it).vel[i];
          if((*it).vel[i]<min[i])min[i] = (*it).vel[i];
          if((*it).vel[i]>max[i])max[i] = (*it).vel[i];
        }
        count++;
      }
    //}
  }

  parallel.min(min,3);
  parallel.max(max,3);
  parallel.sum(count);
  parallel.sum(mean,3);
  for(int i=0;i<3;i++)mean[i] /= count;

  COUT<<"--- "<<text<<" velocity stats ---"<<setprecision(10)<<endl;
  for(int i=0;i<3;i++)
  {
    COUT<<"comp "<<i<< " : min "<<min[i]<< " ; max "<<max[i]<< " ; mean "<<mean[i]<<endl;
  }
  COUT<<"----------------------"<<setprecision(6)<<endl;

}

template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::prepare_RK()
{
  Site  xPart(lat_part_);
  typename std::forward_list<part>::iterator it;

  for(xPart.first() ; xPart.test(); xPart.next())
  {
    //if(field_part_(xPart).size!=0)
    //{
      for (it=(field_part_)(xPart).parts.begin(); it != (field_part_)(xPart).parts.end(); ++it)
      {
        for (int l=0; l<3; l++)
        {
          (*it).pos_in[l] = (*it).pos[l];
          (*it).pos_out[l] = (*it).pos[l];
        }
      }
    //}
  }

}

template <typename part, typename part_info, typename part_dataType>
template <typename mappingClass>
Real Particles<part,part_info,part_dataType>::updateVel(Real (*updateVel_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,mappingClass *,int,double*,double*,int),
               double dtau,
               Field<Real> ** fields,
               int nfields,
               mappingClass * mc,
               double * params,
               double * output,
               int * reduce_type,
               int noutput)
{
  Site  xPart(lat_part_);
  Site * sites = NULL;

  if(nfields!=0)
  {
      sites = new LATfield2::Site[nfields];
      for(int i = 0;i<nfields;i++)
      {
          sites[i].initialize(fields[i]->lattice());
          sites[i].first();
      }
  }

  typename std::forward_list<part>::iterator it;
  double frac[3];
  Real x0;
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

  for(xPart.first() ; xPart.test(); xPart.next())
  {
      //if(field_part_(xPart).size!=0)
      //{
          for (it=(field_part_)(xPart).parts.begin(); it != (field_part_)(xPart).parts.end(); ++it)
          {
              //old fashion uncompatible with runge kutta, frac is in [0,1], but shoud be allowed to be in [-1,2]
              //to take into account displacements during the runge kutta steps....
              //when using runge kutta, one have to use the proper projections methods...!!!
              //what a crap:
              // new method still compatible with previous methodology.... but! frac is still the offeset of the particles
              // in respect to the lowest corner of the cell where the particle is. (but in the runge kutta step, it is the cells
              // where the particle is at the beginning of the runge kutta (NOT THE ONE OF THE RK SUBSTEP!!!!))

              //for (int l=0; l<3; l++)
              //    frac[l] = modf( (*it).pos[l] / lat_resolution_, &x0);
              for (int l=0; l<3; l++)
              {
                frac[l] =  (*it).pos[l]*lat_part_.size(l) - xPart.coord(l);
              }



              v2 = updateVel_funct(dtau,
                         lat_resolution_,
                         &(*it),
                         frac,
                         part_global_info_,
                         fields,
                         sites,
                         mc,
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
      //}
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
  if(nfields>0) delete[] sites;

  return sqrt(maxvel);

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

    Site  xPart(lat_part_);
    Site * sites = NULL;

    if(nfields!=0)
    {
        sites = new LATfield2::Site[nfields];
        for(int i = 0;i<nfields;i++)
        {
            sites[i].initialize(fields[i]->lattice());
            sites[i].first();
        }
    }

    typename std::forward_list<part>::iterator it;
    double frac[3];
    Real x0;
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

    for(xPart.first() ; xPart.test(); xPart.next())
    {
        //if(field_part_(xPart).size!=0)
        //{
            for (it=(field_part_)(xPart).parts.begin(); it != (field_part_)(xPart).parts.end(); ++it)
            {
                //old fashion uncompatible with runge kutta, frac is in [0,1], but shoud be allowed to be in [-1,2]
                //to take into account displacements during the runge kutta steps....
                //when using runge kutta, one have to use the proper projections methods...!!!
                //what a crap:
                // new method still compatible with previous methodology.... but! frac is still the offeset of the particles
                // in respect to the lowest corner of the cell where the particle is. (but in the runge kutta step, it is the cells
                // where the particle is at the beginning of the runge kutta (NOT THE ONE OF THE RK SUBSTEP!!!!))

                //for (int l=0; l<3; l++)
                //    frac[l] = modf( (*it).pos[l] / lat_resolution_, &x0);
                for (int l=0; l<3; l++)
                {
                  frac[l] =  (*it).pos[l]*lat_part_.size(l) - xPart.coord(l);
                }



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
        //}
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
    if(nfields>0) delete[] sites;

    return sqrt(maxvel);


}

/*template <typename part, typename part_info, typename part_dataType>
template<typename mappingClass>
void Particles<part,part_info,part_dataType>::moveParticles( void (*move_funct)(double,double,part*,double *,part_info,Field<Real> **,Site *,mappingClass *,int,double*,double*,int),
                   double dtau,
                   Field<Real> ** fields,
                   int nfields,
                   mappingClass * mc,
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
      LATfield2::Site * sites = NULL;

      typename std::list<part>::iterator it,itTemp;
      //Real b[3];
      double frac[3];
      Real x0;
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
              output[i]=MAX_NUMBER;
              //COUT<<"min"<<endl;
          }
          else if(reduce_type[i] & (MAX | MAX_LOCAL))
          {
              output[i]=-MAX_NUMBER;
              //COUT<<"max"<<endl;
          }
      }

      part partTest;

      for(x.first();x.test();x.next())
      {
          if(field_part_(x).size!=0)
          {
              for(int i=0;i<3;i++)localCoord[i] = x.coordLocal(i);

              for(it=field_part_(x).parts.begin(); it != field_part_(x).parts.end();)
              {
                  itTemp = it;
                  ++it;

                  partTest = (*itTemp);

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
                             mc,
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
                          cout<< "particle position: "<< (*itTemp) <<endl;
                          cout<< "particle position old: "<< partTest <<endl;
                          cout<<"particle : "<<(*itTemp).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
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
                          cout<< "particle position: "<< (*itTemp) <<endl;
                          cout<< "particle position old: "<< partTest <<endl;
                          cout<<"particle : "<<(*itTemp).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
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
                          cout<< "particle : "<<(*itTemp).ID<<" has moved too much (more than 1 proc)."<<endl;
                          cout<< "particle position: "<< (*itTemp) <<endl;
                          cout<< "particle position old: "<< partTest <<endl;
                          cout<<"particle : "<<(*itTemp).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
                      }
                  }
                  else
                  {
                      cout<< "particle : "<<(*itTemp).ID<<" has moved too much (more than 1 proc)."<<endl;
                      cout<< "particle position: "<< (*itTemp) <<endl;
                      cout<< "particle position old: "<< partTest <<endl;
                      cout<<"particle : "<<(*itTemp).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
                  }

              }
          }

          if(nfields!=0) for(int i=0;i<nfields;i++) sites[i].next();

      }




      for(x.first();x.test();x.next())if((field_part_)(x).size!=0)(field_part_)(x).parts.splice((field_part_)(x).parts.end(), (field_part_)(x).partsTemp );


       //cout<<"starting first dim"<<endl;


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
      // this is wrong: sites gets deleted again later
      //    if(nfields!=0) {
      //      if(sites) delete[] sites;
      //      else { std::cerr << "WTF, nfields != 0, but sites is not initialized.\n"; exit(-1); }
      //    }

      //remove the number of part...
      for(int i=0;i<8;i++)
      {
          numParticles_ -=  part_moveProc[i].size();
      }



      ////////////////////////////////////////////////////////////////
      //First send Y direction


      //cout<<"okokok  move firs statrt pack"<<endl;
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


       //cout<<"okokok  move firs statrt comm"<<endl;

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


      //cout<<"okokok  move firs statrt unpack"<<endl;

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

      //cout<<"buffer size: "<<bufferSizeRec[5]<<endl;
      //if(bufferSizeRec[5]!=0)
      for(int i=0;i<bufferSizeRec[5];i++)
      {
          //cout<<i<<endl;

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

      //cout<<"unpack done  "<<endl;

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

    if(nfields!=0 && sites) { delete[] sites; sites = NULL; };
}*/

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
    LATfield2::Site * sites = NULL;

    typename std::forward_list<part>::iterator it,prev;
    //Real b[3];
    double frac[3];
    Real x0;
    int localCoord[3];
    int newLocalCoord[3];

//    bool flag[4];


    part **sendBuffer;
    part **recBuffer;


    sendBuffer = new part*[6];
    recBuffer = new part*[6];

    //part * pTemp;

    long p;
    long bufferSize[6];
    long bufferSizeRec[6];

    std::vector<part>  part_moveProc[8];


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
            output[i]=MAX_NUMBER;
            //COUT<<"min"<<endl;
        }
        else if(reduce_type[i] & (MAX | MAX_LOCAL))
        {
            output[i]=-MAX_NUMBER;
            //COUT<<"max"<<endl;
        }
    }

    //part partTest;
    
    for(x.first();x.test();x.next())
    {
        for(it=field_part_(x).parts.begin(); it != field_part_(x).parts.end(); it++)
        {
        	for (int l=0; l<3; l++)
                    frac[l] = modf( (*it).pos[l] * lat_part_.size(l), &x0);


                move_funct(dtau,
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
        
        if(nfields!=0) for(int i=0;i<nfields;i++) sites[i].next();
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
    
    int partRanks[2];
    int thisRanks[2];
    thisRanks[0] = parallel.grid_rank()[0];
    thisRanks[1] = parallel.grid_rank()[1];

    for(x.first();x.test();x.next())
    {
//        if(field_part_(x).size!=0)
//        {
            for(int i=0;i<3;i++)localCoord[i] = x.coordLocal(i);
            
			prev = field_part_(x).parts.before_begin();
            for(it=field_part_(x).parts.begin(); it != field_part_(x).parts.end(); it++)
            {
/*                itTemp = it;
                ++it;

                partTest = (*itTemp);

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
*/


//                int partRanks[2];
//                int thisRanks[2];
                getPartNewProcess((*it),partRanks);
//                thisRanks[0] = parallel.grid_rank()[0];
//                thisRanks[1] = parallel.grid_rank()[1];

                for(int i=0;i<3;i++)
                {
                    if((*it).pos[i]<0)it->pos[i] +=  boxSize_[i];
                    if((*it).pos[i]>=boxSize_[i])it->pos[i]-=boxSize_[i];
                }


                if(partRanks[0]==thisRanks[0] && partRanks[1]==thisRanks[1])
                {
                    int r[3];
                    getPartCoord(*it, r);
                    if(!xNew.setCoord(r))
                    {
                        cout<<"arg"<< *it <<" ; "<< partRanks[0]<< " , " << thisRanks[0] <<endl;
                    }

                    getPartCoordLocal(*it, newLocalCoord);
                    if(localCoord[0]!=newLocalCoord[0] || localCoord[1]!=newLocalCoord[1] || localCoord[2]!=newLocalCoord[2] )
                    {
                        xNew.setCoordLocal(newLocalCoord);
                        //field_part_(xNew).partsTemp.splice(field_part_(xNew).partsTemp.end(),field_part_(x).parts,itTemp);
                        field_part_(xNew).parts.splice_after(field_part_(xNew).parts.before_begin(),field_part_(x).parts,prev);
                        //field_part_(xNew).size += 1;
                        //field_part_(x).size -= 1;
                        it = prev;
                    }
                    else
                    	prev++;

                }
                else if(partRanks[1]==thisRanks[1]-1)
                {
                    if(partRanks[0]==thisRanks[0])
                    {
                        //part_moveProc[2].splice(part_moveProc[2].end(),field_part_(x).parts,itTemp);
                        part_moveProc[2].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]-1)
                    {
                        //part_moveProc[0].splice(part_moveProc[0].end(),field_part_(x).parts,itTemp);
                        part_moveProc[0].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]+1)
                    {
                        //part_moveProc[1].splice(part_moveProc[1].end(),field_part_(x).parts,itTemp);
                        part_moveProc[1].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else
                    {
                        cout<< "particle : "<<(*it).ID<<" have move to far away (more than 1 proc)."<<endl;
                        cout<< "particle position: "<< (*it) <<endl;
                        //cout<< "particle position old: "<< partTest <<endl;
                        cout<<"particle : "<<(*it).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
                    }
                    field_part_(x).parts.erase_after(prev);
                    it = prev;
                }
                else if(partRanks[1]==thisRanks[1]+1)
                {
                    if(partRanks[0]==thisRanks[0])
                    {
                        //part_moveProc[5].splice(part_moveProc[5].end(),field_part_(x).parts,itTemp);
                        part_moveProc[5].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]-1)
                    {
                        //part_moveProc[3].splice(part_moveProc[3].end(),field_part_(x).parts,itTemp);
                        part_moveProc[3].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]+1)
                    {
                        //part_moveProc[4].splice(part_moveProc[4].end(),field_part_(x).parts,itTemp);
                        part_moveProc[4].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else
                    {
                        cout<< "particle : "<<(*it).ID<<" have move to far away (more than 1 proc)."<<endl;
                        cout<< "particle position: "<< (*it) <<endl;
                        //cout<< "particle position old: "<< partTest <<endl;
                        cout<<"particle : "<<(*it).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
                    }
                    field_part_(x).parts.erase_after(prev);
                    it = prev;
                }
                else if(partRanks[1]==thisRanks[1])
                {
                    if(partRanks[0]==thisRanks[0]-1)
                    {
                        //part_moveProc[6].splice(part_moveProc[6].end(),field_part_(x).parts,itTemp);
                        part_moveProc[6].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else if(partRanks[0]==thisRanks[0]+1)
                    {
                        //part_moveProc[7].splice(part_moveProc[7].end(),field_part_(x).parts,itTemp);
                        part_moveProc[7].push_back(*it);
                        //field_part_(x).size -= 1;
                    }
                    else
                    {
                        cout<< "particle : "<<(*it).ID<<" has moved too much (more than 1 proc)."<<endl;
                        cout<< "particle position: "<< (*it) <<endl;
                        //cout<< "particle position old: "<< partTest <<endl;
                        cout<<"particle : "<<(*it).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
                    }
                    field_part_(x).parts.erase_after(prev);
                    it = prev;
                }
                else
                {
                    cout<< "particle : "<<(*it).ID<<" has moved too much (more than 1 proc)."<<endl;
                    cout<< "particle position: "<< (*it) <<endl;
                    //cout<< "particle position old: "<< partTest <<endl;
                    cout<<"particle : "<<(*it).ID<< " "<< thisRanks[0]<<" , "<< thisRanks[1]<<" , "<< partRanks[0]<<" , "<< partRanks[1]<<endl;
                }

            }
//        }

//        if(nfields!=0) for(int i=0;i<nfields;i++) sites[i].next();

    }




    //for(x.first();x.test();x.next())if((field_part_)(x).size!=0)(field_part_)(x).parts.splice((field_part_)(x).parts.end(), (field_part_)(x).partsTemp );


     //cout<<"starting first dim"<<endl;


    
    /* this is wrong: sites gets deleted again later */
    //    if(nfields!=0) {
    //      if(sites) delete[] sites;
    //      else { std::cerr << "WTF, nfields != 0, but sites is not initialized.\n"; exit(-1); }
    //    }

    //remove the number of part...
    for(int i=0;i<8;i++)
    {
        numParticles_ -=  part_moveProc[i].size();
    }



    ////////////////////////////////////////////////////////////////
    //First send Y direction


    //cout<<"okokok  move firs statrt pack"<<endl;
    //pack data
    for(int i=0;i<6;i++)
    {
        bufferSize[i]=part_moveProc[i].size();
        if( bufferSize[i]!=0 )
        {
            //sendBuffer[i] = new part[bufferSize[i]];
            //for(it=part_moveProc[i].begin(),p=0; it != part_moveProc[i].end();++it,p++)sendBuffer[i][p] = (*it);
            //part_moveProc[i].clear();
            sendBuffer[i] = part_moveProc[i].data();
        }

    }


     //cout<<"okokok  move firs statrt comm"<<endl;

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


    //cout<<"okokok  move firs statrt unpack"<<endl;

    //unpack local list: rec[2] and rec[5]

    //add partnum
    numParticles_ += bufferSizeRec[2] + bufferSizeRec[5];
    for(int i=0;i<bufferSizeRec[2];i++)
    {
        this->getPartCoordLocal(recBuffer[2][i],newLocalCoord);
        x.setCoordLocal(newLocalCoord);

        //field_part_(x).size += 1;
        //field_part_(x).parts.push_back(recBuffer[2][i]);
        field_part_(x).parts.push_front(recBuffer[2][i]);
#ifdef DEBUG_MOVE
        int verif;
        verif=addParticle_global(recBuffer[2][i]);

        if(verif != 1)
        {
            cout<<parallel.rank()<<"; MOVEBUF2 partID "<< recBuffer[2][i].ID<<" is not in the correct proc. " << recBuffer[2][i]<<endl;
        }
#endif
    }

    //cout<<"buffer size: "<<bufferSizeRec[5]<<endl;
    //if(bufferSizeRec[5]!=0)
    for(int i=0;i<bufferSizeRec[5];i++)
    {
        //cout<<i<<endl;

        this->getPartCoordLocal(recBuffer[5][i],newLocalCoord);

        x.setCoordLocal(newLocalCoord);

        //field_part_(x).size += 1;
        //field_part_(x).parts.push_back(recBuffer[5][i]);
        field_part_(x).parts.push_front(recBuffer[5][i]);
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

    //cout<<"unpack done  "<<endl;
    
    for (int i = 0; i < 6; i++)
    	part_moveProc[i].clear();

    //pTemp = sendBuffer[0];
    sendBuffer[0]=recBuffer[0];
    //recBuffer[0]=pTemp;
    //if(bufferSize[0]!=0)delete[] recBuffer[0];
    bufferSize[0]=bufferSizeRec[0];

    //pTemp = sendBuffer[1];
    sendBuffer[1]=recBuffer[3];
    //recBuffer[3]=pTemp;
    //if(bufferSize[1]!=0)delete[] recBuffer[3];
    bufferSize[1]=bufferSizeRec[3];

    //pTemp = sendBuffer[3];
    sendBuffer[3]=recBuffer[1];
    //recBuffer[1]=pTemp;
    //if(bufferSize[3]!=0)delete[] recBuffer[1];
    bufferSize[3]=bufferSizeRec[1];

    //pTemp = sendBuffer[4];
    sendBuffer[4]=recBuffer[4];
    //recBuffer[4]=pTemp;
    //if(bufferSize[4]!=0)delete[] recBuffer[4];
    bufferSize[4]=bufferSizeRec[4];

    //if(bufferSize[2]!=0)delete[] sendBuffer[2];
    if(bufferSizeRec[2]!=0)delete[] recBuffer[2];
    //if(bufferSize[5]!=0)delete[] sendBuffer[5];
    if(bufferSizeRec[5]!=0)delete[] recBuffer[5];



    //pack list 6 & 7 into buffer 2 & 5

    bufferSize[2]=part_moveProc[6].size();
    if( bufferSize[2]!=0 )
    {
        //sendBuffer[2] = new part[bufferSize[2]];
        //for(it=part_moveProc[6].begin(),p=0; it != part_moveProc[6].end();++it,p++)sendBuffer[2][p]=(*it);
        //part_moveProc[6].clear();
        sendBuffer[2] = part_moveProc[6].data();
    }


    bufferSize[5]=part_moveProc[7].size();
    if( bufferSize[5]!=0 )
    {
        //sendBuffer[5] = new part[bufferSize[5]];
        //for(it=part_moveProc[7].begin(),p=0; it != part_moveProc[7].end();++it,p++)sendBuffer[5][p]=(*it);
        //part_moveProc[7].clear();
        sendBuffer[5] = part_moveProc[7].data();
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

            //field_part_(x).size += 1;
            //field_part_(x).parts.push_back(recBuffer[p][i]);
            field_part_(x).parts.push_front(recBuffer[p][i]);
        }
    }


    for(int i=0;i<6;i++)
    {
        if(bufferSize[i]!=0 && i%3 < 2) delete[] sendBuffer[i];
        if(bufferSizeRec[i]!=0)delete[] recBuffer[i];
    }
    part_moveProc[6].clear();
    part_moveProc[7].clear();
    delete[] sendBuffer;
    delete[] recBuffer;

  if(nfields!=0 && sites) { delete[] sites; sites = NULL; };

}



#ifdef HDF5
template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::saveHDF5(string filename_base, int fileNumber)
{
    if(parallel.grid_size()[1] % fileNumber != 0)
    {
        COUT<<"fileNumber need to be a devider of parallel.grid_size()[0], aborting..."<<endl;
        exit(111);
    }
    string filename;



    int numProcPerFile = parallel.size()/fileNumber;
    int numProcPerFileDim1 = parallel.grid_size()[1]/fileNumber;
    int whichFile  = parallel.grid_rank()[1] * fileNumber / parallel.grid_size()[1];
    //int rankInFile;
    long numParts[numProcPerFile];
    int ranksList[numProcPerFile];
    MPI_Comm fileComm;
    MPI_Group fileGroup;
    part * partlist;
    partlist = new part[numParticles_];
    long index;

    LATfield2::Site x(lat_part_);
    typename std::forward_list<part>::iterator it;
    int rang[3];


    rang[0]= whichFile * numProcPerFile ;
    rang[1]= ((whichFile+1) * numProcPerFile) -1;
    rang[2]=1;
    MPI_Group_range_incl(parallel.lat_world_group(),1,&rang,&fileGroup);
    MPI_Comm_create(parallel.lat_world_comm(),fileGroup , &fileComm);

  index=0;
  for(x.first();x.test();x.next())
    {
      //if(field_part_(x).size!=0)
      //  {
	         for(it=field_part_(x).parts.begin(); it != field_part_(x).parts.end();++it)
           {
	             partlist[index]=(*it);
               index++;
           }
      //  }
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
    fileBoxSize[whichFile]=fd.localBoxSize[1];
    parallel.sum_dim1(fileBoxSize,fileNumber);



    Real fileBoxOffset[fileNumber];
    for(int i=0;i<fileNumber;i++)fileBoxOffset[i]=boxSize_[1]+1.;
    fileBoxOffset[whichFile]=fd.localBoxOffset[1];
    parallel.min_dim1(fileBoxOffset,fileNumber);

    fd.fileBoxSize = fileBoxSize[whichFile];
    fd.fileBoxOffset = fileBoxOffset[whichFile];



  if(fileNumber==1)  filename = filename_base +".h5";
  else filename = filename_base + "_" + int2string(whichFile,999)+".h5";

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
            for(int i=0;i<fileNumber;i++)get_fileDsc_global(filename_base + "_" + int2string(i,999)+".h5",fd[i]);
            get_partInfo(filename_base + "_" + int2string(0,999)+".h5",part_info_file,part_datatype_);
        }

    //cout<< fd[0].boxSize[0]<< " "<< fd[0].boxSize[1]<< " "<< fd[0].boxSize[2]<< endl;

    if(fd[0].boxSize[0]!=boxSize_[0] || fd[0].boxSize[1]!=boxSize_[1] || fd[0].boxSize[2]!=boxSize_[2]){
        cout<<"LATfield2::Particles::loadHDF5  :  wrong boxSize, exiting: "<< fd[0].boxSize[0] <<", "<< boxSize_[0] <<endl;
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

        if( !(localBoxOffset[1] >= fd[i].fileBoxOffset+fd[i].fileBoxSize) &&
           !(fd[i].fileBoxOffset >= localBoxOffset[1] + localBoxSize[1]) ) file_list.push_back(i);

    }

    for(it = file_list.begin();it != file_list.end(); it++)
    {

        //load the file block list.
	int nPP = fd[(*it)].numProcPerFile;
        numParts_file = new long[nPP];
        localBoxOffset_file = new Real[3*nPP];
        localBoxSize_file = new Real[3*nPP];

        if(fileNumber ==1)get_fileDsc_local(filename_base + ".h5",numParts_file,
                                            localBoxOffset_file,localBoxSize_file,nPP);
        else get_fileDsc_local(filename_base + "_" + int2string((*it),999)+".h5",
                               numParts_file,localBoxOffset_file,localBoxSize_file,nPP);


        //look if need to reed a block, if yes read it and add particles...
        for(int i=0;i<fd[(*it)].numProcPerFile;i++)
        {
            if( !(localBoxOffset[2] >= localBoxOffset_file[3*i+2] + localBoxSize_file[3*i+2]) &&
                !(localBoxOffset_file[3*i+2] >= localBoxOffset[2] + localBoxSize[2])  &&
               !(localBoxOffset[1] >= localBoxOffset_file[3*i+1] + localBoxSize_file[3*i+1]) &&
               !(localBoxOffset_file[3*i+1] >= localBoxOffset[1] + localBoxSize[1])  ){

                //load the particles list...
                partList_size = numParts_file[i];
                partList_offset=0;
                for(int l=0;l<i;l++)partList_offset += numParts_file[l];
                partList = new part[partList_size];

                //cout<< "list size: "<<partList_size <<endl;

                if(fileNumber ==1)get_part_sublist(filename_base + ".h5",
                                                   partList_offset,partList_size,partList,part_datatype_);
                else get_part_sublist(filename_base + "_" + int2string((*it),999)+".h5",
                                      partList_offset,partList_size,partList,part_datatype_);


                for(int p=0;p<partList_size;p++)
                {
                    this->addParticle_global(partList[p]);
                }


                delete[] partList;
            }
        }

        delete[] numParts_file;
        delete[] localBoxOffset_file;
        delete[] localBoxSize_file;
        //cout<< parallel.grid_rank()[0]<<";"<< parallel.grid_rank()[1] <<" ,  nparts: "<< numParticles_ << endl;


    }


}
#endif

#ifdef EXTERNAL_IO
template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::saveHDF5_server_open(string filename_base)
{
    io_file_ = ioserver.openFile(filename_base.c_str() ,UNSTRUCTURED_H5_FILE, part_datatype_.part_memType, part_datatype_.part_fileType);
}
template <typename part, typename part_info, typename part_dataType>
void Particles<part,part_info,part_dataType>::saveHDF5_server_write(string filename_base)
{
    if(!(io_file_.is_open))
        io_file_ = ioserver.openFile(filename_base.c_str() ,UNSTRUCTURED_H5_FILE, part_datatype_.part_memType, part_datatype_.part_fileType);


    part * partlist;
    partlist = new part[numParticles_];
    LATfield2::Site x(lat_part_);
    typename std::forward_list<part>::iterator it;
    long index=0;
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
    ioserver.sendData(io_file_,(char*)partlist,numParticles_ * H5Tget_size(part_datatype_.part_memType));
    delete[] partlist;

    hsize_t dim;
    hsize_t size[3];
    int file_number = ioserver.io_node_number();
    int numProcPerFile = parallel.size()/file_number;
    int world_size = parallel.size();
    int grid_size[2];
    grid_size[0]=parallel.grid_size()[0];
    grid_size[1]=parallel.grid_size()[1];
    int whichFile = ioserver.my_node();
    Real localBoxSize[3];
    int latSize[3];
    for(int i=0;i<3;i++){
        localBoxSize[i] = lat_resolution_ * (Real)(lat_part_.sizeLocal(i));
        latSize[i]=lat_part_.size(i);
    }

    Real localBoxOffset[3];
    localBoxOffset[0] = 0;
    localBoxOffset[1] = lat_resolution_ * (Real)(lat_part_.coordSkip()[1]);
    localBoxOffset[2] = lat_resolution_ * (Real)(lat_part_.coordSkip()[0]);


    Real fbs[file_number];
    for(int i=0;i<file_number;i++)fbs[i]=0;
    fbs[whichFile]=localBoxSize[1];
    parallel.sum_dim1(fbs,file_number);


    Real fbo[file_number];
    for(int i=0;i<file_number;i++)fbo[i]=boxSize_[1]+1.;
    fbo[whichFile]=localBoxOffset[1];
    parallel.min_dim1(fbo,file_number);


    long numPartsAll[numProcPerFile];
    Real localBoxOffsetAll[numProcPerFile*3];
    Real localBoxSizeAll[numProcPerFile*3];
    int mpi_rank = ioserver.compute_file_rank();


    numPartsAll[mpi_rank] = numParticles_;
    for(int i=0;i<3;i++)
    {
        localBoxOffsetAll[3*mpi_rank+i]=localBoxOffset[i];
        localBoxSizeAll[3*mpi_rank+i]=localBoxSize[i];
    }

    for(int i=0;i<numProcPerFile;i++)
    {
        MPI_Bcast(&numPartsAll[i],1,MPI_LONG,i,ioserver.compute_file_comm());
        MPI_Bcast(&localBoxOffsetAll[i*3],3,MPI_RealC,i,ioserver.compute_file_comm());
        MPI_Bcast(&localBoxSizeAll[i*3],3,MPI_RealC,i,ioserver.compute_file_comm());
    }

    dim =1;
    size[0]=1;
    ioserver.sendDataset(io_file_,"part_info",(char*)&part_global_info_,dim,size,part_datatype_.part_info_memType);
    ioserver.sendATTR(io_file_,"fileNumber",(char*)&file_number,1,H5T_NATIVE_INT);
    ioserver.sendATTR(io_file_,"numProcPerFile",(char*)&numProcPerFile,1,H5T_NATIVE_INT);
    ioserver.sendATTR(io_file_,"world_size",(char*)&world_size,1,H5T_NATIVE_INT);
    ioserver.sendATTR(io_file_,"grid_size",(char*)&grid_size,2,H5T_NATIVE_INT);
    ioserver.sendATTR(io_file_,"boxSize",(char*)&boxSize_,3,REAL_TYPE);
    ioserver.sendATTR(io_file_,"fileBoxSize",(char*)&fbs[whichFile],1,REAL_TYPE);
    ioserver.sendATTR(io_file_,"fileBoxOffset",(char*)&fbo[whichFile],1,REAL_TYPE);
    ioserver.sendATTR(io_file_,"latSize",(char*)latSize,3,H5T_NATIVE_INT);
    dim = 1;
    size[0] = numProcPerFile;
    ioserver.sendDataset(io_file_,"numParts",(char*)numPartsAll,dim,size,H5T_NATIVE_LONG);
    size[0] = 3;
    size[1] = numProcPerFile;
    dim = 2;
    ioserver.sendDataset(io_file_,"localBoxOffset",(char*)localBoxOffsetAll,dim,size,REAL_TYPE);
    ioserver.sendDataset(io_file_,"localBoxSize",(char*)localBoxSizeAll,dim,size,REAL_TYPE);


    ioserver.closeFile(io_file_);

}
#endif
/**@}*/




#endif
