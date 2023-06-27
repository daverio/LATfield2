#ifndef LATFIELD2_PROJECTIONS_HPP
#define LATFIELD2_PROJECTIONS_HPP



/*! \file projections.hpp
 \brief projection function for scalar, vector and tensor particle properties.

 */


/**
 * \addtogroup projections
 * @{
 */



/*! \fn inline long setIndex(long * size, long i, long j,long k)
 \return  i + size[0] *( j+size[1]*k)
 */
inline long setIndex(long * size, long i, long j,long k)
{
    return i + size[0] *( j+size[1]*k) ;
}

/*! \fn void projection_init(Field<Real> * f)
 Set to zero all components of the field f on the entire lattice (including the halo).
 Have to be performed before any projection.
 */
void projection_init(Field<Real> * f)
{
    for(long i=0;i< f->lattice().sitesLocalGross()  * (long)(f->components())  ;i++)(*f)(i)=0.0;
}



static int const FROM_PART = INDIVIDUAL_MASS;
static int const FROM_INFO = GLOBAL_MASS;

#ifndef DOXYGEN_SHOULD_SKIP_THIS

int mapv1v2[4][2][2][2] = { { {{0,6},{2,8}},{{1,7},{3,9}} } , { {{2,8},{4,10}},{{3,9},{5,11}} },
    { {{6,12},{8,14}},{{7,13},{9,15}} } , { {{8,14},{10,16}},{{9,15},{11,17}} } };

int mapv0v1[4][2][2][2] = { { {{36,37},{42,43}},{{38,39},{44,45}} } , { {{38,39},{44,45}},{{40,41},{46,47}} },
    { {{42,43},{48,49}},{{44,45},{50,51}} },{ {{44,45},{50,51} },{{46,47},{52,53}} } };

int mapv0v2[4][2][2][2] =  { { {{18,24},{19,25}},{{20,26},{21,27}} } , { {{20,26},{21,27}},{{22,28},{23,29}} },
    { {{24,30},{25,31}},{{26,32},{27,33}} } , { {{26,32},{27,33}},{{28,34},{29,35}} } };

int CIC_mapv0[2][2][2][2] =   { { {{0,1},{2,3}},{{4,5},{6,7}} } , { {{4,5},{6,7}},{{8,9},{10,11}} } };
int CIC_mapv1[2][2][2][2] =   { { {{12,14},{16,18}},{{13,15},{17,19}} } , { {{16,18},{20,22}},{{17,19},{21,23}} } };
int CIC_mapv2[2][2][2][2] =   { { {{24,28},{25,29}},{{26,30},{27,31}} } , { {{28,32},{29,33}},{{30,34},{31,35}} } };

#endif

//////Scalar projection


/*! \fn template<typename part, typename part_info, typename part_dataType> void scalarProjectionCIC_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * rho, size_t * oset = NULL,int flag_where = FROM_INFO)

 \brief cloud-in-cell scalar projection.

 Perform the projection of a scalar property of the particles. Last 2 argument are used to select the property. After the call of this method, the method scalarProjectionCIC_comm(Field<Real> * rho) have to be called to perform the reduction of the halo.

 \param Particles<part,part_info,part_dataType> * parts: bla
 \param Field<Real> * rho: pointer to the real scalar field on which the property is projected.
 \param size_t * oset: pointer to a integer which store the offset of the scalar property in the structure it belong to. If set to NULL, the projection will project the property "mass" of the particles.
 \param int flag_where: Flag to set if the property is a global or a individual one, i.e. the specify in which structure the property is. FROM_INFO set the property to be global and FROM_PART set the property to be individual.

 \sa scalarProjectionCIC_comm(Field<Real> * rho)
 \sa projection_init(Field<Real> * f)
 */
template<typename part, typename part_info, typename part_dataType>
void scalarProjectionCIC_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * rho, size_t * oset = NULL,int flag_where = FROM_INFO)
{

    if(rho->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Site xPart(parts->lattice());
    Site xField(rho->lattice());

    typename std::forward_list<part>::iterator it;

    //size_t offset;
    //*offset = oset;

    double referPos[3];
    double rescalPos[3];
    double rescalPosDown[3];
    double latresolution = parts->res();

    double mass;
    double cicVol;
    cicVol= latresolution*latresolution*latresolution;
    cicVol *= cicVol;

    Real localCube[8]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111

    int sfrom;

    size_t offset;

    if(oset == NULL)
    {
        sfrom = parts->mass_type();
        offset = parts->mass_offset();
    }
    else
    {
        sfrom =  flag_where;
        offset = *oset;
    }

    if(sfrom == FROM_INFO)
    {
        mass = *(double*)((char*)parts->parts_info() + offset);
        mass /= cicVol;
        //cout << mass<<endl;
    }




    for(xPart.first(),xField.first();xPart.test();xPart.next(),xField.next())
    {


        if(!parts->field()(xPart).parts.empty())
        {
            for(int i=0;i<3;i++)referPos[i]=xPart.coord(i)*latresolution;
            for(int i=0;i<8;i++)localCube[i]=0;

            for (it=(parts->field())(xPart).parts.begin(); it != (parts->field())(xPart).parts.end(); ++it)
            {
                for(int i =0;i<3;i++)
                {
                    rescalPos[i]=(*it).pos[i]-referPos[i];
                    rescalPosDown[i]=latresolution -rescalPos[i];
                }

                if(sfrom==FROM_PART)
                {
                    mass = *(double*)((char*)&(*it)+offset);
                    mass /=cicVol;
                }

                //000
                localCube[0] += rescalPosDown[0]*rescalPosDown[1]*rescalPosDown[2] * mass;
                //001
                localCube[1] += rescalPosDown[0]*rescalPosDown[1]*rescalPos[2] * mass;
                //010
                localCube[2] += rescalPosDown[0]*rescalPos[1]*rescalPosDown[2] * mass;
                //011
                localCube[3] += rescalPosDown[0]*rescalPos[1]*rescalPos[2] * mass;
                //100
                localCube[4] += rescalPos[0]*rescalPosDown[1]*rescalPosDown[2] * mass;
                //101
                localCube[5] += rescalPos[0]*rescalPosDown[1]*rescalPos[2] * mass;
                //110
                localCube[6] += rescalPos[0]*rescalPos[1]*rescalPosDown[2] * mass;
                //111
                localCube[7] += rescalPos[0]*rescalPos[1]*rescalPos[2] * mass;
            }

            (*rho)(xField)+=localCube[0];
            (*rho)(xField+2)+=localCube[1];
            (*rho)(xField+1)+=localCube[2];
            (*rho)(xField+1+2)+=localCube[3];
            (*rho)(xField+0)+=localCube[4];
            (*rho)(xField+0+2)+=localCube[5];
            (*rho)(xField+0+1)+=localCube[6];
            (*rho)(xField+0+1+2)+=localCube[7];
        }
    }


}



/*! \fn scalarProjectionCIC_comm(Field<Real> * rho)
 \brief communication method associated to the cloud-in-cell projection of scalars .


 \sa scalarProjectionCIC_project(...)
*/
void scalarProjectionCIC_comm(Field<Real> * rho)
{

    if(rho->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = rho->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=rho->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;


    iref = sizeLocalGross[0] - halo;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            (*rho)(setIndex(sizeLocalGross,halo,j,k)) += (*rho)(setIndex(sizeLocalGross,iref,j,k));
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1) * (long)sizeLocal[0];
    bufferSizeZ = (long)sizeLocal[0] * (long)sizeLocal[1] ;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }

    //pack data
    imax = sizeLocal[0];
    iref = sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            bufferSend[i+k*imax]=(*rho)(setIndex(sizeLocalGross,i+halo,iref,k+halo));
        }
    }
    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);
    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)(*rho)(setIndex(sizeLocalGross,i+halo,halo,k+halo))+=bufferRec[i+k*imax];

    }

    //send halo in direction Z



    //pack data
    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            bufferSend[i+j*imax]=(*rho)(setIndex(sizeLocalGross,i+halo,j+halo,iref));
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)(*rho)(setIndex(sizeLocalGross,i+halo,j+halo,halo))+=bufferRec[i+j*imax];
    }

    free(bufferRec);
    free(bufferSend);

}

void vertexProjectionCIC_comm(Field<Real> * vel)
{

      if(vel->lattice().halo() == 0)
      {
          cout<< "NBCC_CIC_vel_comm: the field has to have at least a halo of 1" <<endl;
          cout<< "NBCC_CIC_vel_comm: aborting" <<endl;
          exit(-1);
      }

      Real *bufferSend;
      Real *bufferRec;


      long bufferSizeY;
      long bufferSizeZ;


      int sizeLocal[3];
      long sizeLocalGross[3];
      int sizeLocalOne[3];
      int halo = vel->lattice().halo();

      for(int i=0;i<3;i++)
      {
          sizeLocal[i]=vel->lattice().sizeLocal(i);
          sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
          sizeLocalOne[i]=sizeLocal[i]+2;
      }

      int distHaloOne = halo - 1;

      int iref;
      int imax;

      bufferSizeY =  (long)(sizeLocalOne[2]-1) * (long)sizeLocal[0];
      bufferSizeZ = (long)sizeLocal[0] * (long)sizeLocal[1] ;
      if(bufferSizeY>bufferSizeZ)
      {
          bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY*3);
          bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY*3);
      }
      else
      {
          bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ*3);
          bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ*3);
      }

      //send halo in direction X
      iref = sizeLocalGross[0] - halo;
      for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
      {
          for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
          {
              for(int c=0;c<3;c++)(*vel)(setIndex(sizeLocalGross,halo,j,k),c) += (*vel)(setIndex(sizeLocalGross,iref,j,k),c);
          }
      }
      //send halo in direction Y
      imax = sizeLocal[0];
      iref = sizeLocalGross[1]- halo;
      for(int k=0;k<(sizeLocalOne[2]-1);k++)
      {
          for(int i=0;i<imax;i++)
          {
              for(int c=0;c<3;c++)bufferSend[c+3*(i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
          }
      }
      parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY*3);
      for(int k=0;k<(sizeLocalOne[2]-1);k++)
      {
          for(int i=0;i<imax;i++)
          {
            for(int c=0;c<3;c++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+3*(i+k*imax)];
          }

      }

      //send halo in direction Z
      iref=sizeLocalGross[2]-halo;
      for(int j=0;j<(sizeLocalOne[1]-2);j++)
      {
          for(int i=0;i<imax;i++)
          {
              for(int c=0;c<3;c++)bufferSend[c+3*(i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
          }
      }
      parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ*3);
      for(int j=0;j<(sizeLocalOne[1]-2);j++)
      {
          for(int i=0;i<imax;i++)
          {
            for(int c=0;c<3;c++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+3*(i+j*imax)];
          }
      }

      free(bufferRec);
      free(bufferSend);
}

//////vector projection



/*! \fn template<typename part, typename part_info, typename part_dataType> void vectorProjectionCICNGP_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * vel, size_t * oset = NULL,int flag_where = FROM_INFO)

 \brief hybrid cloud-in-cell nearest-grid-point vector projection.

 Perform the projection of a vector defined as \f$\vec{p}= \rho \vec{v}\f$ where \f$\vec{v}\f$ is a vectorial (individual) property of the particles and \f$\rho\f$ a scalar (global or individual) one.



 \param Particles<part,part_info,part_dataType> * parts: bla
 \param Field<Real> * vel: pointer to the real vector field (a field with 3 components) on which the property is projected.
 \param size_t * oset: pointer to array of two elements. First element is the offset of the \f$\rho\f$ property, and second of the \f$\vec{v}\f$ property. If set to NULL, the projection will project the properties "mass" and "vel" of the particles.
 \param int flag_where: Flag to set if the property \f$\rho\f$ is a global or a individual one, i.e. the specify in which structure the property is. FROM_INFO set the property to be global and FROM_PART set the property to be individual. The property \f$\vec{v}\f$ is always individual

 \sa vectorProjectionCICNGP_comm(Field<Real> * rho)
 \sa projection_init(Field<Real> * f)
 */
template<typename part, typename part_info, typename part_dataType>
void vectorProjectionCICNGP_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * vel, size_t * oset = NULL,int flag_where = FROM_INFO)
{
    Site xPart(parts->lattice());
    Site xVel(vel->lattice());

    typename std::forward_list<part>::iterator it;

    double vi[36];//3 * 4 v0:0..3 v1:4..7 v2:8..11

    double mass;
    double latresolution = parts->res();

    double cicVol = latresolution * latresolution * latresolution;

    double weightScalarGridDown[3];
    double weightScalarGridUp[3];

    double referPos[3];


    int sfrom;

    size_t offset_mass;
    size_t offset_vel;



    if(oset == NULL)
    {
        sfrom = parts->mass_type();
        offset_mass = parts->mass_offset();
        offset_vel = offsetof(part,vel);
    }
    else
    {
        sfrom =  flag_where;
        offset_mass = oset[0];
        offset_vel = oset[1];
    }

    if(sfrom == FROM_INFO)
    {
        mass = *(double*)((char*)parts->parts_info() + offset_mass);
        mass /= cicVol;
        //cout << mass<<endl;
    }



    for(xPart.first(),xVel.first();xPart.test();xPart.next(),xVel.next())
    {
        if(parts->field()(xPart).size!=0)
        {
            for(int i=0;i<3;i++)

                referPos[i] = xPart.coord(i)*latresolution;

            for(int i=0;i<12;i++)vi[i]=0.0;

            for (it=(parts->field())(xPart).parts.begin(); it != (parts->field())(xPart).parts.end(); ++it)
            {
                for(int i =0;i<3;i++)
                {
                    weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / latresolution;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
                }


                double massVel;
                Real * v;

                if(sfrom==FROM_PART)
                {
                    mass = *(double*)((char*)&(*it)+offset_mass);
                    mass /= cicVol;
                }
                v = (Real*)((char*)&(*it)+offset_vel);
                //cout<< v[0]<<" , "<<v[1]<<" , "<<v[2]<<endl;

                massVel = mass*v[0];

                vi[0] +=  massVel * weightScalarGridDown[1] * weightScalarGridDown[2] ;
                vi[1] +=  massVel * weightScalarGridUp[1]   * weightScalarGridDown[2] ;
                vi[2] +=  massVel * weightScalarGridDown[1] * weightScalarGridUp[2] ;
                vi[3] +=  massVel * weightScalarGridUp[1]   * weightScalarGridUp[2] ;

                //working on v1

                massVel = mass*v[1];

                vi[4] +=  massVel * weightScalarGridDown[0] * weightScalarGridDown[2] ;
                vi[5] +=  massVel * weightScalarGridUp[0]   * weightScalarGridDown[2] ;
                vi[6] +=  massVel * weightScalarGridDown[0] * weightScalarGridUp[2] ;
                vi[7] +=  massVel * weightScalarGridUp[0]   * weightScalarGridUp[2] ;

                //working on v2

                massVel = mass*v[2];

                vi[8] +=  massVel * weightScalarGridDown[0] * weightScalarGridDown[1] ;
                vi[9] +=  massVel * weightScalarGridUp[0]   * weightScalarGridDown[1] ;
                vi[10]+=  massVel * weightScalarGridDown[0] * weightScalarGridUp[1] ;
                vi[11]+=  massVel * weightScalarGridUp[0]   * weightScalarGridUp[1] ;

            }


            (*vel)(xVel,0)+=vi[0];
            (*vel)(xVel,1)+=vi[4];
            (*vel)(xVel,2)+=vi[8];

            (*vel)(xVel+0,1)+=vi[5];
            (*vel)(xVel+0,2)+=vi[9];

            (*vel)(xVel+1,0)+=vi[1];
            (*vel)(xVel+1,2)+=vi[10];

            (*vel)(xVel+2,0)+=vi[2];
            (*vel)(xVel+2,1)+=vi[6];

            (*vel)(xVel+1+2,0)+=vi[3];
            (*vel)(xVel+0+2,1)+=vi[7];
            (*vel)(xVel+0+1,2)+=vi[11];

        }
    }

}


/*! \fn vectorProjectionCICNGP_comm(Field<Real> * vel)
 \brief communication method associated to the hybrid cloud-in-cell nearest-grid-point projection of vectors .


 \sa vectorProjectionCICNGP_project(...)
 */
void vectorProjectionCICNGP_comm(Field<Real> * vel)
{

    if(vel->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = vel->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=vel->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;




    int comp=3;


    iref =sizeLocalGross[0] - halo ;
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,halo,j,k),c) += (*vel)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }


    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }


    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);


    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

    //send halo in direction Z

    //pack data

    //cout<<"okok"<<endl;

    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
        }
    }

    free(bufferRec);
    free(bufferSend);


}







//////tensor projection


/*! \fn template<typename part, typename part_info, typename part_dataType> void symtensorProjectionCICNGP_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * Tij, size_t * oset = NULL,int flag_where = FROM_INFO)

 \brief hybrid cloud-in-cell nearest-grid-point symetric tensor projection.

 Perform the projection of a symmetric tensor defined as \f$\tau_{ij}= \rho v_i v_j\f$ where \f$\vec{v}\f$ is a vectorial (individual) property of the particles and \f$\rho\f$ a scalar (global or individual) one.



 \param Particles<part,part_info,part_dataType> * parts: bla
 \param Field<Real> * Tij: pointer to the real symmetric tensor field (a field in 3x3 components, and definded a symmetric) on which the property is projected.
 \param size_t * oset: pointer to array of two elements. First element is the offset of the \f$\rho\f$ property, and second of the \f$\vec{v}\f$ property. If set to NULL, the projection will project the properties "mass" and "vel" of the particles.
 \param int flag_where: Flag to set if the property \f$\rho\f$ is a global or a individual one, i.e. the specify in which structure the property is. FROM_INFO set the property to be global and FROM_PART set the property to be individual. The property \f$\vec{v}\f$ is always individual

 \sa symtensorProjectionCICNGP_comm(Field<Real> * rho)
 \sa projection_init(Field<Real> * f)
 */
template<typename part, typename part_info, typename part_dataType>
void symtensorProjectionCICNGP_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * Tij, size_t * oset = NULL,int flag_where = FROM_INFO)
{

    Site xPart(parts->lattice() );
    Site xTij(Tij->lattice() );

    typename std::forward_list<part>::iterator it;

    double mass;
    double latresolution = parts->res();
    double cicVol = latresolution * latresolution * latresolution;

    double  tij[6];
    double  tii[24];

    double weightScalarGridDown[3];
    double weightScalarGridUp[3];

    double referPos[3];


    int sfrom;

    size_t offset_mass;
    size_t offset_vel;



    if(oset == NULL)
    {
        sfrom = parts->mass_type();
        offset_mass = parts->mass_offset();
        offset_vel = offsetof(part,vel);
    }
    else
    {
        sfrom =  flag_where;
        offset_mass = oset[0];
        offset_vel = oset[1];
    }

    if(sfrom == FROM_INFO)
    {
        mass = *(double*)((char*)parts->parts_info() + offset_mass);
        mass /= cicVol;
        //cout << mass<<endl;
    }



    for(xPart.first(),xTij.first();xPart.test();xPart.next(),xTij.next())
    {
        if(parts->field()(xPart).size!=0)
        {
            for(int i=0;i<3;i++)
                referPos[i] = (double)xPart.coord(i)*latresolution;

            for(int i=0;i<6;i++)tij[i]=0.0;
            for(int i=0;i<24;i++)tii[i]=0.0;

            for (it=(parts->field())(xPart).parts.begin(); it != (parts->field())(xPart).parts.end(); ++it)
            {
                for(int i =0;i<3;i++)
                {
                    weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / latresolution;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];
                }

                Real * vel;
                if(sfrom==FROM_PART)
                {
                    mass = *(double*)((char*)&(*it)+offset_mass);
                    mass /= cicVol;
                }
                vel = (Real*)((char*)&(*it)+offset_vel);


                //working on scalars mass * v_i * v_i

                for(int i=0;i<3;i++)
                {
                    double massVel2 = mass * vel[i] * vel[i];
                    //000
                    tii[0+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridDown[2];
                    //001
                    tii[1+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridUp[2];
                    //010
                    tii[2+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridDown[2];
                    //011
                    tii[3+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridUp[2];
                    //100
                    tii[4+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridDown[2];
                    //101
                    tii[5+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridUp[2];
                    //110
                    tii[6+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridDown[2];
                    //111
                    tii[7+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridUp[2];
                }

                double massVelVel;

                massVelVel = mass * vel[0] * vel[1];
                tij[0] +=  massVelVel * weightScalarGridDown[2];
                tij[1] +=  massVelVel * weightScalarGridUp[2];

                massVelVel = mass * vel[0] * vel[2];
                tij[2] +=  massVelVel * weightScalarGridDown[1];
                tij[3] +=  massVelVel * weightScalarGridUp[1];

                massVelVel = mass * vel[1] * vel[2];
                tij[4] +=  massVelVel * weightScalarGridDown[0];
                tij[5] +=  massVelVel * weightScalarGridUp[0];

            }


            for(int i=0;i<3;i++)(*Tij)(xTij,i,i)+=tii[8*i];
            (*Tij)(xTij,0,1)+=tij[0];
            (*Tij)(xTij,0,2)+=tij[2];
            (*Tij)(xTij,1,2)+=tij[4];

            for(int i=0;i<3;i++)(*Tij)(xTij+0,i,i)+=tii[4+8*i];
            (*Tij)(xTij+0,1,2)+=tij[5];

            for(int i=0;i<3;i++)(*Tij)(xTij+1,i,i)+=tii[2+8*i];
            (*Tij)(xTij+1,0,2)+=tij[3];

            for(int i=0;i<3;i++)(*Tij)(xTij+2,i,i)+=tii[1+8*i];
            (*Tij)(xTij+2,0,1)+=tij[1];

            for(int i=0;i<3;i++)(*Tij)(xTij+0+1,i,i)+=tii[6+8*i];
            for(int i=0;i<3;i++)(*Tij)(xTij+0+2,i,i)+=tii[5+8*i];
            for(int i=0;i<3;i++)(*Tij)(xTij+1+2,i,i)+=tii[3+8*i];
            for(int i=0;i<3;i++)(*Tij)(xTij+0+1+2,i,i)+=tii[7+8*i];
        }
    }


}

/*! \fn symtensorProjectionCICNGP_comm(Field<Real> * vel)
 \brief communication method associated to the hybrid cloud-in-cell nearest-grid-point projection of symmetric tensor .


 \sa symtensorProjectionCICNGP_project(...)
 */
void symtensorProjectionCICNGP_comm(Field<Real> * Tij)
{

    if(Tij->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSend;
    Real *bufferRec;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = Tij->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=Tij->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int iref;
    int imax;




    int comp=6;

    iref = sizeLocalGross[0]-halo;
    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,halo,j,k),c) += (*Tij)(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }


    //send halo in direction Y
    bufferSizeY =  (long)(sizeLocalOne[2]-1)*sizeLocal[0] * comp;
    bufferSizeZ = sizeLocal[0] * sizeLocal[1] * comp;
    if(bufferSizeY>bufferSizeZ)
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeY);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeY);
    }
    else
    {
        bufferSend = (Real*)malloc(sizeof(Real)*bufferSizeZ);
        bufferRec = (Real*)malloc(sizeof(Real)*bufferSizeZ);
    }


    //pack data
    imax=sizeLocalGross[0]-2* halo;
    iref=sizeLocalGross[1]- halo;
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
        }
    }


    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);


    //unpack data
    for(int k=0;k<(sizeLocalOne[2]-1);k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[c+comp*(i+k*imax)];
        }

    }

    //send halo in direction Z

    //pack data
    iref=sizeLocalGross[2]-halo;
    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)bufferSend[c+comp*(i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
        }
    }

    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int c=0;c<comp;c++)(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[c+comp*(i+j*imax)];
        }
    }

    free(bufferRec);
    free(bufferSend);



}

#ifndef DOXYGEN_SHOULD_SKIP_THIS
template<typename part, typename part_info, typename part_dataType>
void vectorProjectionCIC_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * vel, size_t * oset = NULL,int flag_where = FROM_INFO)
{

    Site xPart(parts->lattice());
    Site xVel(vel->lattice());

    typename std::forward_list<part>::iterator it;

    double mass;
    double latresolution = parts->res();
    double cicVol = latresolution * latresolution * latresolution;

    double vi[36]; // 3 * 12 v0:0..11 v1:12..23 v2:24..35

    int gridChoice[3];

    double weightScalarGridDown[3];
    double weightScalarGridUp[3];
    double weightTensorGridDown[3];
    double weightTensorGridUp[3];

    double referDistShift;

    double referPos[3];
    double referPosShift[3];



    int sfrom;

    size_t offset_mass;
    size_t offset_vel;



    if(oset == NULL)
    {
        sfrom = parts->mass_type();
        offset_mass = parts->mass_offset();
        offset_vel = offsetof(part,vel);
    }
    else
    {
        sfrom =  flag_where;
        offset_mass = oset[0];
        offset_vel = oset[1];
    }

    if(sfrom == FROM_INFO)
    {
        mass = *(double*)((char*)parts->parts_info() + offset_mass);
        mass /= cicVol;
        //cout << mass<<endl;
    }


    for(xPart.first(),xVel.first();xPart.test();xPart.next(),xVel.next())
    {
        if(parts->field()(xPart).size!=0)
        {
            for(int i=0;i<3;i++)
            {
                referPos[i] = xPart.coord(i)*latresolution;
                referPosShift[i] = referPos[i] +  latresolution* 0.5 ;
            }

            for(int i=0;i<36;i++)vi[i]=0.0;


            for (it=(parts->field())(xPart).parts.begin(); it != (parts->field())(xPart).parts.end(); ++it)
            {

                for(int i =0;i<3;i++)
                {
                    weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / latresolution;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];

                    referDistShift = ((*it).pos[i] - referPosShift[i]) / latresolution;
                    if(referDistShift<0)
                    {
                        gridChoice[i]=0;
                        weightTensorGridDown[i] = -referDistShift;
                        weightTensorGridUp[i] = 1.0l - weightTensorGridDown[i];
                    }
                    else
                    {
                        gridChoice[i]=1;
                        weightTensorGridUp[i] = referDistShift;
                        weightTensorGridDown[i] = 1.0l - weightTensorGridUp[i];
                    }

                }

                Real * v;
                double massVel;


                if(sfrom==FROM_PART)
                {
                    mass = *(double*)((char*)&(*it)+offset_mass);
                    mass /= cicVol;
                }

                v = (Real*)((char*)&(*it)+offset_vel);

                //do the projection into the buffer

                //working on v0
                massVel = mass*v[0];
                vi[ CIC_mapv0[gridChoice[0]][0][0][0] ] +=  massVel * weightTensorGridDown[0] * weightScalarGridDown[1] * weightScalarGridDown[2];
                vi[ CIC_mapv0[gridChoice[0]][0][0][1] ] +=  massVel * weightTensorGridDown[0] * weightScalarGridDown[1] * weightScalarGridUp[2];
                vi[ CIC_mapv0[gridChoice[0]][0][1][0] ] +=  massVel * weightTensorGridDown[0] * weightScalarGridUp[1]   * weightScalarGridDown[2];
                vi[ CIC_mapv0[gridChoice[0]][0][1][1] ] +=  massVel * weightTensorGridDown[0] * weightScalarGridUp[1]   * weightScalarGridUp[2];
                vi[ CIC_mapv0[gridChoice[0]][1][0][0] ] +=  massVel * weightTensorGridUp[0]   * weightScalarGridDown[1] * weightScalarGridDown[2];
                vi[ CIC_mapv0[gridChoice[0]][1][0][1] ] +=  massVel * weightTensorGridUp[0]   * weightScalarGridDown[1] * weightScalarGridUp[2];
                vi[ CIC_mapv0[gridChoice[0]][1][1][0] ] +=  massVel * weightTensorGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridDown[2];
                vi[ CIC_mapv0[gridChoice[0]][1][1][1] ] +=  massVel * weightTensorGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridUp[2];

                //working on v1
                massVel = mass*v[1];
                vi[ CIC_mapv1[gridChoice[1]][0][0][0] ] +=  massVel * weightScalarGridDown[0] * weightTensorGridDown[1] * weightScalarGridDown[2];
                vi[ CIC_mapv1[gridChoice[1]][0][0][1] ] +=  massVel * weightScalarGridDown[0] * weightTensorGridDown[1] * weightScalarGridUp[2];
                vi[ CIC_mapv1[gridChoice[1]][0][1][0] ] +=  massVel * weightScalarGridDown[0] * weightTensorGridUp[1]   * weightScalarGridDown[2];
                vi[ CIC_mapv1[gridChoice[1]][0][1][1] ] +=  massVel * weightScalarGridDown[0] * weightTensorGridUp[1]   * weightScalarGridUp[2];
                vi[ CIC_mapv1[gridChoice[1]][1][0][0] ] +=  massVel * weightScalarGridUp[0]   * weightTensorGridDown[1] * weightScalarGridDown[2];
                vi[ CIC_mapv1[gridChoice[1]][1][0][1] ] +=  massVel * weightScalarGridUp[0]   * weightTensorGridDown[1] * weightScalarGridUp[2];
                vi[ CIC_mapv1[gridChoice[1]][1][1][0] ] +=  massVel * weightScalarGridUp[0]   * weightTensorGridUp[1]   * weightScalarGridDown[2];
                vi[ CIC_mapv1[gridChoice[1]][1][1][1] ] +=  massVel * weightScalarGridUp[0]   * weightTensorGridUp[1]   * weightScalarGridUp[2];


                //working on v2
                massVel = mass*v[2];
                vi[ CIC_mapv2[gridChoice[2]][0][0][0] ] +=  massVel * weightScalarGridDown[0] * weightScalarGridDown[1] * weightTensorGridDown[2];
                vi[ CIC_mapv2[gridChoice[2]][0][0][1] ] +=  massVel * weightScalarGridDown[0] * weightScalarGridDown[1] * weightTensorGridUp[2];
                vi[ CIC_mapv2[gridChoice[2]][0][1][0] ] +=  massVel * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightTensorGridDown[2];
                vi[ CIC_mapv2[gridChoice[2]][0][1][1] ] +=  massVel * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightTensorGridUp[2];
                vi[ CIC_mapv2[gridChoice[2]][1][0][0] ] +=  massVel * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightTensorGridDown[2];
                vi[ CIC_mapv2[gridChoice[2]][1][0][1] ] +=  massVel * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightTensorGridUp[2];
                vi[ CIC_mapv2[gridChoice[2]][1][1][0] ] +=  massVel * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightTensorGridDown[2];
                vi[ CIC_mapv2[gridChoice[2]][1][1][1] ] +=  massVel * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightTensorGridUp[2];

            }

            //copy to field

            (*vel)(xVel -0       ,0) += vi[0];//v0

            (*vel)(xVel -0    +2 ,0) += vi[1];//v0

            (*vel)(xVel -0 +1    ,0) += vi[2];//v0

            (*vel)(xVel -0 +1 +2 ,0) += vi[3];//v0

            //--------------------------------------------

            (*vel)(xVel    -1    ,1) += vi[12];//v1

            (*vel)(xVel    -1 +2 ,1) += vi[14];//v1

            (*vel)(xVel       -2 ,2) += vi[24];//v2

            (*vel)(xVel          ,0) += vi[4];//v0
            (*vel)(xVel          ,1) += vi[16];//v1
            (*vel)(xVel          ,2) += vi[28];//v2

            (*vel)(xVel       +2 ,0) += vi[5];//v0
            (*vel)(xVel       +2 ,1) += vi[18];//v1
            (*vel)(xVel       +2 ,2) += vi[32];//v2

            (*vel)(xVel    +1 -2 ,2) += vi[25];//v2

            (*vel)(xVel    +1    ,0) += vi[6];//v0
            (*vel)(xVel    +1    ,1) += vi[20];//v1
            (*vel)(xVel    +1    ,2) += vi[29];//v2

            (*vel)(xVel    +1 +2 ,0) += vi[7];//v0
            (*vel)(xVel    +1 +2 ,1) += vi[22];//v1
            (*vel)(xVel    +1 +2 ,2) += vi[33];//v2

            //--------------------------------------------

            (*vel)(xVel +0 -1    ,1) += vi[13];//v1

            (*vel)(xVel +0 -1 +2 ,1) += vi[15];//v1

            (*vel)(xVel +0    -2 ,2) += vi[26];//v2

            (*vel)(xVel +0       ,0) += vi[8];//v0
            (*vel)(xVel +0       ,1) += vi[17];//v1
            (*vel)(xVel +0       ,2) += vi[30];//v2

            (*vel)(xVel +0    +2 ,0) += vi[9];//v0
            (*vel)(xVel +0    +2 ,1) += vi[19];//v1
            (*vel)(xVel +0    +2 ,2) += vi[34];//v2

            (*vel)(xVel +0 +1 -2 ,2) += vi[27];//v2

            (*vel)(xVel +0 +1    ,0) += vi[10];//v0
            (*vel)(xVel +0 +1    ,1) += vi[21];//v1
            (*vel)(xVel +0 +1    ,2) += vi[31];//v2

            (*vel)(xVel +0 +1 +2 ,0) += vi[11];//v0
            (*vel)(xVel +0 +1 +2 ,1) += vi[23];//v1
            (*vel)(xVel +0 +1 +2 ,2) += vi[35];//v2

        }
    }


}


void vectorProjectionCIC_comm(Field<Real> * vel)
{


    if(vel->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSendUp;
    Real *bufferSendDown;

    Real *bufferRecUp;
    Real *bufferRecDown;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = vel->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=vel->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;

    int imax;
    //update from halo data in X;

    int iref1=sizeLocalGross[0]-halo;
    int iref2=sizeLocalGross[0]-halo-1;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int comp=0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,halo,j,k),comp) += (*vel)(setIndex(sizeLocalGross,iref1,j,k),comp);

            for(int comp=0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,iref2,j,k),comp) += (*vel)(setIndex(sizeLocalGross,distHaloOne,j,k),comp);
        }
    }

    //communication



    //build buffer size

    bufferSizeY = (long)(sizeLocalOne[2]) * sizeLocal[0];
    bufferSizeZ = sizeLocal[0] * sizeLocal[1];

    if(bufferSizeY>bufferSizeZ)
    {
        bufferSendUp = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferRecUp = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
    }
    else
    {
        bufferSendUp = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferRecUp = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
    }


    //working on Y dimension


    imax=sizeLocal[0];
    iref1 = sizeLocalGross[1]- halo;


    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendUp[comp + 3l * (i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),comp);
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendDown[comp + 3l * (i+k*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),comp);
        }
    }

    parallel.sendUpDown_dim1(bufferSendUp,bufferRecUp,bufferSizeY * 3l,bufferSendDown,bufferRecDown,bufferSizeY * 3l);

    //unpack data
    iref1 = sizeLocalGross[1]- halo - 1;

    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,halo,k+distHaloOne),comp) += bufferRecUp[comp + 3l * (i+k*imax)];
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),comp) += bufferRecDown[comp + 3l * (i+k*imax)];

        }
    }

    //workin on dim z

    //pack data

    iref1=sizeLocalGross[2]- halo;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendUp[comp + 3l * (i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),comp);
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)bufferSendDown[comp + 3l * (i+j*imax)]=(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),comp);
        }

    }

    parallel.sendUpDown_dim0(bufferSendUp,bufferRecUp,bufferSizeY * 3l,bufferSendDown,bufferRecDown,bufferSizeY * 3l);

    //unpack data

    iref1 = sizeLocalGross[2]- halo - 1;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,halo),comp) += bufferRecUp[comp + 3l * (i+j*imax)];
        }
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<3;comp++)(*vel)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),comp) += bufferRecDown[comp + 3l * (i+j*imax)];
        }

    }


    free(bufferSendUp);
    free(bufferSendDown);
    free(bufferRecUp);
    free(bufferRecDown);


}



template<typename part, typename part_info, typename part_dataType>
void VecVecProjectionCIC_project(Particles<part,part_info,part_dataType> * parts,Field<Real> * Tij, size_t * oset = NULL,int flag_where = FROM_INFO)
{
    Site xPart(parts->lattice());
    Site xTij(Tij->lattice());
    typename std::forward_list<part>::iterator it;


    double mass;
    double latresolution = parts->res();
    double cicVol = latresolution * latresolution * latresolution;


    double  tij[54];
    double  tii[24];

    int gridChoice[3];

    double weightScalarGridDown[3];
    double weightScalarGridUp[3];
    double weightTensorGridDown[3];
    double weightTensorGridUp[3];

    double referDistShift;

    double referPos[3];
    double referPosShift[3];

    int sfrom;

    size_t offset_mass;
    size_t offset_vel;



    if(oset == NULL)
    {
        sfrom = parts->mass_type();
        offset_mass = parts->mass_offset();
        offset_vel = offsetof(part,vel);
    }
    else
    {
        sfrom =  flag_where;
        offset_mass = oset[0];
        offset_vel = oset[1];
    }

    if(sfrom == FROM_INFO)
    {
        mass = *(double*)((char*)parts->parts_info() + offset_mass);
        mass /= cicVol;
        //cout << mass<<endl;
    }



    for(xPart.first(),xTij.first();xPart.test();xPart.next(),xTij.next())
    {
        if(parts->field()(xPart).size!=0)
        {
            for(int i=0;i<3;i++)
            {
                referPos[i] = (double)xPart.coord(i)*latresolution;
                referPosShift[i] = referPos[i] +  latresolution* 0.5 ;
            }

            for(int i=0;i<54;i++)tij[i]=0.0;
            for(int i=0;i<24;i++)tii[i]=0.0;

            for (it=(parts->field())(xPart).parts.begin(); it != (parts->field())(xPart).parts.end(); ++it)
            {
                for(int i =0;i<3;i++)
                {
                    weightScalarGridUp[i] = ((*it).pos[i] - referPos[i]) / latresolution;
                    weightScalarGridDown[i] = 1.0l - weightScalarGridUp[i];

                    referDistShift = ((*it).pos[i] - referPosShift[i]) / latresolution;
                    if(referDistShift<0)
                    {
                        gridChoice[i]=0;
                        weightTensorGridDown[i] = -referDistShift;
                        weightTensorGridUp[i] = 1.0l - weightTensorGridDown[i];
                    }
                    else
                    {
                        gridChoice[i]=1;
                        weightTensorGridUp[i] = referDistShift;
                        weightTensorGridDown[i] = 1.0l - weightTensorGridUp[i];
                    }

                }
                Real * vel;
                if(sfrom==FROM_PART)
                {
                    mass = *(double*)((char*)&(*it)+offset_mass);
                    mass /= cicVol;
                }
                vel = (Real*)((char*)&(*it)+offset_vel);


                //working on scalars mass * v_i * v_i

                for(int i=0;i<3;i++)
                {
                    double massVel2 = mass * vel[i] * vel[i];
                    //000
                    tii[0+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridDown[2];
                    //001
                    tii[1+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridDown[1] * weightScalarGridUp[2];
                    //010
                    tii[2+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridDown[2];
                    //011
                    tii[3+i*8] += massVel2 * weightScalarGridDown[0] * weightScalarGridUp[1]   * weightScalarGridUp[2];
                    //100
                    tii[4+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridDown[2];
                    //101
                    tii[5+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridDown[1] * weightScalarGridUp[2];
                    //110
                    tii[6+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridDown[2];
                    //111
                    tii[7+i*8] += massVel2 * weightScalarGridUp[0]   * weightScalarGridUp[1]   * weightScalarGridUp[2];
                }

                //working on plaquette mass * v_i * v_j
                int whichCube;
                double massVelVel;

                //1) m * v_1 * v_2
                whichCube  = gridChoice[1] + 2 * gridChoice[2];
                massVelVel = mass * vel[1] * vel[2];

                tij[ mapv1v2[whichCube][0][0][0] ] +=  massVelVel * weightScalarGridDown[0] * weightTensorGridDown[1] * weightTensorGridDown[2];
                tij[ mapv1v2[whichCube][0][0][1] ] +=  massVelVel * weightScalarGridDown[0] * weightTensorGridDown[1] * weightTensorGridUp[2];
                tij[ mapv1v2[whichCube][0][1][0] ] +=  massVelVel * weightScalarGridDown[0] * weightTensorGridUp[1]   * weightTensorGridDown[2];
                tij[ mapv1v2[whichCube][0][1][1] ] +=  massVelVel * weightScalarGridDown[0] * weightTensorGridUp[1]   * weightTensorGridUp[2];
                tij[ mapv1v2[whichCube][1][0][0] ] +=  massVelVel * weightScalarGridUp[0]   * weightTensorGridDown[1] * weightTensorGridDown[2];
                tij[ mapv1v2[whichCube][1][0][1] ] +=  massVelVel * weightScalarGridUp[0]   * weightTensorGridDown[1] * weightTensorGridUp[2];
                tij[ mapv1v2[whichCube][1][1][0] ] +=  massVelVel * weightScalarGridUp[0]   * weightTensorGridUp[1]   * weightTensorGridDown[2];
                tij[ mapv1v2[whichCube][1][1][1] ] +=  massVelVel * weightScalarGridUp[0]   * weightTensorGridUp[1]   * weightTensorGridUp[2];

                //2) m * v_0 * v_1
                whichCube  = gridChoice[0] + 2 * gridChoice[1];
                massVelVel = mass * vel[0] * vel[1];

                tij[ mapv0v1[whichCube][0][0][0] ] +=  massVelVel * weightTensorGridDown[0] * weightTensorGridDown[1] * weightScalarGridDown[2];
                tij[ mapv0v1[whichCube][0][0][1] ] +=  massVelVel * weightTensorGridDown[0] * weightTensorGridDown[1] * weightScalarGridUp[2];
                tij[ mapv0v1[whichCube][0][1][0] ] +=  massVelVel * weightTensorGridDown[0] * weightTensorGridUp[1]   * weightScalarGridDown[2];
                tij[ mapv0v1[whichCube][0][1][1] ] +=  massVelVel * weightTensorGridDown[0] * weightTensorGridUp[1]   * weightScalarGridUp[2];
                tij[ mapv0v1[whichCube][1][0][0] ] +=  massVelVel * weightTensorGridUp[0]   * weightTensorGridDown[1] * weightScalarGridDown[2];
                tij[ mapv0v1[whichCube][1][0][1] ] +=  massVelVel * weightTensorGridUp[0]   * weightTensorGridDown[1] * weightScalarGridUp[2];
                tij[ mapv0v1[whichCube][1][1][0] ] +=  massVelVel * weightTensorGridUp[0]   * weightTensorGridUp[1]   * weightScalarGridDown[2];
                tij[ mapv0v1[whichCube][1][1][1] ] +=  massVelVel * weightTensorGridUp[0]   * weightTensorGridUp[1]   * weightScalarGridUp[2];


                //3) m * v_0 * v_2
                whichCube = gridChoice[0] + 2 * gridChoice[2];
                massVelVel = mass * vel[0] * vel[2];

                tij[ mapv0v2[whichCube][0][0][0] ] +=  massVelVel * weightTensorGridDown[0] * weightScalarGridDown[1] * weightTensorGridDown[2];
                tij[ mapv0v2[whichCube][0][0][1] ] +=  massVelVel * weightTensorGridDown[0] * weightScalarGridDown[1] * weightTensorGridUp[2];
                tij[ mapv0v2[whichCube][0][1][0] ] +=  massVelVel * weightTensorGridDown[0] * weightScalarGridUp[1]   * weightTensorGridDown[2];
                tij[ mapv0v2[whichCube][0][1][1] ] +=  massVelVel * weightTensorGridDown[0] * weightScalarGridUp[1]   * weightTensorGridUp[2];
                tij[ mapv0v2[whichCube][1][0][0] ] +=  massVelVel * weightTensorGridUp[0]   * weightScalarGridDown[1] * weightTensorGridDown[2];
                tij[ mapv0v2[whichCube][1][0][1] ] +=  massVelVel * weightTensorGridUp[0]   * weightScalarGridDown[1] * weightTensorGridUp[2];
                tij[ mapv0v2[whichCube][1][1][0] ] +=  massVelVel * weightTensorGridUp[0]   * weightScalarGridUp[1]   * weightTensorGridDown[2];
                tij[ mapv0v2[whichCube][1][1][1] ] +=  massVelVel * weightTensorGridUp[0]   * weightScalarGridUp[1]   * weightTensorGridUp[2];

            }

            //add to the field T_ij

            (*Tij)(xTij -0 -1    ,0,1) += tij[36];

            (*Tij)(xTij -0 -1 +2 ,0,1) += tij[37];

            (*Tij)(xTij -0    -2 ,0,2) += tij[18];

            (*Tij)(xTij -0       ,0,1) += tij[42];
            (*Tij)(xTij -0       ,0,2) += tij[24];

            (*Tij)(xTij -0    +2 ,0,1) += tij[43];
            (*Tij)(xTij -0    +2 ,0,2) += tij[30];

            (*Tij)(xTij -0 +1 -2 ,0,2) += tij[19];

            (*Tij)(xTij -0 +1    ,0,1) += tij[48];
            (*Tij)(xTij -0 +1    ,0,2) += tij[25];

            (*Tij)(xTij -0 +1 +2 ,0,1) += tij[49];
            (*Tij)(xTij -0 +1 +2 ,0,2) += tij[31];

            //--------------------------------------------

            (*Tij)(xTij    -1 -2 ,1,2) += tij[0];

            (*Tij)(xTij    -1    ,0,1) += tij[38];
            (*Tij)(xTij    -1    ,1,2) += tij[6];

            (*Tij)(xTij    -1 +2 ,0,1) += tij[39];
            (*Tij)(xTij    -1 +2 ,1,2) += tij[12];

            (*Tij)(xTij       -2 ,0,2) += tij[20];
            (*Tij)(xTij       -2 ,1,2) += tij[2];

            (*Tij)(xTij          ,0,1) += tij[44];
            (*Tij)(xTij          ,0,2) += tij[26];
            (*Tij)(xTij          ,1,2) += tij[8];
            for(int i=0;i<3;i++)(*Tij)(xTij,i,i)+=tii[8*i];

            (*Tij)(xTij       +2 ,0,1) += tij[45];
            (*Tij)(xTij       +2 ,0,2) += tij[32];
            (*Tij)(xTij       +2 ,1,2) += tij[14];
            for(int i=0;i<3;i++)(*Tij)(xTij+2,i,i)+=tii[1+8*i];

            (*Tij)(xTij    +1 -2 ,0,2) += tij[21];
            (*Tij)(xTij    +1 -2 ,1,2) += tij[4];

            (*Tij)(xTij    +1    ,0,1) += tij[50];
            (*Tij)(xTij    +1    ,0,2) += tij[27];
            (*Tij)(xTij    +1    ,1,2) += tij[10];
            for(int i=0;i<3;i++)(*Tij)(xTij+1,i,i)+=tii[2+8*i];

            (*Tij)(xTij    +1 +2 ,0,1) += tij[51];
            (*Tij)(xTij    +1 +2 ,0,2) += tij[33];
            (*Tij)(xTij    +1 +2 ,1,2) += tij[16];
            for(int i=0;i<3;i++)(*Tij)(xTij+1+2,i,i)+=tii[3+8*i];

            //--------------------------------------------

            (*Tij)(xTij +0 -1 -2 ,1,2) += tij[1];

            (*Tij)(xTij +0 -1    ,0,1) += tij[40];
            (*Tij)(xTij +0 -1    ,1,2) += tij[7];

            (*Tij)(xTij +0 -1 +2 ,0,1) += tij[41];
            (*Tij)(xTij +0 -1 +2 ,1,2) += tij[13];

            (*Tij)(xTij +0    -2 ,0,2) += tij[22];
            (*Tij)(xTij +0    -2 ,1,2) += tij[3];

            (*Tij)(xTij +0       ,0,1) += tij[46];
            (*Tij)(xTij +0       ,0,2) += tij[28];
            (*Tij)(xTij +0       ,1,2) += tij[9];
            for(int i=0;i<3;i++)(*Tij)(xTij+0,i,i)+=tii[4+8*i];

            (*Tij)(xTij +0    +2 ,0,1) += tij[47];
            (*Tij)(xTij +0    +2 ,0,2) += tij[34];
            (*Tij)(xTij +0    +2 ,1,2) += tij[15];
            for(int i=0;i<3;i++)(*Tij)(xTij+0+2,i,i)+=tii[5+8*i];

            (*Tij)(xTij +0 +1 -2 ,0,2) += tij[23];
            (*Tij)(xTij +0 +1 -2 ,1,2) += tij[5];

            (*Tij)(xTij +0 +1    ,0,1) += tij[52];
            (*Tij)(xTij +0 +1    ,0,2) += tij[29];
            (*Tij)(xTij +0 +1    ,1,2) += tij[11];
            for(int i=0;i<3;i++)(*Tij)(xTij+0+1,i,i)+=tii[6+8*i];

            (*Tij)(xTij +0 +1 +2 ,0,1) += tij[53];
            (*Tij)(xTij +0 +1 +2 ,0,2) += tij[35];
            (*Tij)(xTij +0 +1 +2 ,1,2) += tij[17];
            for(int i=0;i<3;i++)(*Tij)(xTij+0+1+2,i,i)+=tii[7+8*i];

        }
    }

}

void VecVecProjectionCIC_comm(Field<Real> * Tij)
{


    if(Tij->lattice().halo() == 0)
    {
        cout<< "LATfield2::scalarProjectionCIC_proj: the field has to have at least a halo of 1" <<endl;
        cout<< "LATfield2::scalarProjectionCIC_proj: aborting" <<endl;
        exit(-1);
    }

    Real *bufferSendUp;
    Real *bufferSendDown;

    Real *bufferRecUp;
    Real *bufferRecDown;


    long bufferSizeY;
    long bufferSizeZ;


    int sizeLocal[3];
    long sizeLocalGross[3];
    int sizeLocalOne[3];
    int halo = Tij->lattice().halo();

    for(int i=0;i<3;i++)
    {
        sizeLocal[i]=Tij->lattice().sizeLocal(i);
        sizeLocalGross[i] = sizeLocal[i] + 2 * halo;
        sizeLocalOne[i]=sizeLocal[i]+2;
    }

    int distHaloOne = halo - 1;








    int imax;

    int iref1=sizeLocalGross[0]-halo;
    int iref2=sizeLocalGross[0]-halo-1;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            for(int comp=0;comp<6;comp++)(*Tij)(setIndex(sizeLocalGross,halo,j,k),comp) += (*Tij)(setIndex(sizeLocalGross,iref1,j,k),comp);

            (*Tij)(setIndex(sizeLocalGross,iref2,j,k),0,1) += (*Tij)(setIndex(sizeLocalGross,distHaloOne,j,k),0,1);
            (*Tij)(setIndex(sizeLocalGross,iref2,j,k),0,2) += (*Tij)(setIndex(sizeLocalGross,distHaloOne,j,k),0,2);
            (*Tij)(setIndex(sizeLocalGross,iref2,j,k),1,2) += (*Tij)(setIndex(sizeLocalGross,distHaloOne,j,k),1,2);
        }
    }



    //build buffer size

    bufferSizeY = (long)(sizeLocalOne[2]) * (long)sizeLocal[0];
    bufferSizeZ = ((long)sizeLocal[0]) * ((long)sizeLocal[1]);



    if(bufferSizeY >= bufferSizeZ)
    {
        bufferSendUp = (Real*)malloc(6*bufferSizeY * sizeof(Real));
        bufferRecUp = (Real*)malloc(6*bufferSizeY * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeY * sizeof(Real));
    }
    else //if(bufferSizeZ> bufferSizeY )
    {
        bufferSendUp = (Real*)malloc(6*bufferSizeZ * sizeof(Real));
        bufferRecUp = (Real*)malloc(6*bufferSizeZ * sizeof(Real));
        bufferSendDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
        bufferRecDown = (Real*)malloc(3*bufferSizeZ * sizeof(Real));
    }


    //working on Y dimension;


    ///verif:



    imax=sizeLocal[0];
    iref1 = sizeLocalGross[1]- halo;

    //COUT<<(*Tij)(setIndex(sizeLocalGross,6+halo,iref1,6+halo),0)<<endl;

    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)bufferSendUp[comp + 6l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),comp);
        }

        for(int i=0;i<imax;i++)
        {
            bufferSendDown[    3l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),0,1);
            bufferSendDown[1l + 3l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),0,2);
            bufferSendDown[2l + 3l * (i+k*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,distHaloOne,k+distHaloOne),1,2);
        }

    }

    //COUT<<bufferSendUp[0l + 6l * (6+7*imax)]<<endl;

    parallel.sendUpDown_dim1(bufferSendUp,bufferRecUp,bufferSizeY * 6l,bufferSendDown,bufferRecDown,bufferSizeY * 3l);

    //if(parallel.rank()==2)cout<<bufferRecUp[0l + 6l * (6+7*imax)]<<endl;

    //unpack data
    iref1 = sizeLocalGross[1]- halo - 1;


    for(int k=0;k<sizeLocalOne[2];k++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)(*Tij)(setIndex(sizeLocalGross,i+halo,halo,k+distHaloOne),comp) += bufferRecUp[comp + 6l * (i+k*imax)];
        }

        for(int i=0;i<imax;i++)
        {
            (*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),0,1) += bufferRecDown[     3l * (i+k*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),0,2) += bufferRecDown[1l + 3l * (i+k*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,iref1,k+distHaloOne),1,2) += bufferRecDown[2l + 3l * (i+k*imax)];
        }

    }

    //if(parallel.rank()==2)cout<<(*Tij)(setIndex(sizeLocalGross,6+halo,halo ,6+halo),0)<<endl;

    //working on Z dimension;

    //pack data

    iref1=sizeLocalGross[2]- halo;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)bufferSendUp[comp + 6l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),comp);
        }
        for(int i=0;i<imax;i++)
        {
            bufferSendDown[     3l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),0,1);
            bufferSendDown[1l + 3l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),0,2);
            bufferSendDown[2l + 3l * (i+j*imax)]=(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,distHaloOne),1,2);
        }

    }

    parallel.sendUpDown_dim0(bufferSendUp,bufferRecUp,bufferSizeZ * 6l,bufferSendDown,bufferRecDown,bufferSizeZ * 3l);

    //unpack data

    iref1 = sizeLocalGross[2]- halo - 1;

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)
        {
            for(int comp =0;comp<6;comp++)(*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,halo),comp) += bufferRecUp[comp + 6l * (i+j*imax)];
        }

        for(int i=0;i<imax;i++)
        {

            (*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),0,1) += bufferRecDown[     3l * (i+j*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),0,2) += bufferRecDown[1l + 3l * (i+j*imax)];
            (*Tij)(setIndex(sizeLocalGross,i+halo,j+halo,iref1),1,2) += bufferRecDown[2l + 3l * (i+j*imax)];
        }


    }

    free(bufferSendUp);
    free(bufferSendDown);
    free(bufferRecUp);
    free(bufferRecDown);
}

#endif

/**@}*/
#endif
