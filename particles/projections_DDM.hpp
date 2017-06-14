#ifndef LATFIELD2_PROJECTIONS_DDM_HPP
#define LATFIELD2_PROJECTIONS_DDM_HPP


template<typename part, typename part_info, typename part_dataType>
void DDM_CIC_project(Particles<part,part_info,part_dataType> * parts,
                     Field<Real> * gij,
                     Field<Real> * E,
                     Field<Real> * Pi,
                     Field<Real> * Sij,
                     size_t * oset = NULL,int flag_where = FROM_INFO)
{
  if(E->lattice().halo() == 0 || Pi->lattice().halo() == 0 || Sij->lattice().halo() == 0)
  {
      cout<< "LATfield2::DDM_CIC_project: the field has to have at least a halo of 1" <<endl;
      cout<< "LATfield2::DDM_CIC_project: aborting" <<endl;
      exit(-1);
  }

  Site xPart(parts->lattice());
  Site xField(gij->lattice());

  typename std::list<part>::iterator it;

  double referPos[3];
  double rescalPos[3];
  double rescalPosDown[3];
  double latresolution = parts->res();

  double mass;
  double Epart;
  double gij_onPart[6];
  double cicVol;
  cicVol= latresolution*latresolution*latresolution;
  double cicVol2 = cicVol*cicVol;

  int vectorSize = gij->lattice()->vectorSize();

  LFVector<partList<part>> vecPart;
  LFVector<Real> vecgij[8][6];
  LFVector<Real> scalarCube[8]; // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110 | 111
  LFVector<Real> vectorCube[8][3];
  LFVector<Real> tensorCube[8][6];


      int sfrom;

      size_t offset_m;
      size_t offset_p;



      if(oset == NULL)
      {
          sfrom = parts->mass_type();
          offset_m = parts->mass_offset();
          offset_p = offsetof(part,vel);
      }
      else
      {
          sfrom =  flag_where;
          offset_m = oset[0];
          offset_p = oset[1];
      }

      if(sfrom == FROM_INFO)
      {
          mass = *(double*)((char*)parts->parts_info() + offset_m);
          mass /= cicVol;
          //cout << mass<<endl;
      }

      for(xPart.first(),xField.first();xPart.test();xPart.next(),xField.next())
      {
        vecPart.nocopy(parts->field().(xPart));
        // XYZ = 000 | 001 | 010 | 011 | 100 | 101 | 110tensor
        scalarCube[0].nocopy( (*E)(xField) );
        scalarCube[1].nocopy( (*E)(xField+2) );
        scalarCube[2].nocopy( (*E)(xField+1) );
        scalarCube[3].nocopy( (*E)(xField+1+2) );
        scalarCube[4].nocopy( (*E)(xField+0) );
        scalarCube[5].nocopy( (*E)(xField+0+2) );
        scalarCube[6].nocopy( (*E)(xField+0+1) );
        scalarCube[7].nocopy( (*E)(xField+0+1+2) );
        for(int i=0;i<3;i++)
        {
          vectorCube[0][i].nocopy( (*Pi)(xField,i) );
          vectorCube[1][i].nocopy( (*Pi)(xField+2,i) );
          vectorCube[2][i].nocopy( (*Pi)(xField+1,i) );
          vectorCube[3][i].nocopy( (*Pi)(xField+1+2,i) );
          vectorCube[4][i].nocopy( (*Pi)(xField+0,i) );
          vectorCube[5][i].nocopy( (*Pi)(xField+0+2,i) );
          vectorCube[6][i].nocopy( (*Pi)(xField+0+1,i) );
          vectorCube[7][i].nocopy( (*Pi)(xField+0+1+2,i) );
        }
        for(int i=0;i<6;i++)
        {
          tensorCube[0][i].nocopy( (*Sij)(xField,i) );
          tensorCube[1][i].nocopy( (*Sij)(xField+2,i) );
          tensorCube[2][i].nocopy( (*Sij)(xField+1,i) );
          tensorCube[3][i].nocopy( (*Sij)(xField+1+2,i) );
          tensorCube[4][i].nocopy( (*Sij)(xField+0,i) );
          tensorCube[5][i].nocopy( (*Sij)(xField+0+2,i) );
          tensorCube[6][i].nocopy( (*Sij)(xField+0+1,i) );
          tensorCube[7][i].nocopy( (*Sij)(xField+0+1+2,i) );
        }
        for(int i=0;i<6;i++)
        {
          vecgij[0][i].nocopy( (*gij)(xField,i) );
          vecgij[1][i].nocopy( (*gij)(xField+2,i) );
          vecgij[2][i].nocopy( (*gij)(xField+1,i) );
          vecgij[3][i].nocopy( (*gij)(xField+1+2,i) );
          vecgij[4][i].nocopy( (*gij)(xField+0,i) );
          vecgij[5][i].nocopy( (*gij)(xField+0+2,i) );
          vecgij[6][i].nocopy( (*gij)(xField+0+1,i) );
          vecgij[7][i].nocopy( (*gij)(xField+0+1+2,i) );
        }

        for(int i=0;i<3;i++)referPos[i]=xPart.coord(i)*latresolution;

        for(int v = 0;v<vectorSize;v++)
        {
          if(vecPart[v].size!=0)
          {

            for (vecPart[v].parts.begin(); it != vecPart[v].parts.end(); ++it)
            {
              for(int i =0;i<3;i++)
              {
                  rescalPos[i]=(*it).pos[i]-referPos[i];
                  rescalPosDown[i]=latresolution -rescalPos[i];
              }

              if(sfrom==FROM_PART)
              {
                  mass = *(double*)((char*)&(*it)+offset_m);
                  //mass /=cicVol;
              }

              //building gij_onPart;
              for(int i=0;i<6;i++)
              {
                gij_onPart[i] =  rescalPosDown[0]*rescalPosDown[1]*rescalPosDown[2] * vectgij[0][i][v] / cicVol;
                gij_onPart[i] += rescalPosDown[0]*rescalPosDown[1]*rescalPos[2] * vectgij[1][i][v] / cicVol;
                gij_onPart[i] += rescalPosDown[0]*rescalPos[1]*rescalPosDown[2] * vectgij[2][i][v] / cicVol;
                gij_onPart[i] += rescalPosDown[0]*rescalPos[1]*rescalPos[2] * vectgij[3][i][v] / cicVol;
                gij_onPart[i] += rescalPos[0]*rescalPosDown[1]*rescalPosDown[2] * vectgij[4][i][v] / cicVol;
                gij_onPart[i] += rescalPos[0]*rescalPosDown[1]*rescalPos[2] * vectgij[5][i][v] / cicVol;
                gij_onPart[i] += rescalPos[0]*rescalPos[1]*rescalPosDown[2] * vectgij[6][i][v] / cicVol;
                gij_onPart[i] += rescalPos[0]*rescalPos[1]*rescalPos[2] * vectgij[7][i][v] / cicVol;
              }

              Epart = mass*mass;
              for(int i=0;i<3;i++)
              {
                for(int j=0;j<3;j++)
                {
                  Epart += gij_onPart[getCompSym(i,j)] * (*it).vel[i] * (*it).vel[j] ;
                }
              }
              Epart = sqrt(Epart);

              //project E
              //000
              scalarCube[0] += rescalPosDown[0]*rescalPosDown[1]*rescalPosDown[2] * Epart / cicVol;
              //001
              scalarCube[1] += rescalPosDown[0]*rescalPosDown[1]*rescalPos[2] * Epart / cicVol;
              //010
              scalarCube[2] += rescalPosDown[0]*rescalPos[1]*rescalPosDown[2] * Epart / cicVol;
              //011
              scalarCube[3] += rescalPosDown[0]*rescalPos[1]*rescalPos[2] * Epart / cicVol;
              //100
              scalarCube[4] += rescalPos[0]*rescalPosDown[1]*rescalPosDown[2] * Epart / cicVol;
              //101
              scalarCube[5] += rescalPos[0]*rescalPosDown[1]*rescalPos[2] * Epart / cicVol;
              //110
              scalarCube[6] += rescalPos[0]*rescalPos[1]*rescalPosDown[2] * Epart / cicVol;
              //111
              scalarCube[7] += rescalPos[0]*rescalPos[1]*rescalPos[2] * Epart / cicVol;

              //project Pi
              for(int i=0;i<3;i++)
              {
                //000
                vectorCube[0][i] += rescalPosDown[0]*rescalPosDown[1]*rescalPosDown[2] * (*it).vel[i] / cicVol;
                //001
                vectorCube[1][i] += rescalPosDown[0]*rescalPosDown[1]*rescalPos[2] * (*it).vel[i] / cicVol;
                //010
                vectorCube[2][i] += rescalPosDown[0]*rescalPos[1]*rescalPosDown[2] * (*it).vel[i] / cicVol;
                //011
                vectorCube[3][i] += rescalPosDown[0]*rescalPos[1]*rescalPos[2] * (*it).vel[i] / cicVol;
                //100
                vectorCube[4][i] += rescalPos[0]*rescalPosDown[1]*rescalPosDown[2] * (*it).vel[i] / cicVol;
                //101
                vectorCube[5][i] += rescalPos[0]*rescalPosDown[1]*rescalPos[2] * (*it).vel[i] / cicVol;
                //110
                vectorCube[6][i] += rescalPos[0]*rescalPos[1]*rescalPosDown[2] * (*it).vel[i] / cicVol;
                //111
                vectorCube[7][i] += rescalPos[0]*rescalPos[1]*rescalPos[2] * (*it).vel[i] / cicVol;
              }


              //project Sij
              for(int i=0;i<3;i++)
              {
                for(int j=i; j<3;j++)
                {
                  //000
                  tensorCube[0][getCompSym(i,j)] += rescalPosDown[0]*rescalPosDown[1]*rescalPosDown[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //001
                  tensorCube[1][getCompSym(i,j)] += rescalPosDown[0]*rescalPosDown[1]*rescalPos[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //010
                  tensorCube[2][getCompSym(i,j)] += rescalPosDown[0]*rescalPos[1]*rescalPosDown[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //011
                  tensorCube[3][getCompSym(i,j)] += rescalPosDown[0]*rescalPos[1]*rescalPos[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //100
                  tensorCube[4][getCompSym(i,j)] += rescalPos[0]*rescalPosDown[1]*rescalPosDown[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //101
                  tensorCube[5][getCompSym(i,j)] += rescalPos[0]*rescalPosDown[1]*rescalPos[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //110
                  tensorCube[6][getCompSym(i,j)] += rescalPos[0]*rescalPos[1]*rescalPosDown[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                  //111
                  tensorCube[7][getCompSym(i,j)] += rescalPos[0]*rescalPos[1]*rescalPos[2] * (*it).vel[i]* (*it).vel[j] / cicVol / Epart;
                }
              }




            }
          }

          referPos[0]+=latresolution;
        }
      }
}



void DDM_CIC_project(Field<Real> * E,
                     Field<Real> * Pi,
                     Field<Real> * Sij)
{
  if(E->lattice().halo() == 0 || Pi->lattice().halo() == 0 || Sij->lattice().halo() == 0)
  {
      cout<< "LATfield2::DDM_CIC_project: the field has to have at least a halo of 1" <<endl;
      cout<< "LATfield2::DDM_CIC_project: aborting" <<endl;
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
  int jmax;
  int kmax;
  int cbuffer;


    iref = sizeLocalGross[0] - halo;

    for(int k=distHaloOne;k<sizeLocalOne[2]+distHaloOne;k++)
    {
        for(int j=distHaloOne;j<sizeLocalOne[1]+distHaloOne;j++)
        {
            (*E).value(setIndex(sizeLocalGross,halo,j,k)) += (*E).value(setIndex(sizeLocalGross,iref,j,k));
            for(int c=0;c<3;c++)(*Pi).value(setIndex(sizeLocalGross,halo,j,k),c) += (*Pi).value(setIndex(sizeLocalGross,iref,j,k),c);
            for(int c=0;c<6;c++)(*Sij).value(setIndex(sizeLocalGross,halo,j,k),c) += (*Sij).value(setIndex(sizeLocalGross,iref,j,k),c);
        }
    }

    bufferSizeY =  (long)(sizeLocalOne[2]-1) * (long)sizeLocal[0] * 10;
    bufferSizeZ = (long)sizeLocal[0] * (long)sizeLocal[1] * 10;
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
    kmax = sizeLocalOne[2]-1
    iref = sizeLocalGross[1]- halo;

    for(int k=0;k<kmax;k++)
    {
        for(int i=0;i<imax;i++)
        {
          bufferSend[i+k*imax]=(*E).value(setIndex(sizeLocalGross,i+halo,iref,k+halo));
        }
    }
    for(c = 0;c<3;c++)
    {
      cbuffer = c+1;
      for(int k=0;k<kmax;k++)
      {
          for(int i=0;i<imax;i++)
          {
            bufferSend[i+ imax*(k + kmax * cbuffer) ]=(*Pi).value(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
          }
      }
    }
    for(c = 0;c<6;c++)
    {
      cbuffer = c+4;
      for(int k=0;k<kmax;k++)
      {
          for(int i=0;i<imax;i++)
          {
            bufferSend[i+ imax*(k + kmax * cbuffer) ]=(*Sij).value(setIndex(sizeLocalGross,i+halo,iref,k+halo),c);
          }
      }
    }
    parallel.sendUp_dim1(bufferSend,bufferRec,bufferSizeY);
    //unpack data
    for(int k=0;k<kmax;k++)
    {
        for(int i=0;i<imax;i++)(*E).value(setIndex(sizeLocalGross,i+halo,halo,k+halo))+=bufferRec[i+k*imax];

    }
    for(c = 0;c<3;c++)
    {
      cbuffer = c+1;
      for(int k=0;k<kmax;k++)
      {
          for(int i=0;i<imax;i++)
            (*Pi).value(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[i+ imax*(k + kmax * cbuffer)];
      }
    }
    for(c = 0;c<6;c++)
    {
      cbuffer = c+4;
      for(int k=0;k<kmax;k++)
      {
          for(int i=0;i<imax;i++)
            (*Sij).value(setIndex(sizeLocalGross,i+halo,halo,k+halo),c)+=bufferRec[i+ imax*(k + kmax * cbuffer)];
      }
    }



    iref=sizeLocalGross[2]-halo;
    jmax = sizeLocalOne[1]-2
    for(int j=0;j<jmax;j++)
    {
        for(int i=0;i<imax;i++)
          bufferSend[i+j*imax]=(*E).value(setIndex(sizeLocalGross,i+halo,j+halo,iref));
    }
    for(c = 0;c<3;c++)
    {
      cbuffer = c+1;
      for(int j=0;j<jmax;j++)
      {
          for(int i=0;i<imax;i++)
            bufferSend[i+imax*(j+jmax*cbuffer)]=(*Pi).value(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
      }
    }
    for(c = 0;c<6;c++)
    {
      cbuffer = c+4;
      for(int j=0;j<jmax;j++)
      {
          for(int i=0;i<imax;i++)
            bufferSend[i+imax*(j+jmax*cbuffer)]=(*Sij).value(setIndex(sizeLocalGross,i+halo,j+halo,iref),c);
      }
    }



    parallel.sendUp_dim0(bufferSend,bufferRec,bufferSizeZ);


    //unpack data

    for(int j=0;j<(sizeLocalOne[1]-2);j++)
    {
        for(int i=0;i<imax;i++)(*E).value(setIndex(sizeLocalGross,i+halo,j+halo,halo))+=bufferRec[i+j*imax];
    }
    for(c = 0;c<3;c++)
    {
      cbuffer = c+1;
      for(int j=0;j<(sizeLocalOne[1]-2);j++)
      {
          for(int i=0;i<imax;i++)(*Pi).value(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[i+imax*(j+jmax*cbuffer)];
      }
    }
    for(c = 0;c<6;c++)
    {
      cbuffer = c+4;
      for(int j=0;j<(sizeLocalOne[1]-2);j++)
      {
          for(int i=0;i<imax;i++)(*Sij).value(setIndex(sizeLocalGross,i+halo,j+halo,halo),c)+=bufferRec[i+imax*(j+jmax*cbuffer)];
      }

    }


    free(bufferRec);
    free(bufferSend);
}


#endif
