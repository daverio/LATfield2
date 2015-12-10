#ifndef LATFIELD2_PARTICLE_DEF_HPP
#define LATFIELD2_PARTICLE_DEF_HPP


/*! \file LATfield2_particle_description_template.hpp
 \brief Template to create a particle description
 
 */

/*type table
 
 datatype      memType     fileType     size
 int            H5T_NATIVE_INT   INT_TYPE_FILE
 long           H5T_NATIVE_LONG  INT_TYPE_FILE
 float          H5T_NATIVE_FLOAT  FLOAT_TYPE_FILE
 double         H5T_NATIVE_DOUBLE  DOUBLE_TYPE_FILE
 bool           H5T_NATIVE_HBOOL  BOOL_TYPE_FILE
 
 
 
 */


/**
 * \addtogroup partdesc
 * @{
 */


/**
 * \addtogroup partTemplate Particle description templates
 * @{
 */



/*! \struct particles_name
 \brief individual properties structure.
 
 Structure which contain all individual properties. The minimal list of individual properties is: ID, position,velocity,
 */
struct particles_name{
  
  long ID;
  LATfield2::Real pos[3];
  LATfield2::Real vel[3];
};

/*!
    \brief overloading of the << operator for individual property strucutre.
    \return ostream containing the ID, position and velocity of the particle.
 */
ostream& operator<<(ostream& os, const particles_name& p)
{
    os << "ID: "<<p.ID<<" , Pos: ("<< p.pos[0]<<","<< p.pos[1]<<","<< p.pos[2]<<") , Vel: (" << p.vel[0]<<","<< p.vel[1]<<","<< p.vel[2]<<")"; 
    return os;
}

/*! \struct particles_name_info
 \brief global properties structure.
 
 Structure which contain all global properties. The minimal list of global properties is: type_name and size of type_name.
 */
struct particles_name_info{
    
    char * type_name;
    int type_name_size;
};

#ifdef HDF5
/*! \struct particles_name_dataType
 \brief properties datatype structure.
 
 Structure which contain the HDF5 datatype of all properties.
 */
struct particles_name_dataType{
  hid_t part_memType;
  hid_t part_fileType;
  hid_t part_info_memType;
  hid_t part_info_fileType;

  particles_name_dataType(){

    hid_t strtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (strtype, H5T_VARIABLE);
    
//memory datatype
    part_memType = H5Tcreate(H5T_COMPOUND, sizeof (particles_name)); 
    H5Tinsert(part_memType, "ID", HOFFSET (particles_name, ID), H5T_NATIVE_LONG);
    H5Tinsert(part_memType, "positionX", HOFFSET (particles_name, pos[0]), REAL_TYPE);
    H5Tinsert(part_memType, "positionY", HOFFSET (particles_name, pos[1]), REAL_TYPE);
    H5Tinsert(part_memType, "positionZ", HOFFSET (particles_name, pos[2]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityX", HOFFSET (particles_name, vel[0]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityY", HOFFSET (particles_name, vel[1]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityZ", HOFFSET (particles_name, vel[2]), REAL_TYPE);
    

    part_info_memType = H5Tcreate(H5T_COMPOUND, sizeof (particles_name_info));
    H5Tinsert(part_info_memType, "type_name", HOFFSET (particles_name_info, type_name), strtype);


//file datatype
#ifdef SINGLE
    part_fileType = H5Tcreate (H5T_COMPOUND, 8 + (3*4) + (3*4) );
    H5Tinsert(part_fileType, "ID"       ,0          ,H5T_STD_I64BE);
    H5Tinsert(part_fileType, "positionX",8          ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "positionY",8+4        ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "positionZ",8+4+4      ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "velocityX",8+4+4+4    ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "velocityY",8+4+4+4+4  ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "velocityZ",8+4+4+4+4+4,H5T_IEEE_F32BE);
      
    part_info_fileType = H5Tcreate(H5T_COMPOUND, sizeof(hvl_t));
    H5Tinsert(part_info_fileType, "type_name",0, strtype);
#else
    part_fileType = H5Tcreate (H5T_COMPOUND, 8 + (3*8) + (3*8) );
    H5Tinsert(part_fileType, "ID"       ,0          ,H5T_STD_I64BE);
    H5Tinsert(part_fileType, "positionX",8          ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "positionY",8+8        ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "positionZ",8+8+8      ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "velocityX",8+8+8+8    ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "velocityY",8+8+8+8+8  ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "velocityZ",8+8+8+8+8+8,H5T_IEEE_F64BE);
    
    part_info_fileType = H5Tcreate(H5T_COMPOUND, sizeof(hvl_t));
    H5Tinsert(part_info_fileType, "type_name",0, strtype);
#endif

    H5Tclose (strtype);
  }
    
};
#endif








/**@}*/


#endif
