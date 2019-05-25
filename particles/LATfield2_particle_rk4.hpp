#ifndef LATFIELD2_PARTICLE_RK4_DEF_HPP
#define LATFIELD2_PARTICLE_RK4_DEF_HPP

/*! \file LATfield2_particle_simple.hpp
 \brief Description of the particle type: "part_simple"

 */


/*type table

 datatype      memType         fileType
 int            H5T_NATIVE_INT   INT_TYPE_FILE
 long           H5T_NATIVE_LONG  LONG_TYPE_FILE
 float          H5T_NATIVE_FLOAT  FLOAT_TYPE_FILE
 double         H5T_NATIVE_DOUBLE  DOUBLE_TYPE_FILE
 bool           H5T_NATIVE_HBOOL  BOOL_TYPE_FILE
 */


/**
 * \addtogroup partdesc
 * @{
 */

/**
 * \addtogroup part RK4 "part_RK4"
 * @{
 */

// part RK4 definition

/*! \struct part_rk4
 \brief individual properties of the particle type "part_rk4".

 "part_simple" contains the minimal list of individual properties is: ID, position,velocity,
 */
struct part_rk4{

  long ID;
  Real mass;
  LATfield2::Real pos[3];
  LATfield2::Real vel[3];
  LATfield2::Real pos_out[3];
  LATfield2::Real vel_out[3];
  LATfield2::Real pos_in[3];
  LATfield2::Real vel_in[3];
  // interesting one need in fact only 3 fields, the point is that here we cannot save computation by using an 4th field, pos_dot and vel_dot.
};
/*!
 \brief overloading of the << operator for individual property strucutre.
 \return ostream containing the ID, position and velocity of the particle.
 */
ostream& operator<<(ostream& os, const part_rk4& p)
{
    os << "ID: "<<p.ID<<", mass: "<< p.mass<<" , Pos: ("<< p.pos[0]<<","<< p.pos[1]<<","<< p.pos[2]<<") , Vel: (" << p.vel[0]<<","<< p.vel[1]<<","<< p.vel[2]<<")";
    return os;
}

/*! \struct part_simple_info
 \brief individual properties of the particle type "part_simple".

 "part_simple" contains the minimal list of global properties is: type_name and size of type_name. But also contain two additionnal global properties: "mass" and "relativistic".
 */
struct part_rk4_info{
    int type_name_size;
    char  type_name[64];

};


#ifdef HDF5
/*! \struct part_simple_dataType
 \brief properties datatype structure.

 Structure which contain the HDF5 datatype of all properties.
 */
struct part_rk4_dataType{
  hid_t part_memType;
  hid_t part_fileType;
  hid_t part_info_memType;
  hid_t part_info_fileType;

  part_rk4_dataType(){

    hid_t strtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (strtype, 64);

    part_memType = H5Tcreate(H5T_COMPOUND, sizeof (part_rk4));
    H5Tinsert(part_memType, "ID", HOFFSET (part_rk4, ID), H5T_NATIVE_LONG);
    H5Tinsert(part_memType, "mass", HOFFSET (part_rk4, mass), REAL_TYPE);
    H5Tinsert(part_memType, "positionX", HOFFSET (part_rk4, pos[0]), REAL_TYPE);
    H5Tinsert(part_memType, "positionY", HOFFSET (part_rk4, pos[1]), REAL_TYPE);
    H5Tinsert(part_memType, "positionZ", HOFFSET (part_rk4, pos[2]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityX", HOFFSET (part_rk4, vel[0]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityY", HOFFSET (part_rk4, vel[1]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityZ", HOFFSET (part_rk4, vel[2]), REAL_TYPE);

    part_info_memType = H5Tcreate(H5T_COMPOUND, sizeof (part_rk4_info));
    H5Tinsert(part_info_memType, "type_name_size", HOFFSET (part_rk4_info, type_name_size),INT_TYPE_FILE );
    H5Tinsert(part_info_memType, "type_name", HOFFSET (part_rk4_info, type_name), strtype);


    part_fileType = H5Tcreate (H5T_COMPOUND, sizeof(long) + 7 * sizeof(Real) );
    H5Tinsert(part_fileType, "ID"       ,0  ,LONG_TYPE_FILE);
    H5Tinsert(part_fileType, "mass"     ,sizeof(long)                    ,LONG_TYPE_FILE);
    H5Tinsert(part_fileType, "positionX",sizeof(long) + 1 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "positionY",sizeof(long) + 2 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "positionZ",sizeof(long) + 3 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "velocityX",sizeof(long) + 4 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "velocityY",sizeof(long) + 5 * sizeof(Real) ,REAL_TYPE_FILE);
    H5Tinsert(part_fileType, "velocityZ",sizeof(long) + 6 * sizeof(Real) ,REAL_TYPE_FILE);

    part_info_fileType = H5Tcreate(H5T_COMPOUND, sizeof(int)+ 64);
    H5Tinsert(part_info_fileType, "type_name_size", 0 ,INT_TYPE_FILE );
    H5Tinsert(part_info_fileType, "type_name", sizeof(int), strtype);

    H5Tclose (strtype);
  }

};
#endif




/**@}*/

/**@}*/


#endif
