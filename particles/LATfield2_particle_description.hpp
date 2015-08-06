#ifndef LATFIELD2_PARTICLE_DEF_HPP
#define LATFIELD2_PARTICLE_DEF_HPP

/*type table
 
 datatype      memType         fileType
 
 
 
 
 */



// part_simple definition
struct part_simple{
  
  long ID;
  LATfield2::Real pos[3];
  LATfield2::Real vel[3];
};
ostream& operator<<(ostream& os, const part_simple& p)
{
    os << "ID: "<<p.ID<<" , Pos: ("<< p.pos[0]<<","<< p.pos[1]<<","<< p.pos[2]<<") , Vel: (" << p.vel[0]<<","<< p.vel[1]<<","<< p.vel[2]<<")"; 
    return os;
}

struct part_simple_info{
    double  mass;
    bool relativistic;
    char * type_name;
};

#ifdef HDF5
struct part_simple_dataType{
  hid_t part_memType;
  hid_t part_fileType;
  hid_t part_info_memType;
  hid_t part_info_fileType;

  part_simple_dataType(){

    hid_t strtype = H5Tcopy (H5T_C_S1);
    H5Tset_size (strtype, H5T_VARIABLE);
    
    
    part_memType = H5Tcreate(H5T_COMPOUND, sizeof (part_simple)); 
    H5Tinsert(part_memType, "ID", HOFFSET (part_simple, ID), H5T_NATIVE_LONG);
    H5Tinsert(part_memType, "positionX", HOFFSET (part_simple, pos[0]), REAL_TYPE);
    H5Tinsert(part_memType, "positionY", HOFFSET (part_simple, pos[1]), REAL_TYPE);
    H5Tinsert(part_memType, "positionZ", HOFFSET (part_simple, pos[2]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityX", HOFFSET (part_simple, vel[0]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityY", HOFFSET (part_simple, vel[1]), REAL_TYPE);
    H5Tinsert(part_memType, "velocityZ", HOFFSET (part_simple, vel[2]), REAL_TYPE);
    


    part_info_memType = H5Tcreate(H5T_COMPOUND, sizeof (part_simple_info));
    H5Tinsert(part_info_memType, "mass", HOFFSET (part_simple_info, mass), H5T_NATIVE_DOUBLE);
    H5Tinsert(part_info_memType, "relativistic", HOFFSET (part_simple_info, relativistic), H5T_NATIVE_HBOOL);
    H5Tinsert(part_info_memType, "type_name", HOFFSET (part_simple_info, type_name), strtype);

    
#ifdef SINGLE
    part_fileType = H5Tcreate (H5T_COMPOUND, 8 + (3*4) + (3*4) );
    H5Tinsert(part_fileType, "ID"       ,0          ,H5T_STD_I64BE);
    H5Tinsert(part_fileType, "positionX",8          ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "positionY",8+4        ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "positionZ",8+4+4      ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "velocityX",8+4+4+4    ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "velocityY",8+4+4+4+4  ,H5T_IEEE_F32BE);
    H5Tinsert(part_fileType, "velocityZ",8+4+4+4+4+4,H5T_IEEE_F32BE);
      
    part_info_fileType = H5Tcreate(H5T_COMPOUND, 8 + sizeof(hvl_t)+1);
    H5Tinsert(part_info_fileType, "mass", 0 ,H5T_IEEE_F64BE );
    H5Tinsert(part_info_fileType, "relativistic",8, H5T_STD_I8BE);
    H5Tinsert(part_info_fileType, "type_name",8+1, strtype);
#else
    part_fileType = H5Tcreate (H5T_COMPOUND, 8 + (3*8) + (3*8) );
    H5Tinsert(part_fileType, "ID"       ,0          ,H5T_STD_I64BE);
    H5Tinsert(part_fileType, "positionX",8          ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "positionY",8+8        ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "positionZ",8+8+8      ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "velocityX",8+8+8+8    ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "velocityY",8+8+8+8+8  ,H5T_IEEE_F64BE);
    H5Tinsert(part_fileType, "velocityZ",8+8+8+8+8+8,H5T_IEEE_F64BE);
    
    part_info_fileType = H5Tcreate(H5T_COMPOUND, 8 + sizeof(hvl_t)+1);
    H5Tinsert(part_info_fileType, "mass", 0 ,H5T_IEEE_F64BE );
    H5Tinsert(part_info_fileType, "relativistic",8, H5T_STD_I8BE);
    H5Tinsert(part_info_fileType, "type_name",8+1, strtype);
#endif

    H5Tclose (strtype);
  }
    
};
#endif









#endif
