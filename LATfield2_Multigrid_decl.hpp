#ifndef LATFIELD2_MULTIGRID_DECL_HPP
#define LATFIELD2_MULTIGRID_DECL_HPP



#include "LATfield2_MultiLAT.hpp"
#include "LATfield2_MultiField.hpp"



class MultiGrid
{
public:
      MultiGrid();
      ~MultiGrid();

      void initialize(Lattice * lat_top, int levelNumber_max, int minGridperProc);

      template<class FieldType>
      void intitialize_Field(Field<FieldType> * FieldBase, MultiField<FieldType> *& field);


      /*
      restrict level to level + 1;
      */
      template<class FieldType>
      void restrict(MultiField<FieldType> *& field, int level);

      template<class FieldType>
      void restrict3d_spl(MultiField<FieldType> *& field, int level);

      template<class FieldType>
      void restrict3d_dpl(MultiField<FieldType> *& field, int level);


private:

  MultiLAT * lattice_;

  int  npl_, nl_;
  int * lLayer_;

  int dim_;
  int ** lat_size_;
};





#endif
