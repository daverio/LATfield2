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
      void initialize_Field(Field<FieldType> * FieldBase, MultiField<FieldType> *& field);


      /*
      restrict level to level + 1;
      */
      template<class FieldType>
      void restrict(MultiField<FieldType> *& field, int level, int method=1){this->restrict(field,field,level,method);}

      template<class FieldType>
      void restrict(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level, int method = 1);

      template<class FieldType>
      void restrict3d_spl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);

      template<class FieldType>
      void restrict3d_dpl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);

      template<class FieldType>
      void restrict3d_spl_fw(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);

      template<class FieldType>
      void restrict3d_dpl_fw(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);

      template<class FieldType>
      void prolong(MultiField<FieldType> *& field, int level){this->prolong(field,field,level);}

      template<class FieldType>
      void prolong(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);

      template<class FieldType>
      void prolong3d_spl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);

      template<class FieldType>
      void prolong3d_dpl(MultiField<FieldType> *& src, MultiField<FieldType> *& dst, int level);


      int nl(){return nl_;}
      int npl(){return npl_;}
      int player(int level){return pLayer_[level];}
      bool isPartLayer(int i){return parallel.layer(pLayer_[i]).isPartLayer();}
      MultiLAT lattice(int i){return lattice_[i];}
private:

  MultiLAT * lattice_;

  int  npl_, nl_;
  int * pLayer_;
  bool * plr0_;
  bool * plr1_;

  int dim_;
  int ** lat_size_;
};





#endif
