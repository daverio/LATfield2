#ifndef LATFIELD2_MULTIGRID_DECL_HPP
#define LATFIELD2_MULTIGRID_DECL_HPP


class MultiGrid
{
public:
      MultiGrid();
      ~MultiGrid();

      void initialize(Lattice * lat_top, int levelNumber_max, int minGridperProc);

private:

  Lattice * lattice_;

  int  pl_number_, nl_;


  int dim_;
  int ** lat_size_;



};



#endif
