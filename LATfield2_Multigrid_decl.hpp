#ifndef LATFIELD2_MULTIGRID_DECL_HPP
#define LATFIELD2_MULTIGRID_DECL_HPP



class MultiLAT : public Lattice
{
private:
    int parallel_layer_;
public:
    void initialize(int dim, const int* size, int halo, int parallel_layer);
    void initialize(int dim, const int size, int halo, int parallel_layer);

    int parallel_layer(){return parallel_layer_;}

};


class MultiGrid
{
public:
      MultiGrid();
      ~MultiGrid();

      void initialize(Lattice * lat_top, int levelNumber_max, int minGridperProc);

private:

  MultiLAT * lattice_;

  int  npl_, nl_;
  int * lLayer_;

  int dim_;
  int ** lat_size_;
};





#endif
