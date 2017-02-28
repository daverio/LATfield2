#ifndef LATFIELD2_PLANFFT_TRANSPOSE_HPP
#define LATFIELD2_PLANFFT_TRANSPOSE_HPP


//transposition function

#ifndef FFTW_CPLX
#ifdef SINGLE
#define FFTW_CPLX fftwf_complex
#else
#define FFTW_CPLX fftw_complex
#endif
#endif

/////
void transpose_0_2( FFTW_CPLX * in, FFTW_CPLX * out,int dim_i,int dim_j ,int dim_k)
{
  int i,j,k;
  for(i=0;i<dim_i;i++)
  {
    for(j=0;j<dim_j;j++)
    {
      for(k=0;k<dim_k;k++)
      {

        out[k+dim_k*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[k+dim_k*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];

      }
    }
  }

}

void transpose_0_2_last_proc( FFTW_CPLX * in, FFTW_CPLX * out,int dim_i,int dim_j ,int dim_k)
{
  int i,j,k;
  for(i=0;i<dim_i;i++)
  {
    for(j=0;j<dim_j;j++)
    {
      for(k=0;k<dim_k;k++)
      {
        out[k+(dim_k+1)*(j+i*dim_j)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[k+(dim_k+1)*(j+i*dim_j)][1]=in[i+dim_i*(j+k*dim_j)][1];
      }
    }
  }
}

void implement_local_0_last_proc( FFTW_CPLX * in, FFTW_CPLX * out,int proc_dim_i,int proc_dim_j,int proc_dim_k,int proc_size)
{
  int i_in,i_out,j,rank;
  for(i_in=0;i_in<proc_dim_i;i_in++)
  {
    for(j=0;j<proc_dim_j;j++)
    {
      for(rank=0;rank<proc_size;rank++)
      {
        i_out=i_in+rank*proc_dim_i;
        out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][0]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][0];
        out[proc_dim_k + (proc_dim_k+1)*(j+i_out*proc_dim_j)][1]=in[i_in+proc_dim_i*(j+proc_dim_j*rank)][1];
      }
    }
  }
}


void transpose_1_2(FFTW_CPLX * in , FFTW_CPLX * out ,int dim_i,int dim_j ,int dim_k )
{
  int i,j,k;
  for(i=0;i<dim_i;i++)
  {
    for(j=0;j<dim_j;j++)
    {
      for(k=0;k<dim_k;k++)
      {
        out[i+dim_i*(k+j*dim_k)][0]=in[i+dim_i*(j+k*dim_j)][0];
        out[i+dim_i*(k+j*dim_k)][1]=in[i+dim_i*(j+k*dim_j)][1];
      }
    }
  }
}

void transpose_back_0_3( FFTW_CPLX * in, FFTW_CPLX * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size,int halo,int components, int comp)
{
  int i,j,k,l, i_t, j_t, k_t;
  int r2c_halo = r2c + 2*halo;
  int local_size_k_halo = local_size_k + 2*halo;
  for (i=0;i<local_r2c;i++)
  {
    for(k=0;k<local_size_k;k++)
    {
      for(j=0;j<local_size_j;j++)
      {
        for(l=0;l<proc_size;l++)
        {
          i_t = i + l*local_r2c;
          j_t = j ;
          k_t = k ;
          out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
          out[comp+components*(i_t + r2c_halo * (k_t + local_size_k_halo * j_t))][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
        }
      }
    }
  }
}

void implement_0(FFTW_CPLX * in, FFTW_CPLX * out,int r2c_size,int local_size_j,int local_size_k, int halo,int components, int comp)
{
  int i,j,k;
  i=r2c_size-1;
  int r2c_halo = r2c_size + 2*halo;
  int local_size_k_halo = local_size_k + 2*halo;

  for(j=0;j<local_size_j;j++)
  {
    for(k=0;k<local_size_k;k++)
    {
      out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][0]=in[j + local_size_j *k][0];
      out[comp+components*(i + r2c_halo * (k + local_size_k_halo *j))][1]=in[j + local_size_j *k][1];
    }
  }

}


void b_arrange_data_0(FFTW_CPLX *in, FFTW_CPLX * out,int dim_i,int dim_j ,int dim_k, int khalo, int components, int comp)
{
  int i,j,k;
  int jump_i=(dim_i+ 2 *khalo);
  int jump_j=dim_j+ 2 *khalo;
  for(i=0;i<dim_i;i++)
  {
    for(j=0;j<dim_j;j++)
    {
      for(k=0;k<dim_k;k++)
      {
        out[j + dim_j * (k + dim_k * i)][0]=in[comp+components*(i + jump_i * (j + jump_j*k))][0];
        out[j + dim_j * (k + dim_k * i)][1]=in[comp+components*(i + jump_i * (j + jump_j*k))][1];
      }
    }
  }

}

void b_transpose_back_0_1( FFTW_CPLX * in, FFTW_CPLX * out,int r2c,int local_r2c,int local_size_j,int local_size_k,int proc_size)
{
  int i,j,k,l, i_t, j_t, k_t;

  for (i=0;i<local_r2c;i++)
  {
    for(k=0;k<local_size_k;k++)
    {
      for(j=0;j<local_size_j;j++)
      {
        for(l=0;l<proc_size;l++)
        {
          i_t = i + l*local_r2c;
          j_t = j ;
          k_t = k ;
          out[i_t + r2c * (k_t + local_size_k * j_t)][0]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][0];
          out[i_t + r2c * (k_t + local_size_k * j_t)][1]=in[i + local_r2c * (j + local_size_j * (k + local_size_k *l)) ][1];
        }
      }
    }
  }
}

void b_implement_0(FFTW_CPLX * in, FFTW_CPLX * out,int r2c_size,int local_size_j,int local_size_k)
{
  int i,j,k;
  i=r2c_size-1;


  for(j=0;j<local_size_j;j++)
  {
    for(k=0;k<local_size_k;k++)
    {
      out[i + r2c_size * (k + local_size_k *j)][0]=in[j + local_size_j *k][0];
      out[i + r2c_size * (k + local_size_k *j)][1]=in[j + local_size_j *k][1];
    }
  }

}




#endif
