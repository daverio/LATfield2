#ifndef LATFIELD2_vector_double_HPP
#define LATFIELD2_vector_double_HPP

template<>
LFvector<double>::LFvector(const LFvector& other)
{
  size_ = other.size_;
  //posix_memalign((void**)&data_,ALIGNEMENT,size_*sizeof(double));
  data_ = (double*)malloc(size_*sizeof(double));
  allocated_=true;
  for(int i = 0;i<size_;i+=NUM_DOUBLES) data_[i] = other.data_[i];
  linked_=false;
}


template<>
LFvector<double>::LFvector(Lattice & lat, double * source)
{
  size_=lat.vectorSize();
  if(source!=NULL)
  {
    data_=source;
    allocated_=false;
    linked_=true;
  }
  else
  {
    //posix_memalign((void**)&data_,ALIGNEMENT,size_*sizeof(double));
    data_ = (double*)malloc(size_*sizeof(double));
    allocated_=true;
    linked_=false;
  }
}

template<>
void LFvector<double>::initialize(Lattice & lat, double * source)
{
  size_=lat.vectorSize();
  if(source!=NULL)
  {
    data_=source;
    allocated_=false;
    linked_=true;
  }
  else
  {
    //posix_memalign((void**)data_,ALIGNEMENT,size_*sizeof(double));
    data_ = (double*)malloc(size_*sizeof(double));
    allocated_=true;
    linked_=false;
  }
}

template<>
LFvector<double>::LFvector(int size, double * source)
{
  this->size_=size;
  if(source!=NULL)
  {
    this->data_=source;
    this->allocated_=false;
    linked_=true;
  }
  else
  {
    //cout<<" vector allocation "<<endl;
    //posix_memalign((void**)&(this->data_),ALIGNEMENT,size_*sizeof(double));
    data_ = (double*)malloc(size_*sizeof(double));
    this->allocated_=true;
    linked_=false;
  }
}
template<>
void LFvector<double>::initialize(int size, double * source)
{
  this->size_=size;
  if(source!=NULL)
  {
    this->data_=source;
    this->allocated_=false;
    linked_=true;
  }
  else
  {
    //cout<<" vector allocation "<<endl;
    //posix_memalign((void**)&(this->data_),ALIGNEMENT,size_*sizeof(double));
    data_ = (double*)malloc(size_*sizeof(double));
    this->allocated_=true;
    linked_=false;
  }
}



template<>
LFvector<double>& LFvector<double>::operator=(const LFvector<double>& source)
{

  VECT_DOUBLE a;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    a.load(source.data_+i);
    a.store(this->data_+i);
    //this->data_[i] = source.data_[i];
  }
  return *this;
}

template<>
LFvector<double>& LFvector<double>::operator=(LFvector<double>&& source)
{
  if(this->linked_)
  {
    VECT_DOUBLE a;
    for(int i = 0;i<size_;i+=NUM_DOUBLES)
    {
      a.load(source.data_+i);
      a.store(this->data_+i);
      //this->data_[i] = source.data_[i];
    }
  }
  else
  {
    if(allocated_)free(this->data_);
    this->size_ = source.size_;
    this->data_ = source.data_;
    this->allocated_ = source.allocated_;
    source.allocated_ = false;
  }
  return *this;
}

template<>
LFvector<double>& LFvector<double>::operator=(const double& a)
{
  VECT_DOUBLE av(a);

  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    av.store(this->data_+i);
  }
  return *this;
}


template<>
LFvector<double>& LFvector<double>::operator+=(const LFvector<double>& source)
{
  VECT_DOUBLE r,s;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r += s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator+=(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE r;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    r.load(this->data_+i);
    r += av;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator-=(const LFvector<double>& source)
{
  VECT_DOUBLE r,s;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r -= s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator-=(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE r;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    r.load(this->data_+i);
    r -= av;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator*=(const LFvector<double>& source)
{
  VECT_DOUBLE r,s;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r *= s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator*=(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE r;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    r.load(this->data_+i);
    r *= av;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator/=(const LFvector<double>& source)
{
  VECT_DOUBLE r,s;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r /= s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<double>& LFvector<double>::operator/=(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE r;
  for(int i = 0;i<size_;i+=NUM_DOUBLES)
  {
    r.load(this->data_+i);
    r /= av;
    r.store(this->data_+i);
  }
  return *this;
}


template<>
LFvector<double> LFvector<double>::operator+(const LFvector<double>& v1)
{
  VECT_DOUBLE a,b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a += b;
    a.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<double> LFvector<double>::operator+(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<result.size_;i+=NUM_DOUBLES)
  {
    b.load(this->data_ + i);
    b += av;
    b.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<double> operator+( const double& a, const LFvector<double>& v1 )
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b;
  LFvector<double> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    b.load(v1.data_ + i);
    b += av;
    b.store(result.data_ + i);
  }
  return result;
}




template<>
LFvector<double> LFvector<double>::operator-(const LFvector<double>& v1)
{
  VECT_DOUBLE a,b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a -= b;
    a.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<double> LFvector<double>::operator-(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<result.size_;i+=NUM_DOUBLES)
  {
    b.load(this->data_ + i);
    b -= a;
    b.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<double> operator-( const double& a, const LFvector<double>& v1 )
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b,r;
  LFvector<double> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    b.load(v1.data_ + i);
    r = av - b;
    r.store(result.data_ + i);
  }
  return result;
}


template<>
LFvector<double> LFvector<double>::operator*(const LFvector<double>& v1)
{
  VECT_DOUBLE a,b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a *= b;
    a.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<double> LFvector<double>::operator*(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<result.size_;i+=NUM_DOUBLES)
  {
    b.load(this->data_ + i);
    b *= a;
    b.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<double> operator*( const double& a, const LFvector<double>& v1 )
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b;
  LFvector<double> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    b.load(v1.data_ + i);
    b *= av;
    b.store(result.data_ + i);
  }
  return result;
}


template<>
LFvector<double> LFvector<double>::operator/(const LFvector<double>& v1)
{
  VECT_DOUBLE a,b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a /= b;
    a.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<double> LFvector<double>::operator/(const double& a)
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b;
  LFvector<double> result(this->size_,NULL);
  for(int i = 0;i<this->size_;i+=NUM_DOUBLES)
  {
    b.load(this->data_ + i);
    b /= a;
    b.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<double> operator/( const double& a, const LFvector<double>& v1 )
{
  VECT_DOUBLE av(a);
  VECT_DOUBLE b,r;
  LFvector<double> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    b.load(v1.data_ + i);
    r = av/b;
    r.store(result.data_ + i);
  }
  return result;
}


template<>
double& operator+=(double& res, const LFvector<double>& v1)
{
  VECT_DOUBLE b;
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    b.load(v1.data_ + i);
    res += horizontal_add(b);
  }
  return res;
}

template<>
double vsum(const LFvector<double>& v1)
{
  double result = 0.0;
  VECT_DOUBLE b;
  for(int i = 0;i<v1.size_;i+=NUM_DOUBLES)
  {
    b.load(v1.data_ + i);
    result += horizontal_add(b);
  }
  return result;
}

template<>
double vmax(const LFvector<double>& v1)
{
  double result = 0.0;
  for(int i = 0;i<v1.size_;i++)
  {
    if(result<v1.data_[i])result = v1.data_[i];
  }
  return result;
}



template<>
double vmin(const LFvector<double>& v1)
{
  double result = 100000000000000000.0;
  for(int i = 0;i<v1.size_;i++)
  {
    if(result>v1.data_[i])result = v1.data_[i];
  }
  return result;
}



LFvector<double> vpow(const LFvector<double> v, double n)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::pow(v.data_[i],n);
  }
  return result;
}

LFvector<double> vsqrt(const LFvector<double> v)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::sqrt(v.data_[i]);
  }
  return result;
}


LFvector<double> vexp(const LFvector<double> v)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::exp(v.data_[i]);
  }
  return result;
}

LFvector<double> vabs(const LFvector<double> v)
{
  LFvector<double> result(v.size_,NULL);
  VECT_DOUBLE b;
  for(int i = 0;i<v.size_;i+=NUM_DOUBLES)
  {
    b.load(v.data_ + i);
    b = abs(b);
    b.store(result.data_ + i);
  }
  return result;
}


#endif
