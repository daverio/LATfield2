#ifndef LATFIELD2_vector_float_HPP
#define LATFIELD2_vector_float_HPP

template<>
LFvector<float>::LFvector(const LFvector& other)
{
  size_ = other.size_;
  //posix_memalign((void**)&data_,ALIGNEMENT,size_*sizeof(float));
  data_ = (float*)malloc(size_*sizeof(float));
  allocated_=true;
  for(int i = 0;i<size_;i+=NUM_FLOATS) data_[i] = other.data_[i];
  linked_=false;
}


template<>
LFvector<float>::LFvector(Lattice & lat, float * source)
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
    //posix_memalign((void**)&data_,ALIGNEMENT,size_*sizeof(float));
    data_ = (float*)malloc(size_*sizeof(float));
    allocated_=true;
    linked_=false;
  }
}

template<>
void LFvector<float>::initialize(Lattice & lat, float * source)
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
    //posix_memalign((void**)data_,ALIGNEMENT,size_*sizeof(float));
    data_ = (float*)malloc(size_*sizeof(float));
    allocated_=true;
    linked_=false;
  }
}

template<>
LFvector<float>::LFvector(int size, float * source)
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
    //posix_memalign((void**)&(this->data_),ALIGNEMENT,size_*sizeof(float));
    data_ = (float*)malloc(size_*sizeof(float));
    this->allocated_=true;
    linked_=false;
  }
}
template<>
void LFvector<float>::initialize(int size, float * source)
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
    //posix_memalign((void**)&(this->data_),ALIGNEMENT,size_*sizeof(float));
    data_ = (float*)malloc(size_*sizeof(float));
    this->allocated_=true;
    linked_=false;
  }
}



template<>
LFvector<float>& LFvector<float>::operator=(const LFvector<float>& source)
{

  VECT_FLOAT a;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    a.load(source.data_+i);
    a.store(this->data_+i);
    //this->data_[i] = source.data_[i];
  }
  return *this;
}

template<>
LFvector<float>& LFvector<float>::operator=(LFvector<float>&& source)
{
  if(this->linked_)
  {
    VECT_FLOAT a;
    for(int i = 0;i<size_;i+=NUM_FLOATS)
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
LFvector<float>& LFvector<float>::operator=(const float& a)
{
  VECT_FLOAT av(a);

  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    av.store(this->data_+i);
  }
  return *this;
}


template<>
LFvector<float>& LFvector<float>::operator+=(const LFvector<float>& source)
{
  VECT_FLOAT r,s;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r += s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator+=(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT r;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    r.load(this->data_+i);
    r += av;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator-=(const LFvector<float>& source)
{
  VECT_FLOAT r,s;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r -= s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator-=(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT r;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    r.load(this->data_+i);
    r -= av;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator*=(const LFvector<float>& source)
{
  VECT_FLOAT r,s;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r *= s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator*=(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT r;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    r.load(this->data_+i);
    r *= av;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator/=(const LFvector<float>& source)
{
  VECT_FLOAT r,s;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    s.load(this->data_+i);
    r.load(source.data_+i);
    r /= s;
    r.store(this->data_+i);
  }
  return *this;
}
template<>
LFvector<float>& LFvector<float>::operator/=(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT r;
  for(int i = 0;i<size_;i+=NUM_FLOATS)
  {
    r.load(this->data_+i);
    r /= av;
    r.store(this->data_+i);
  }
  return *this;
}


template<>
LFvector<float> LFvector<float>::operator+(const LFvector<float>& v1)
{
  VECT_FLOAT a,b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a += b;
    a.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<float> LFvector<float>::operator+(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<result.size_;i+=NUM_FLOATS)
  {
    b.load(this->data_ + i);
    b += av;
    b.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<float> operator+( const float& a, const LFvector<float>& v1 )
{
  VECT_FLOAT av(a);
  VECT_FLOAT b;
  LFvector<float> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    b.load(v1.data_ + i);
    b += av;
    b.store(result.data_ + i);
  }
  return result;
}




template<>
LFvector<float> LFvector<float>::operator-(const LFvector<float>& v1)
{
  VECT_FLOAT a,b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a -= b;
    a.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<float> LFvector<float>::operator-(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<result.size_;i+=NUM_FLOATS)
  {
    b.load(this->data_ + i);
    b -= a;
    b.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<float> operator-( const float& a, const LFvector<float>& v1 )
{
  VECT_FLOAT av(a);
  VECT_FLOAT b,r;
  LFvector<float> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    b.load(v1.data_ + i);
    r = av - b;
    r.store(result.data_ + i);
  }
  return result;
}


template<>
LFvector<float> LFvector<float>::operator*(const LFvector<float>& v1)
{
  VECT_FLOAT a,b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a *= b;
    a.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<float> LFvector<float>::operator*(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<result.size_;i+=NUM_FLOATS)
  {
    b.load(this->data_ + i);
    b *= a;
    b.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<float> operator*( const float& a, const LFvector<float>& v1 )
{
  VECT_FLOAT av(a);
  VECT_FLOAT b;
  LFvector<float> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    b.load(v1.data_ + i);
    b *= av;
    b.store(result.data_ + i);
  }
  return result;
}


template<>
LFvector<float> LFvector<float>::operator/(const LFvector<float>& v1)
{
  VECT_FLOAT a,b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    a.load(this->data_ + i);
    b.load(v1.data_ + i);
    a /= b;
    a.store(result.data_ + i);
  }
  return result;
}
template<>
LFvector<float> LFvector<float>::operator/(const float& a)
{
  VECT_FLOAT av(a);
  VECT_FLOAT b;
  LFvector<float> result(this->size_,NULL);
  for(int i = 0;i<this->size_;i+=NUM_FLOATS)
  {
    b.load(this->data_ + i);
    b /= a;
    b.store(result.data_ + i);
  }
  return result;
}

template<>
LFvector<float> operator/( const float& a, const LFvector<float>& v1 )
{
  VECT_FLOAT av(a);
  VECT_FLOAT b,r;
  LFvector<float> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    b.load(v1.data_ + i);
    r = av/b;
    r.store(result.data_ + i);
  }
  return result;
}


template<>
float& operator+=(float& res, const LFvector<float>& v1)
{
  VECT_FLOAT b;
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    b.load(v1.data_ + i);
    res += horizontal_add(b);
  }
  return res;
}

template<>
float vsum(const LFvector<float>& v1)
{
  float result = 0.0;
  VECT_FLOAT b;
  for(int i = 0;i<v1.size_;i+=NUM_FLOATS)
  {
    b.load(v1.data_ + i);
    result += horizontal_add(b);
  }
  return result;
}

template<>
float vmax(const LFvector<float>& v1)
{
  float result = 0.0;
  for(int i = 0;i<v1.size_;i++)
  {
    if(result<v1.data_[i])result = v1.data_[i];
  }
  return result;
}



template<>
float vmin(const LFvector<float>& v1)
{
  float result = 100000000000000000.0;
  for(int i = 0;i<v1.size_;i++)
  {
    if(result>v1.data_[i])result = v1.data_[i];
  }
  return result;
}



LFvector<float> vpow(const LFvector<float> v, float n)
{
  LFvector<float> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::pow(v.data_[i],n);
  }
  return result;
}

LFvector<float> vsqrt(const LFvector<float> v)
{
  LFvector<float> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::sqrt(v.data_[i]);
  }
  return result;
}


LFvector<float> vexp(const LFvector<float> v)
{
  LFvector<float> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::exp(v.data_[i]);
  }
  return result;
}

LFvector<float> vabs(const LFvector<float> v)
{
  LFvector<float> result(v.size_,NULL);
  VECT_FLOAT b;
  for(int i = 0;i<v.size_;i+=NUM_FLOATS)
  {
    b.load(v.data_ + i);
    b = abs(b);
    b.store(result.data_ + i);
  }
  return result;
}


#endif
