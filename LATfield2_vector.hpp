#ifndef LATFIELD2_vector_HPP
#define LATFIELD2_vector_HPP




#include <stdio.h>
#include "vectorclass.h"


template<class T>
class LFvector
{
public:
  LFvector();
  LFvector(const LFvector& other);
  LFvector(LFvector&& other);
  LFvector(Lattice & lat, T * source);
  LFvector(int size, T * source);
  LFvector(Lattice & lat);
  LFvector(int size);

  void initialize(Lattice & lat);
  void initialize(Lattice & lat, T * source);
  void initialize(int size);
  void initialize(int size, T * source);
  void nocopy(const LFvector& other);
  void setData(T * source);

  ~LFvector();

  T& operator[](int i);
  LFvector<T>& operator=(const LFvector<T>& source);
  LFvector<T>& operator=(LFvector<T>&& source);
  LFvector<T>& operator=(const T& a);
  LFvector<T>& operator+=(const LFvector<T>& source);
  LFvector<T>& operator+=(const T& a);
  LFvector<T>& operator-=(const LFvector<T>& source);
  LFvector<T>& operator-=(const T& a);
  LFvector<T>& operator*=(const LFvector<T>& source);
  LFvector<T>& operator*=(const T& a);
  LFvector<T>& operator/=(const LFvector<T>& source);
  LFvector<T>& operator/=(const T& a);

  LFvector<T> operator+(const LFvector<T>& v1);
  LFvector<T> operator+(const T& a);
  LFvector<T> operator-(const LFvector<T>& v1);
  LFvector<T> operator-(const T& a);
  LFvector<T> operator*(const LFvector<T>& v1);
  LFvector<T> operator*(const T& a);
  LFvector<T> operator/(const LFvector<T>& v1);
  LFvector<T> operator/(const T& a);

  T * data_;
  int size_;
  bool allocated_;
  bool linked_;
};


template<class T>
LFvector<T>::~LFvector()
{
  if(allocated_)free(data_);
}

template<class T>
LFvector<T>::LFvector()
{
  size_=0;
  data_=NULL;
  allocated_=false;
  linked_=false;
}

template<class T>
LFvector<T>::LFvector(const LFvector& other)
{
  size_ = other.size_;
  data_= (T*)malloc(size_*sizeof(T));
  allocated_=true;
  for(int i = 0;i<size_;i++) *(data_+i) = other.data_[i];
  linked_=false;
}

template<class T>
LFvector<T>::LFvector(LFvector&& other)
{
  if(allocated_)free(data_);

  size_ = other.size_;
  data_ = other.data_;
  linked_=other.linked_;
  allocated_ = other.allocated_;
  other.allocated_ = false;
  other.data_ = nullptr;

}

template<class T>
LFvector<T>::LFvector(Lattice & lat, T * source)
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
    data_= (T*)malloc(size_*sizeof(T));
    allocated_=true;
    linked_=false;
  }
}

template<class T>
LFvector<T>::LFvector(Lattice & lat)
{
  this->size_ = lat.vectorSize();
  allocated_=false;
  linked_=false;
}

template<class T>
void LFvector<T>::initialize(Lattice & lat)
{
  size_=lat.vectorSize();
  allocated_=false;
  linked_=false;
}

template<class T>
void LFvector<T>::initialize(Lattice & lat, T * source)
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
    data_= (T*)malloc(size_*sizeof(T));
    allocated_=true;
    linked_=false;
  }
}

template<class T>
LFvector<T>::LFvector(int size)
{
  this->size_ = size;
  allocated_=false;
  linked_=false;
}

template<class T>
void LFvector<T>::initialize(int size)
{
  this->size_=size;
  allocated_=false;
  linked_=false;
}

template<class T>
LFvector<T>::LFvector(int size, T * source)
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
    this->data_= (T*)malloc(size_*sizeof(T));
    this->allocated_=true;
    linked_=false;
  }
}
template<class T>
void LFvector<T>::initialize(int size, T * source)
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
    this->data_= (T*)malloc(size_*sizeof(T));
    this->allocated_=true;
    linked_=false;
  }
}




template<class T>
void LFvector<T>::nocopy(const LFvector& other)
{
  this->size_ = other->size_;
  this->data_ = other->data_;
  this->allocated_ = false;
  linked_=true;
}

template<class T>
void LFvector<T>::setData(T * source)
{
  if(allocated_)free(this->data_);
  this->data_ = source;
  this->allocated_ = false;
  linked_=true;
}

template<class T>
inline T& LFvector<T>::operator[](int i)
{
  return data_[i];
}


template<class T>
LFvector<T>& LFvector<T>::operator=(const LFvector<T>& source)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] = source.data_[i];
  }
  return *this;
}

template<class T>
LFvector<T>& LFvector<T>::operator=(LFvector<T>&& source)
{
  if(this->linked_)
  {
    for(int i = 0;i<size_;i++)
    {
      this->data_[i] = source.data_[i];
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

template<class T>
LFvector<T>& LFvector<T>::operator=(const T& a)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] = a;
  }
  return *this;
}


template<class T>
LFvector<T>& LFvector<T>::operator+=(const LFvector<T>& source)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] += source.data_[i];
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator+=(const T& a)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] += a;
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator-=(const LFvector<T>& source)
{

  for(int i = 0;i<size_;i++)
  {
    this->data_[i] -= source.data_[i];
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator-=(const T& a)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] -= a;
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator*=(const LFvector<T>& source)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] *= source.data_[i];
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator*=(const T& a)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] *= a;
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator/=(const LFvector<T>& source)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] /= source.data_[i];
  }
  return *this;
}
template<class T>
LFvector<T>& LFvector<T>::operator/=(const T& a)
{
  for(int i = 0;i<size_;i++)
  {
    this->data_[i] /= a;
  }
  return *this;
}


template<class T>
LFvector<T> LFvector<T>::operator+(const LFvector<T>& v1)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<v1.size_;i++)
  {
    result.data_[i] = this->data_[i] + v1.data_[i];
  }
  return result;
}

template<class T>
LFvector<T> LFvector<T>::operator+(const T& a)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i] + a;
  }
  return result;
}

template<class T>
LFvector<T> operator+( const T& a, const LFvector<T>& v1 ) {
  LFvector<T> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i++)
  {
    result.data_[i] = v1.data_[i] + a;
  }
  return result;
}




template<class T>
LFvector<T> LFvector<T>::operator-(const LFvector<T>& v1)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i]-v1.data_[i];
  }
  return result;
}
template<class T>
LFvector<T> LFvector<T>::operator-(const T& a)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i] - a;
  }
  return result;
}
template<class T>
LFvector<T> operator-( const T& a, const LFvector<T>& v1 ) {
  LFvector<T> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i++)
  {
    result.data_[i] = a - v1.data_[i];
  }
  return result;
}


template<class T>
LFvector<T> LFvector<T>::operator*(const LFvector<T>& v1)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i]*v1.data_[i];
  }
  return result;
}
template<class T>
LFvector<T> LFvector<T>::operator*(const T& a)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i] * a;
  }
  return result;
}

template<class T>
LFvector<T> operator*( const T& a, const LFvector<T>& v1 ) {
  LFvector<T> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i++)
  {
    result.data_[i] = v1.data_[i] * a;
  }
  return result;
}


template<class T>
LFvector<T> LFvector<T>::operator/(const LFvector<T>& v1)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i]/v1.data_[i];
  }
  return result;
}
template<class T>
LFvector<T> LFvector<T>::operator/(const T& a)
{
  LFvector<T> result(this->size_,NULL);
  for(int i = 0;i<size_;i++)
  {
    result.data_[i] = this->data_[i] / a;
  }
  return result;
}

template<class T>
LFvector<T> operator/( const T& a, const LFvector<T>& v1 ) {
  LFvector<T> result(v1.size_,NULL);
  for(int i = 0;i<v1.size_;i++)
  {
    result.data_[i] =  a / v1.data_[i];
  }
  return result;
}


template<class T>
T& operator+=(T& res, const LFvector<T>& v1){
      for(int i = 0;i<v1.size_;i++)
      {
        res += v1.data_[i];
      }
      return res;
}

template<class T>
T vsum(const LFvector<T>& v1)
{
  T result = 0.0;
  for(int i = 0;i<v1.size_;i++)
  {
    result += v1.data_[i];
  }
  return result;
}

template<class T>
T vmax(const LFvector<T>& v1)
{
  T result = 0.0;
  for(int i = 0;i<v1.size_;i++)
  {
    if(result<v1.data_[i])result = v1.data_[i];
  }
  return result;
}



template<class T>
T vmin(const LFvector<T>& v1)
{
  T result = 1000000000000000.0;
  for(int i = 0;i<v1.size_;i++)
  {
    if(result>v1.data_[i])result = v1.data_[i];
  }
  return result;
}

ostream& operator<<(ostream& os, const LFvector<double> &v)
{
  for(int i = 0;i<v.size_-1;i++)
  {
    os <<  v.data_[i] << " , ";
  }
  os << v.data_[v.size_-1];
  return os;
}

#ifdef NO_VECTORIZATION

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
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = std::fabs(v.data_[i]);
  }
  return result;
}

#else
  #include "LATfield2_vector_double.hpp"
  #include "LATfield2_vector_float.hpp"
#endif
#endif
