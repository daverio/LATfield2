#ifndef LATFIELD2_vector_HPP
#define LATFIELD2_vector_HPP

#include <stdio.h>
//#include "vectorclass.h"


template<class T>
class LFvector
{
public:
  LFvector();
  LFvector(const LFvector& other);
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

};


template<class T>
LFvector<T>::~LFvector()
{
  if(allocated_)
  {
    free(data_);
    allocated_=false;
  }
}

template<class T>
LFvector<T>::LFvector()
{
  size_=0;
  data_=NULL;
  allocated_=false;
}

template<class T>
LFvector<T>::LFvector(const LFvector& other)
{
  size_ = other.size_;
  data_= (T*)malloc(size_*sizeof(T));
  allocated_=true;
  for(int i = 0;i<size_;i++) *(data_+i) = other.data_[i];
}

template<class T>
LFvector<T>::LFvector(Lattice & lat, T * source)
{
  size_=lat.vectorSize();
  if(source!=NULL)
  {
    data_=source;
    allocated_=false;
  }
  else
  {
    data_= (T*)malloc(size_*sizeof(T));
    allocated_=true;
  }
}

template<class T>
LFvector<T>::LFvector(Lattice & lat)
{
  this->size_ = lat.vectorSize();
}

template<class T>
void LFvector<T>::initialize(Lattice & lat)
{
  size_=lat.vectorSize();
  allocated_=false;
}

template<class T>
void LFvector<T>::initialize(Lattice & lat, T * source)
{
  size_=lat.vectorSize();
  if(source!=NULL)
  {
    data_=source;
    allocated_=false;
  }
  else
  {
    data_= (T*)malloc(size_*sizeof(T));
    allocated_=true;
  }
}

template<class T>
LFvector<T>::LFvector(int size)
{
  this->size_ = size;
}

template<class T>
void LFvector<T>::initialize(int size)
{
  size_=size;
  allocated_=false;
}

template<class T>
LFvector<T>::LFvector(int size, T * source)
{
  size_=size;
  if(source!=NULL)
  {
    data_=source;
    allocated_=false;
  }
  else
  {
    //cout<<" vector allocation "<<endl;
    data_= (T*)malloc(size_*sizeof(T));
    allocated_=true;
  }
}
template<class T>
void LFvector<T>::initialize(int size, T * source)
{
  size_=size;
  if(source!=NULL)
  {
    data_=source;
    allocated_=false;
  }
  else
  {
    //cout<<" vector allocation "<<endl;
    data_= (T*)malloc(size_*sizeof(T));
    allocated_=true;
  }
}




template<class T>
void LFvector<T>::nocopy(const LFvector& other)
{
  this->size_ = other->size_;
  this->data_ = other->data_;
}

template<class T>
void LFvector<T>::setData(T * source)
{
  this->data_ = source;
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



LFvector<double> vpow(const LFvector<double> v, double n)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = pow(v.data_[i],n);
  }
  return result;
}

LFvector<double> vsqrt(const LFvector<double> v)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = sqrt(v.data_[i]);
  }
  return result;
}


LFvector<double> vexp(const LFvector<double> v)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = exp(v.data_[i]);
  }
  return result;
}

LFvector<double> vabs(const LFvector<double> v)
{
  LFvector<double> result(v.size_,NULL);
  for(int i = 0;i<v.size_;i++)
  {
    result.data_[i] = fabs(v.data_[i]);
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

#endif
