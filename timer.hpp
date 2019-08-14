#ifndef TIMER_HPP
#define TIMER_HPP


string second2time(double in)
{
  int days;
  stringstream ss;
  double num;
  days = floor(in / 86400);
  if(days != 0) ss<<days<<"-";
  num = fmod(in,86400);
  ss << floor(num/3600)<<":";
  ss << floor(fmod(num,3600)/60)<<":"<<fmod(num,60);
  return ss.str();
}

double time2second(string time_str)
{
  size_t pos;
  int days=0;
  int hms[3]={0,0,0};
  int ptr;
  size_t size = time_str.length();
  string snum;

  pos = time_str.find('-');
  if(pos!=string::npos)
  {
      snum = "";
      snum.insert(0,&time_str[0],pos);

      stringstream ss(snum);
      ss>>days;
      pos++;
  }
  else
  {
    pos = 0;
  }

  ptr = 0;
  for(int i=pos;i<size;i++)
  {
    if(time_str[i]==':')
    {
      snum = "";
      snum.insert(0,&time_str[pos],i-pos);
      pos = i+1;
      stringstream ss(snum);
      ss>>hms[ptr++];
    }
  }

  snum = "";
  snum.insert(0,&time_str[pos],size-pos);
  stringstream ss(snum);
  ss>>hms[ptr];

  return ((days*24 + hms[0])*60 + hms[1])*60 + hms[2];
}


class MPI_timer{

public:
  MPI_timer();
  MPI_timer(int n);

  ~MPI_timer();

  void start(int i);
  void stop(int i);
  double getTime(int i);

  void initialize(int n);

  double timer(int i){return timers_[i];}
  double aveTimer(int i){if(parallel.rank()==0)return timers_[i]/counts_[i];else return 0;}
  double count(int i){if(parallel.rank()==0) return counts_[i];else return 0;}

private:
  int n_;
  double * starts_;
  double * timers_;
  int * counts_;
  bool * runs_;
};
MPI_timer::MPI_timer()
{
  n_ = 0;
}
MPI_timer::~MPI_timer()
{
  if(n_>0 && parallel.rank()==0)
  {
    delete[] starts_;
    delete[] timers_;
    delete[] counts_;
  }
}
MPI_timer::MPI_timer(int n)
{
  this->initialize(n);
}
void MPI_timer::initialize(int n)
{
  if(n>0 && parallel.rank()==0)
  {
    n_ = n;
    starts_ = new double[n_];
    timers_ = new double[n_];
    counts_ = new int[n_];
    runs_   = new bool[n_];
    for(int i = 0;i<n_;i++)
    {
      timers_[i] = 0;
      counts_[i] = 0;
      runs_[i] = false;
    }
  }
  else
  {
    n_ = 0;
  }
}

void MPI_timer::start(int i)
{
  if(parallel.rank()==0)
  {
    if(runs_[i])
    {
      COUT<<"trying to starts_ the timer "<< i <<" which already runs_!"<<endl;
    }
    else
    {
      starts_[i] = MPI_Wtime();
      runs_[i] = true;
    }
  }
}
void MPI_timer::stop(int i)
{
  if(parallel.rank()==0)
  {
    if(!runs_[i])
    {
      COUT<<"trying to stop the timer "<< i <<" which does not run!"<<endl;
    }
    else
    {
      timers_[i] += MPI_Wtime() - starts_[i];
      counts_[i]++;
      runs_[i] = false;
    }
  }
}
double MPI_timer::getTime(int i)
{
  if(parallel.rank()==0)
  {
    if(!runs_[i])
    {
      COUT<<"trying to get time from timer "<< i <<" which does not run!"<<endl;
      return -1;
    }
    else
    {
      return MPI_Wtime() - starts_[i];
    }
  }
  else
  {
    return -1;
  }
}

#endif
