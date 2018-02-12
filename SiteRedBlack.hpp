#ifndef LATFIELD2_SITEREDBLACK_HPP
#define LATFIELD2_SITEREDBLACK_HPP

/// works only if the halo is >= 1
/// no reason to do a red-black loop without halo!!!
class SiteRedBlack3d : public Site
{
public:
  SiteRedBlack3d();
  SiteRedBlack3d(Lattice& lattice);
  SiteRedBlack3d(Lattice& lattice, long index);
  void first();
  bool test();
  void next();
  void nextRed();
  void nextBlack();
  int color(){return color_;}

private:
  int color_;
};

SiteRedBlack3d::SiteRedBlack3d() {index_=0; lattice_ = NULL;}
SiteRedBlack3d::SiteRedBlack3d(Lattice& lattice) { initialize(lattice); }
SiteRedBlack3d::SiteRedBlack3d(Lattice& lattice, long index) { initialize(lattice, index); }

void SiteRedBlack3d::first()
{
  index_=lattice_->siteFirst();
  color_ = 0;
}

bool SiteRedBlack3d::test()
{
  if(color_ == 1)
  {
    return index_ <= lattice_->siteLast();
  }
  else
  {
    if(index_ >= lattice_->siteLast())
    {
      index_ = lattice_->siteFirst()+1;
      color_ = 1;
    }
    return 1;
  }
}

void SiteRedBlack3d::next()
{
  if(color_)nextBlack();
  else nextRed();
}

void SiteRedBlack3d::nextRed()
{
  index_ += 2;
  if( this->coordLocal(0) < lattice_->sizeLocal(0) && this->coordLocal(0) >= 0) { return; }
  else
  {
    if(this->coordLocal(0) != lattice_->sizeLocal(0)) index_--;
    index_ -= lattice_->sizeLocal(0);

    index_ += lattice_->jump(1);
    if( this->coordLocal(1) != lattice_->sizeLocal(1) )
    {
      if(this->coordLocal(1)%2 != 0)
      {
        if(this->coordLocal(2)%2 == 0) index_++;
      }
      else
      {
        if(this->coordLocal(2)%2 == 1) index_++;
      }
      return;
    }
    index_ -= lattice_->sizeLocal(1) * lattice_->jump(1);
    index_ += lattice_->jump(2);
    if( this->coordLocal(2) != lattice_->sizeLocal(2) )
    {
      if(this->coordLocal(1)%2 != 0)
      {
        if(this->coordLocal(2)%2 == 0) index_++;
      }
      else
      {
        if(this->coordLocal(2)%2 == 1) index_++;
      }
      return;
    }
    index_ = lattice_->siteLast() + 1;
  }

}

void SiteRedBlack3d::nextBlack()
{
  index_ += 2;
  if( this->coordLocal(0) < lattice_->sizeLocal(0) && this->coordLocal(0) >= 0) { return; }
  else
  {
    if(this->coordLocal(0) != lattice_->sizeLocal(0)) index_--;
    index_ -= lattice_->sizeLocal(0);

    index_ += lattice_->jump(1);
    if( this->coordLocal(1) != lattice_->sizeLocal(1) )
    {
      if(this->coordLocal(1)%2 != 0)
      {
        if(this->coordLocal(2)%2 == 1) index_++;
      }
      else
      {
        if(this->coordLocal(2)%2 == 0) index_++;
      }
      return;
    }
    index_ -= lattice_->sizeLocal(1) * lattice_->jump(1);
    index_ += lattice_->jump(2);
    if( this->coordLocal(2) != lattice_->sizeLocal(2) )
    {
      if(this->coordLocal(1)%2 != 0)
      {
        if(this->coordLocal(2)%2 == 1) index_++;
      }
      else
      {
        if(this->coordLocal(2)%2 == 0) index_++;
      }
      return;
    }
    index_ = lattice_->siteLast() + 1;
  }
}




#endif
