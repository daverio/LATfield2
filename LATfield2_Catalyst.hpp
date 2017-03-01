#ifndef CATALYST_HPP
#define CATALYST_HPP

#include <vtkCellData.h>
#include <vtkCellType.h>
#include <vtkCommunicator.h>
#include <vtkCPDataDescription.h>
#include <vtkCPInputDataDescription.h>
#include <vtkCPProcessor.h>
#include <vtkCPPythonScriptPipeline.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkImageData.h>
#include <vtkMultiBlockDataSet.h>
#include <vtkMultiPieceDataSet.h>
#include <vtkMultiProcessController.h>
#include <vtkNew.h>
#include <vtkPoints.h>
#include <vtkPointData.h>


class CataLAT
{
public:
CataLAT();
void initialize(int paraview_scripts_number, int paraview_scripts_offset,char** argv);
void finalize();
void createVTKgrid(Lattice *lat,double dx);
void addField(string arrayName, Field<double> *source);
void copyField(string arrayName, Field<double> *source);
void coProcess(double time, unsigned int timeStep, Lattice * lat);

private:

vtkCPProcessor* processor_; // static data
vtkMultiBlockDataSet* vtkGrid_;
//vtkCPDataDescription* dataDescription_;

};

CataLAT::CataLAT()
{
  processor_ = NULL;
  vtkGrid_ = NULL;
}

void CataLAT::initialize(int paraview_scripts_number, int paraview_scripts_offset,char** argv)
{


  if(processor_ == NULL)
  {
    processor_ = vtkCPProcessor::New();
    processor_->Initialize();
  }
  else
  {
    processor_->RemoveAllPipelines();
  }
  // scripts are passed in as command line arguments
  for(int i=paraview_scripts_offset;i<paraview_scripts_number+paraview_scripts_offset;i++)
  {
    vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
    pipeline->Initialize(argv[i]);
    processor_->AddPipeline(pipeline);
    pipeline->Delete();
  }

}

void CataLAT::createVTKgrid(Lattice *lat, double dx)
{

    vtkGrid_ = vtkMultiBlockDataSet::New();
    vtkNew<vtkImageData> imageData;

    imageData->SetExtent(0,lat->sizeLocal(0),
                         lat->coordSkip()[1],lat->coordSkip()[1]+lat->sizeLocal(1),
                         lat->coordSkip()[0],lat->coordSkip()[0]+lat->sizeLocal(2));
    imageData->SetSpacing(dx,dx,dx);
    imageData->SetOrigin(0,0,0);

    vtkNew<vtkMultiPieceDataSet> multiPiece;
    multiPiece->SetNumberOfPieces(parallel.size());
    multiPiece->SetPiece(parallel.rank(), imageData.GetPointer());
    vtkGrid_->SetNumberOfBlocks(1);
    vtkGrid_->SetBlock(0, multiPiece.GetPointer());

}
void CataLAT::addField(string arrayName, Field<double> *source)
{
  //method called by both server and compute;

  if(vtkGrid_==NULL)
  {
    COUT<<"CataLAT::addField : the grid must have been created before calling CataLat::addField... aborting"<<endl;
    exit(1);
  }
  else
  {

    vtkMultiPieceDataSet* multiPiece = vtkMultiPieceDataSet::SafeDownCast(vtkGrid_->GetBlock(0));
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(multiPiece->GetPiece(parallel.rank()));

    if(dataSet->GetCellData()->HasArray(arrayName.c_str())==0)
    {

      vtkNew<vtkDoubleArray> newarray;
      newarray->SetName(arrayName.c_str());
      newarray->SetNumberOfComponents(source->components());
      newarray->SetNumberOfTuples(static_cast<vtkIdType>(source->lattice().sitesLocal()));
      dataSet->GetCellData()->AddArray(newarray.GetPointer());

    }
    else
    {
      COUT<<"CataLAT::addField : trying to add a field which does already exist: "<<arrayName<<endl;
    }
  }


}
void CataLAT::copyField(string arrayName, Field<double> *source)
{

  //method called by both server and compute;
  if(vtkGrid_==NULL)
  {
    COUT<<"CataLAT::copyField : the grid must have been created before calling CataLat::copyField... aborting"<<endl;
    if(parallel.root())parallel.abortForce();
  }
  else
  {
    vtkMultiPieceDataSet* multiPiece = vtkMultiPieceDataSet::SafeDownCast(vtkGrid_->GetBlock(0));
    vtkDataSet* dataSet = vtkDataSet::SafeDownCast(multiPiece->GetPiece(parallel.rank()));

    if(dataSet->GetCellData()->HasArray(arrayName.c_str())==1)
    {

      //copying the field
      Site x(source->lattice());
      long i;
      vtkDoubleArray* outarray = vtkDoubleArray::SafeDownCast(dataSet->GetCellData()->GetArray(arrayName.c_str()));
      double values[source->components()];

      for(x.first(),i=0;x.test();x.next(),i++)
      {
        for(int c=0;c<source->components();c++)values[c]=(*source)(x,c);
        outarray->SetTupleValue(i, values);
      }

    }
    else
    {
      COUT<<"CataLAT::copyField : trying to copy a field which does not already exist"<<endl;
      COUT<<"CataLAT::copyField : creating the field: "<<arrayName<<endl;
      this->addField(arrayName,source);
      this->copyField(arrayName,source);
    }

  }

}

void CataLAT::coProcess(double time, unsigned int timeStep, Lattice *lat)
{

  //cout<<"blabla"<<endl;
  vtkNew<vtkCPDataDescription> dataDescription;
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time, timeStep);

  dataDescription->ForceOutputOn();

  if(processor_->RequestDataDescription(dataDescription.GetPointer()) != 0)
  {

      dataDescription->GetInputDescriptionByName("input")->SetGrid(vtkGrid_);
      dataDescription->GetInputDescriptionByName("input")->SetWholeExtent(0,lat->size(0),
                                                                          0,lat->size(1),
                                                                          0,lat->size(2));


      processor_->CoProcess(dataDescription.GetPointer());
  }

}


void CataLAT::finalize()
{
  if(processor_)
    {
    processor_->Delete();
    processor_ = NULL;
    }

  if(vtkGrid_)
  {
      vtkGrid_->Delete();
      vtkGrid_ = NULL;
  }



}


/*
vtkCPProcessor* Processor = NULL; // static data

void CatalystInit(int paraview_scripts_number, int paraview_scripts_offset,char** argv)
{
  if(Processor == NULL)
  {
    Processor = vtkCPProcessor::New();
    Processor->Initialize();
  }
  // scripts are passed in as command line arguments
  for(int i=paraview_scripts_offset;i<paraview_scripts_number;i++)
  {
    vtkCPPythonScriptPipeline* pipeline = vtkCPPythonScriptPipeline::New();
    pipeline->Initialize(argv[i]);
    Processor->AddPipeline(pipeline);
    pipeline->Delete();
  }
}

void CatalystFinalize()
{
  if(Processor)
    {
    Processor->Delete();
    Processor = NULL;
    }
}


void CatalystCoProcess(int timeStep, double time, Field<double> * rho)
{
  vtkCPDataDescription* dataDescription = vtkCPDataDescription::New();
  dataDescription->AddInput("input");
  dataDescription->SetTimeData(time, timeStep);
  if(Processor->RequestDataDescription(dataDescription) != 0)
  {
// Catalyst needs to output data
// Create an axis-aligned, uniform grid
    vtkUniformGrid* grid = vtkUniformGrid::New();
    grid->SetWholeExtent(0, rho->lattice().size(0),
                         0, rho->lattice().size(1),
                         0, rho->lattice().size(2));
    grid->SetExtents(0,rho->lattice().localSize(0),
                     rho->lattice().coordSkip[1],rho->lattice().coordSkip[1]+rho->lattice().localSize(1),
                     rho->lattice().coordSkip[0],rho->lattice().coordSkip[0]+rho->lattice().localSize(2));




    dataDescription->GetInputDescriptionByName("input")->SetGrid(grid);
    grid->Delete();
// Create a field associated with points
    vtkDoubleArray* array = vtkDoubleArray::New(); array->SetName("data");
    //array->SetArray(field, grid->GetNumberOfPoints(), 1); grid->GetPointData()->AddArray(array);
    //array->Delete();
    Processor->CoProcess(dataDescription);
  }
  dataDescription->Delete();
}


*/
#endif
