#ifndef LATFIELD2_SETTINGSFILE_HPP
#define LATFIELD2_SETTINGSFILE_HPP


/*! \file LATfield2_SettingsFile.hpp
 \brief LATfield2_SettingsFile.hpp contain the class SettingsFile definition.
 \author N. Bevis
 */ 

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

//CLASS PROTOTYPE======================


/*! \class SettingsFile  
 \brief A utility class designed to make reading in runtime parameter values easier.
 

 If the command-line arguments are input via optional inputs on
 either the constructor or open member function, then these
 take preceident: they are effectively last in the file. 
 
 Note, when used with std::string objects, only one word
 is allowed per setting, ie. spaces are not allowed. This
 is because of the way that the >> operator works for
 this class. This fits nicely with the command-line override, however.
 
 Note that the string specified followed by = is searched
 for in the file and then the input read. If one setting
 name is also the end of another that precedes it in the file
 then the wrong one will be read.
 
 Only the primary MPI process is able to create or add to the setting
 file. Further processes will be sent the file contents via
 MPI. To use this class in serial code the preprocessor defintion SERIAL must be set. This flag have not been remove to allow users to use it outside LATfield2. LATfield2 have no serial version, therefor setting prepocessor flag -DSERIAL should never be used with LATfield2.
 
 */
class SettingsFile
{
private:
  std::string filename_;
  std::fstream file_;
  std::stringstream stream_;
  int mode_;
  bool isRoot_;    //is process the root one (false if non-parallel)                     
  //Search=================================
  bool search(const std::string search_string);
	
public:
  static int noCreate;
  static int autoCreate;


	//Constructors=======================
	//! Constructor
    SettingsFile();
    
    /*!
     Constructor + open a file
     \param filename : path to the file.
     \param mode     : noCreate (the read method will exit if the parameter does not exist) or autoCreate (read will add the missing parameter).
     \param argc     : additionnal argument number.
     \param argv     : pointer to the additionnal arguments.
     */
	SettingsFile(const std::string filename, const int mode, const int argc = 0, char** argv = NULL);
	
	//Destructor=========================
	//! desctructor
    ~SettingsFile();
	
	//File open / close / create ========
	/*!
     Open an existinge settings file
     \param filename : path to the file
     \param mode     : noCreate (the read method will exit if the parameter does not exist) or autoCreate (read will add the missing parameter).
     \param argc     : additional argument number.
     \param argv     : pointer to the additional arguments.
     */
    void open(const std::string filename, const int mode, const int argc = 0, char** argv = NULL);
    /*!
     Close the current settings file
     */
	void close();
    
    /*!
     Create a new settings file and open it.
     \param filename: path to the file.
     */
	void create(const std::string filename);
	
	//Settings read / write==================
    
    /*!
     Method to read a parameter.
     
     \param parameter_name : string containing the name of the parameter. If the parameter does not exite and the mode autocreate is set, this method will add the parameter to the settings file with the current value of "parameter". In the case the mode is set to nocreate, then read will exit for security, for this reason in it is always advise to set the mode to nocreate for production runs.
     \param parameter     : pointer to the variable where the parameter will be assigned. 
     */
	template<class TemplateClass>
	void read(const std::string parameter_name, TemplateClass& parameter);
    
    /*!
     Method to add a parameter to the settings file. The new parameter will be just added to the end of the file, even if it already exists.
     
     \param parameter_name: string containing the name of the parameter.
     \param parameter: pointer to the value of the parameter.
     */
	template<class TemplateClass>
	void add(const std::string parameter_name, const TemplateClass& parameter);
	
    /*!
     Method to write a parameter in the settings file. If the parameter_name exist, it will overwrite the parameter. And if it does not exist in the file, it will be added at the end of the file.
     
     \param parameter_name: string containing the name of the parameter
     \param parameter: pointer to the value of the parameter.
     */
    template<class TemplateClass>
	void write(const std::string parameter_name, const TemplateClass& parameter);


 
};


//CONSTANTS===========================
int SettingsFile::noCreate = 1;
int SettingsFile::autoCreate = 0;


//CONSTRUCTORS========================
SettingsFile::SettingsFile() 
{

#ifndef SERIAL
  isRoot_=parallel.isRoot();
#else
  isRoot_=true;
#endif

}

SettingsFile::SettingsFile(const std::string filename, const int mode, const int argc, char** argv) 
{

#ifndef SERIAL
  isRoot_=parallel.isRoot();
#else
  isRoot_=true;
#endif

  this->open(filename, mode, argc, argv);
}

//DESTRUCTOR==========================
SettingsFile::~SettingsFile() {this->close();}

//OPEN================================
void SettingsFile::open(const std::string filename, const int mode, const int argc, char** argv)
{
  char c;

  filename_=filename;
  mode_=mode;

  if(isRoot_)
    {      
      //Open file  
      file_.open(filename_.c_str(), std::fstream::in);
      if(!file_.is_open())
	  {
		  if((mode_ & SettingsFile::noCreate) == 0)
		  {
			  std::cout<<"SettingsFile: "<<filename_<<" not found."<<std::endl;
			  std::cout<<"SettingsFile: Creating..."<<std::endl;
			  this->create(filename_);
			  std::cout<<"creating ok"<<std::endl;
		  }
		  else
		  {
			  std::cout<<"SettingsFile: "<<filename_<<" not found and auto create off."<<std::endl;
			  std::cout<<"SettingsFile: Exiting..." << std::endl;
#ifndef SERIAL
			  parallel.abortRequest();
#else
			  exit(555);
#endif
		  }
	  }	
  
      //Read command line into stringstream
      for(int i=0; i<argc; i++)
		{
			for(int j=0; argv[i][j]!='\0'; j++)
			{
				stream_<<argv[i][j];
			}
			stream_<<endl;    
		} 
  
      //Read file into stringstream
      while(!file_.eof())
	  {
		  c=file_.get();
		  if(c=='#')
		  {
			  while(!file_.eof() && c!='\n') { c=file_.get(); }    
			  if(file_.eof()) { break; }
		  }
		  stream_.put(c);       
	  }
	  file_.close();  
    } 


#ifndef SERIAL
  //Broadcast results to all processes 
	
	
	parallel.barrier();
	
  if(parallel.size()>1)
    {
      if(parallel.isRoot())   
	  {
		  int len = stream_.str().length();
		  char* streamString = new char[len+1];
		  for(int i=0;i<=len;i++) { streamString[i]=stream_.str()[i]; }
		  parallel.broadcast(len, parallel.root());
		  parallel.broadcast(streamString, len+1, parallel.root());
	  }    
      else
	  {
		  int len;
		  char* streamString;
		  parallel.broadcast(len, parallel.root());
		  streamString=new char[len+1];
		  parallel.broadcast(streamString, len+1, parallel.root());
		  stream_<<streamString;
	  }    
	}
    
#endif
}

//FILE CLOSE============================
void SettingsFile::close()
{
  if(isRoot_) { filename_="."; }
}

//FILE CREATE===========================
void SettingsFile::create(const std::string filename)
{
    
  if(isRoot_)
    {
      filename_=filename;
      mode_=autoCreate;
      
      file_.open(filename_.c_str(), std::fstream::out);
      if(!file_.is_open())
	  {
		  std::cout<<"SettingsFile: Cannot create: "<<filename<<std::endl;
		  std::cout<<"SettingsFile: Exiting..."<<std::endl;
#ifndef SERIAL
		  parallel.abortRequest();
#else
		  exit(555);
#endif	
	  }
      else
	  {
		  file_.close();
		  file_.clear();
		  file_.open(filename.c_str(), std::fstream::in);
	  }
    }

  //parallel.barrier();
}

//PARAMETER READ===========================
template<class TemplateClass>
void SettingsFile::read(const std::string parameterName, TemplateClass &parameter)
{
  if(this->search(parameterName+'='))
    {
      stream_>>parameter;
    }
  else
    { 
      if(isRoot_)
	  {
		  //verifiy that the parameter name is no autocreate
		  if(parameterName!="autocreate")
		  {
			  std::cout<<"SettingsFile: "<<filename_<<" has no parameter: "<<parameterName<<std::endl;
			  std::cout<<"SettingsFile: No command-line override given either"<<std::endl;  
			  
			  if((mode_ & SettingsFile::noCreate) == 0 )
			  {
				  std::cout << "SettingsFile: Adding with current value: " << parameter << std::endl; 
				  this->add(parameterName, parameter);
			  }
			  else
			  {
				  std::cout << "SettingsFile: Auto create off. Exiting..." << std::endl;
#ifndef SERIAL
				  parallel.abortRequest();
#else
				  exit(555);
#endif
			  }
		  }		  
		  
	  }
#ifndef SERIAL
		parallel.barrier();
#endif
    }
}

//PARAMETER WRITE===========================
template<class TemplateClass>
void SettingsFile::add(const std::string parameter_name, const TemplateClass &parameter)
{
  if(isRoot_)
    {
      file_.clear();
      file_.open(filename_.c_str(), std::ios::out | std::ios::app);
      file_ << parameter_name << '=' << parameter << std::endl;
      if(!file_.good())
	  {
	  std::cout << "SettingsFile: Could not write to file: " << filename_ << std::endl;
	  std::cout << "SettingsFile: Exiting... " << std::endl;
#ifndef SERIAL
		parallel.abortRequest();
#else
		exit(555);
#endif
	  }
      file_.close();
    }

}
template<class TemplateClass>
void SettingsFile::write(const std::string parameter_name, const TemplateClass &parameter)
{
	if(isRoot_)
    {
		unsigned int i=0;
		unsigned int l=0;
		int line_number=0;
		
		char c;
		string line_temp;
		
        if(file_.good())file_.close();
		file_.clear();
		file_.open(filename_.c_str(), std::ios::in);
		
		
		//get number of non void line
		file_.seekg(0);
		while (!file_.eof())
		{
			getline(file_,line_temp);
			line_number++;
		}
		
		//creat line array and this_param array
        string * lines;
        lines = new string[line_number];
        
		bool this_param[line_number];
		file_.seekg(0);
		file_.close();
		file_.open(filename_.c_str(), std::ios::in);
		
		//implement non void lines
		l=0;
		while (!file_.eof())
		{
			getline(file_,line_temp);
			
            lines[l]=line_temp;
            l++;
			
		}
		for(l=0;l<line_number;l++)
		{
			for(i=0;i<parameter_name.length();i++)
			{
				c = lines[l][i];
				if(c==parameter_name[i])this_param[l]=true;
				else 
				{
					this_param[l]=false;
					i=parameter_name.length();
				}
			}
		}
		file_.close();
		file_.open(filename_.c_str(), std::ios::out | std::ios::trunc);
		
		
		bool writen;
		for(l=0;l<line_number;l++)
		{
			if(!this_param[l])file_<<lines[l]<<endl;
			else 
			{
				file_<< parameter_name <<'='<<parameter<<endl;
				writen = true;
			}
		}
		
		if(!writen)file_<< parameter_name <<'='<<parameter<<endl;
			
		
		file_.close();
		file_.clear();

    }

#ifndef SERIAL
	parallel.barrier();
#endif
	
}


//SEARCH=====================================
bool SettingsFile::search(const std::string searchString)
{
  unsigned int i=0;
  char c;
  //Set to beginning of file
  stream_.clear(); //clear any errors from having failed a previous search
  stream_.seekg(0);
 
  //Search
  while(stream_.good() && i<searchString.length())
    {
      c=stream_.get();
      if(c==searchString[i]){i++;}
      else{i=0;}
    }
  return stream_.good();
}

#endif
