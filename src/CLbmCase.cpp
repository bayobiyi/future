#include <iostream>
#include <string.h>
#include "CLbmCase.h"
#include "CLbmSolve.h"
#include "CLbmSolve3D.h"


CLbmCase::CLbmCase(void)
{
	m_pDomain= NULL;
}
//tecplot utility functions
void 
CLbmCase::Lbm_write_tecplot(char* tecplotFileName)
{
	//write case information to file attached to m_exportfs file stream

	//char mybuf[10];
	//m_exportfs.rdbuf()->pubsetbuf(mybuf,10);
	m_exportfs.open(tecplotFileName, std::ios_base::out);
	if(!m_exportfs)
	{
		std::cerr<<"File could not be opened\n";
		std::exit(1);
	}

	if(m_model.m_phaseModel == CLbmModel::MULTIPHASE)
		setFileHeaderMulti(m_exportfs, tecplotFileName);
	else
		setFileHeader(m_exportfs, tecplotFileName);

	setZoneHeader(m_exportfs);
	setCoordData(m_exportfs);

	if(m_model.m_phaseModel == CLbmModel::MULTIPHASE)
		setVarDataMulti(m_exportfs);
	else
		setVarData(m_exportfs);

	setConnData(m_exportfs);
	std::cout<<"Finished exporting tecplot file\n";
	m_exportfs.close();
}

//tecplot utility functions
void 
CLbmCase::setFileHeader(std::fstream& outputFile, char* tecplotFileName)
{
	outputFile<<"#HEADER FILE\n"
						<<"TITLE= "<<"\""<<tecplotFileName<<"\"\n"
						<<"VARIABLES="<<"\"X\""<<" "<<"\"Y\""<<" "<<"\"Z\""<<" "
						<<"\"X Velocity\""<<" "<<"\"Y Velocity\""<<" "<<"\"Z Velocity\""<<" "
						<<"\"Velocity Magnitude\""<<" "<<"\"Density\""<<" "<<"\"Uoverall-X\""<<" "<<"\"Uoverall-Y\""<<" "<<"\"Uoverall-Z\""<<" "<<"\"Pressure\"\n";
	outputFile.flush();
	std::cout<<"Finished setting zone header file\n";
}

void 
CLbmCase::setFileHeaderMulti(std::fstream& outputFile, char* tecplotFileName)
{
	outputFile<<"#HEADER FILE\n"
						<<"TITLE= "<<"\""<<tecplotFileName<<"\"\n"
						<<"VARIABLES="<<"\"X\""<<" "<<"\"Y\""<<" "<<"\"Z\""<<" "
						<<"\"X Velocity-0\""<<" "<<"\"X Velocity-1\""<<" "<<"\"Y Velocity-0\""<<" "<<"\"Y Velocity-1\""<<" "<<"\"Z Velocity-0\""<<" "<<"\"Z Velocity-1\""<<" "
						<<"\"Velocity Magnitude-0\""<<" "<<"\"Velocity magnitude-1\""<<" " <<"\"Density-0\""<<" "<<"\"Density-1\""<<" "
						<<"\"Uoverall-X\""<<" "<<"\"Uoverall-Y\""<<" "<<"\"Uoverall-Z\""<<" "<<"\"Pressure\"\n";
	outputFile.flush();
	std::cout<<"Finished setting zone header file\n";
}

void 
CLbmCase::setZoneHeader(std::fstream& outputFile)
{
	cgsize_t twoD(2);
	if(m_base.cell_dim == twoD)
	{
		/*Zone Record*/
		outputFile<<"#Zone record\n"
						<<"ZONE T= "<<"\""<<m_zone.m_pzoneName<<"\""<<","<<" "<<"N="<<m_zone.m_size[0]<<","<<" "
						<<"E="<<m_zone.m_size[1]<<","<<" "<<"ZONETYPE="<<"FEQUADRILATERAL"<<","<<" "<<"DATAPACKING="<<"BLOCK"<<","
						<<" "<<"AUXDATA"<<" "<< m_zone.m_pzoneName <<"= "<<"\""<<m_zone.m_pzoneName<<"_variables\"\n";
	}

	else
	{
		/*Zone Record*/
		outputFile<<"#Zone record\n"
							<<"ZONE T= "<<"\""<<m_zone.m_pzoneName<<"\""<<","<<" "<<"N="<<m_zone.m_size[0]<<","<<" "
							<<"E="<<m_zone.m_size[1]<<","<<" "<<"ZONETYPE="<<"FEBRICK"<<","<<" "<<"DATAPACKING="<<"BLOCK"<<","
							<<" "<<"AUXDATA"<<" "<< m_zone.m_pzoneName <<"= "<<"\""<<m_zone.m_pzoneName<<"_variables\"\n";
	}
	outputFile.flush();
}
void 
CLbmCase::setCoordData(std::fstream& outputFile)
{
	//write coordinates
	for(cgsize_t j=1; j <= m_grid.m_nCoords; j++)
	{
		for(cgsize_t i=1; i<= m_grid.m_NVertex; i++)
		{
			outputFile<<*(*(m_grid.m_coord + j-1) + i-1)<<" ";
			outputFile.flush();
			if(!(i % 2000))
				outputFile<<"\n";	
		}
		outputFile<<"\n";
	}
}
void 
CLbmCase::setVarData(std::fstream& outputFile)
{
	//write Variables
	cgsize_t node;
	for(Domain::VariableType_t::size_type i=0; i < m_pDomain->GetDomainVariables()->size(); i++)
	{
		for(CVariable::MaterialVarType_t::size_type j=0; j < m_pDomain->GetDomainVariables()->at(i)->GetVariable()->size(); j++)
		{
			node = 1;
			for(CVariable::NodeType_t::size_type k=0; k < m_pDomain->GetDomainVariables()->at(i)->GetVariable()->at(j)->size(); k++)
			{
				outputFile<< m_pDomain->GetDomainVariables()->at(i)->GetVariable()->at(j)->at(k)->m_nodeVal<<" ";
				outputFile.flush();

				if(!(node % 2000))
					outputFile<<"\n";
				node++;
			}
			outputFile<<"\n";
		}
		outputFile<<"\n";
	}
	outputFile<<"\n";
}

void 
CLbmCase::setVarDataMulti(std::fstream& outputFile)
{
	//write Variables
	cgsize_t node(0), secondary(0);
	for(Domain::VariableType_t::size_type i=0; i < m_pDomain->GetDomainVariables()->size()-6; i++)
	{
		for(CVariable::MaterialVarType_t::size_type j=0; j < m_pDomain->GetDomainVariables()->at(i)->GetVariable()->size(); j++)
		{
			node = 1;
			for(CVariable::NodeType_t::size_type k=0; k < m_pDomain->GetDomainVariables()->at(i)->GetVariable()->at(j)->size(); k++)
			{
				outputFile<< m_pDomain->GetDomainVariables()->at(i)->GetVariable()->at(j)->at(k)->m_nodeVal<<" ";
				outputFile.flush();

				if(!(node % 2000))
					outputFile<<"\n";
				node++;
			}
			outputFile<<"\n";
		}
		outputFile<<"\n";
	}
	outputFile<<"\n";

	for(Domain::VariableType_t::size_type i=m_pDomain->GetDomainVariables()->size()-6; i < m_pDomain->GetDomainVariables()->size()-3; i++)
	{
		node = 1;
		for(CVariable::NodeType_t::size_type k=0; k < m_pDomain->GetDomainVariables()->at(i)->GetVariable()->at(0)->size(); k++)
			{
				outputFile<< m_pDomain->GetDomainVariables()->at(i)->GetVariable()->at(secondary)->at(k)->m_nodeVal<<" ";
				outputFile.flush();

				if(!(node % 2000))
					outputFile<<"\n";
				node++;
			}
			outputFile<<"\n";
	}
	outputFile<<"\n";
}

void 
CLbmCase::setConnData(std::fstream& outputFile)
{
	cgsize_t twoD(2);
	if(m_base.cell_dim == twoD)
	{
		cgsize_t quad(4);
		cgsize_t rows(0); 
		for(std::map  <cgsize_t,CLbmElement::MaterialRecord_t>::const_iterator iter = this->m_element.m_sectionsMaterial.begin(); 
			iter != this->m_element.m_sectionsMaterial.end(); iter++)
		{
			rows  = (*(this->m_element.m_elementDataSize + iter->first)) / quad;
			for(cgsize_t j = 1; j <= rows; j++)
			{
				for (cgsize_t n=1; n <= quad; n++)
				{
					outputFile<< *(*(this->m_element.m_conn + iter->first) + (j-1)*quad + (n-1))<<" ";
					outputFile.flush();	
				}
				outputFile<<"\n";
			}
		}
	}

	else
	{
		cgsize_t hexa(8);
		cgsize_t rows(0); 
		for(std::map  <cgsize_t,CLbmElement::MaterialRecord_t>::const_iterator iter = this->m_element.m_sectionsMaterial.begin(); 
			iter != this->m_element.m_sectionsMaterial.end(); iter++)
		{
			rows  = (*(this->m_element.m_elementDataSize + iter->first)) / hexa;
			for(cgsize_t j = 1; j <= rows; j++)
			{
				for (cgsize_t n=1; n <= hexa; n++)
				{
					outputFile<< *(*(this->m_element.m_conn + iter->first) + (j-1)*hexa + (n-1))<<" ";
					outputFile.flush();	
				}
				outputFile<<"\n";
			}
		}
	}
	//////////////////////////////////////////////////////////////////////////
	outputFile<<
     "DATASETAUXDATA  EXPERIMENTDATE= "<<"\"April 25, 2014, 1:16am\"\n"<< 
     "VARAUXDATA 3"<<" "<<m_zone.m_pzoneName<<"= "<<"\"x-velocity\"\n"<<
     "VARAUXDATA 4"<<" "<<m_zone.m_pzoneName<<"= "<<"\"y-velocity\""<<"\n";
	std::cout<<"Finished setting connectivity file\n";
}

void
CLbmCase::ExportCaseInfo(char* transferName)
{
	//write case information to file attached to m_exportfs file stream
	//strcat(transferName,".trf");
	m_exportfs.open(transferName, std::ios::out | std::ios::trunc | std::ios::binary);
	if(!m_exportfs)
	{
		std::cerr<<"File could not be opened for export\n";
		std::exit(1);
	}
	m_base.Lbm_write(m_exportfs);
	m_model.Lbm_write(m_exportfs);
	m_zone.Lbm_write(m_exportfs);
	m_grid.Lbm_write(m_exportfs);
	m_element.Lbm_write(m_exportfs);
	m_boco.Lbm_write(m_exportfs);
	m_patchedNodes.Lbm_write(m_exportfs);
	m_sol.Lbm_write(m_exportfs);
	/*if (m_model.m_initStatus == 'i')
	{
		m_pDomain->Lbm_write(m_exportfs);
	}
	else;*/
	std::cout<<"Finished exporting transfer file\n"; 
	m_exportfs.close();
}

void
CLbmCase::ImportCaseInfo(char* transferName)
{
	m_importfs.open(transferName,std::ios::in | std::ios::binary);
	if(!m_importfs)
	{
		std::cerr<<"File could not be opened for import\n";
		std::exit(1);
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_base \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_base.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_model \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_model.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_zone \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_zone.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_grid \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_grid.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_element \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_element.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_boco \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_boco.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_patchedNodes \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_patchedNodes.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to read m_sol \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_sol.Lbm_read(m_importfs);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" now initializing block \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

CLbmCase::~CLbmCase()
{
	std::cerr<<"CLbmCase destructor called\n";
}

///////////////////////////////////////////////////////////
CLbmBase::CLbmBase()
{
	index_file       = 0;
	index_base       = 0;
  	nBases           = 0; 
	cell_dim         = 0;
	phys_dim         = 0;
	m_pbaseName      = new  char [5000];
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmBase \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void  CLbmBase::Lbm_write (std::fstream& outputFile)
{
	outputFile << index_file <<' '<< index_base <<' '<< nBases <<' '<< cell_dim <<' '<< phys_dim <<' '<< m_pbaseName <<' ';
}

void  CLbmBase::Lbm_read (std::fstream& inputFile)
{
	inputFile >> index_file >> index_base >>  nBases >> cell_dim >> phys_dim >> m_pbaseName;
}

CLbmBase::~CLbmBase()
{
	std::cerr<<"CLbmBase destructor called\n";
}


CLbmModel::CLbmModel()
{
	m_solver     = D2Q7;
	m_phaseModel = OFF;
	m_relaxTime  = SRT;
	m_energyEqn  = false;
	m_gravity    = false;
	m_initStatus = 'n';
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmBModel \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void  CLbmModel::Lbm_write (std::fstream& outputFile)
{
	outputFile <<static_cast<int>(m_solver) <<' '<< static_cast<int>(m_phaseModel) <<' '<< static_cast<int>(m_relaxTime) 
		<<' '<< m_energyEqn <<' '<< m_gravity <<' '<< m_initStatus <<' ';
}

void  CLbmModel::Lbm_read (std::fstream& inputFile)
{
	int solver     = 0;
	int phaseModel = 0;
	int relaxTime  = 0;
	inputFile >> solver >> phaseModel >> relaxTime >> m_energyEqn >> m_gravity >> m_initStatus;	
	m_solver     = static_cast<SOLVER>     (solver);
	m_phaseModel = static_cast<PHASEMODEL> (phaseModel);
	m_relaxTime  = static_cast<RELAXTIME>  (relaxTime);
}

CLbmModel::~CLbmModel()
{
	std::cerr<<"CLbmModel destructor called\n";
}
///////////////////////////////////////////////////////////////////////////////////

CLbmZone::CLbmZone()
{
	int size = 3;
	m_size = new cgsize_t[size];
	m_index_dim  = 0;
	m_index_zone = 0;
	m_pzoneName   = new char[5000];
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmZone \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void  CLbmZone::Lbm_write (std::fstream& outputFile)
{
	outputFile <<m_pzoneName<<' ';
	outputFile << m_index_dim <<' '<<m_nZones<<' '<< m_index_zone <<' '<< static_cast<int> ( m_zoneType) <<' ';
	int size = 3;
	for(cgsize_t i=1; i <= size; i++)
	{
		outputFile << m_size[i-1] <<' ' ;
	}
}

void  CLbmZone::Lbm_read (std::fstream& inputFile)
{
	int zoneType = 0;
	inputFile >> m_pzoneName;
	inputFile >> m_index_dim >> m_nZones >> m_index_zone >> zoneType;
	m_zoneType = static_cast<ZoneType_t> (zoneType);
	int size = 3;
	for(cgsize_t i=1; i <= size; i++)
	{
		inputFile >> m_size[i-1];
	}
}

CLbmZone::~CLbmZone()
{
	SAFE_DELETE_ARRAY(m_size);
	SAFE_DELETE_ARRAY(m_pzoneName);
	std::cerr<<"CLbmZone destructor called\n";
}

////////////////////////////////////////////////////////////////////////////
CLbmGrid::CLbmGrid()
{
	m_nGrids    = 0;
	m_NVertex   = 0;
	m_gridIndex  = 0;
  	m_pgridName  = new char [5000];
	m_nCoords   = 0;
	m_dataType  = NULL;
	m_coord     = NULL;
	m_coordIndex = NULL;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmGrid \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

void CLbmGrid::CreateGridMemory(int size)
{
	m_coord      = new double*    [size];
	m_pcoordName = new char*      [size];
	m_dataType   = new DataType_t [size];
	m_coordIndex = new cgsize_t   [size];
}
void  CLbmGrid::Lbm_write (std::fstream& outputFile)
{
	outputFile << m_nGrids <<' ' << m_nCoords <<' '<<m_NVertex <<' ';
	outputFile << m_pgridName <<' ';
	outputFile << m_gridIndex <<' ';
	for(cgsize_t i=1; i <= m_nCoords; i++)
	{
		outputFile << m_pcoordName[i-1] <<' ';
		outputFile << m_coordIndex[i-1] <<' ' << static_cast<int>( m_dataType[i-1]) <<' ';

		for(cgsize_t j=1; j <= m_NVertex; j++)
		{
				outputFile << *(*(m_coord + (i-1)) + (j-1)) <<' ';
		}
	}
}

void  CLbmGrid::Lbm_read (std::fstream& inputFile)
{
	inputFile >> m_nGrids >> m_nCoords >> m_NVertex;
	inputFile >> m_pgridName;
	inputFile >> m_gridIndex;
	int dataType = 0;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" # of vertices (CLbmGrid) " <<m_NVertex<<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	CreateGridMemory(m_nCoords);

	for(cgsize_t i=1; i <= m_nCoords; i++)
	{
		m_pcoordName[i-1]  = new char [5000];
		inputFile >> m_pcoordName[i-1];
		inputFile >> m_coordIndex[i-1] >> dataType;
	  	m_dataType[i-1]   = static_cast<DataType_t> (dataType);

		*(m_coord + (i-1)) = new double [m_NVertex];

		for(cgsize_t j=1; j <= m_NVertex; j++)
		{
			inputFile >> *(*(m_coord + (i-1)) + (j-1));
		}
	}
}

CLbmGrid::~CLbmGrid()
{
	SAFE_DELETE_ARRAY(m_pcoordName);
	SAFE_DELETE_ARRAY(m_coord);
	SAFE_DELETE_ARRAY(m_pgridName);
	SAFE_DELETE_ARRAY(m_dataType);
	SAFE_DELETE_ARRAY(m_coordIndex);
	std::cerr<<"CLbmGrid destructor called\n";
}

/////////////////////////////////////////////////////////////////////////////
CLbmProperties::CLbmProperties(char* _materialName , NodeValueType_t _density,	NodeValueType_t _k_viscosity,
			NodeValueType_t _Lc, NodeValueType_t _Re, NodeValueType_t _tau, NodeValueType_t _vel, NodeValueType_t _gx, 
			NodeValueType_t _gy, NodeValueType_t _gz)
{
	materialName = _materialName;
	density = _density;
	k_viscosity = _k_viscosity;
	Lc  = _Lc;
	Re  = _Re;
	tau = _tau;
	vel = _vel;
	gx = _gx;
	gy = _gy;
	gz = _gz;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmProperties \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

///////////////////////////////////////////////////////////////////////////////
CLbmElement::CLbmElement()
{
	m_nSections       = 0;
	m_elementType     = NULL;
	m_start           = NULL;
	m_end             = NULL;
	m_elementDataSize = NULL;
	m_nbndry          = NULL;
	m_parentFlag      = NULL;
	m_npe             = NULL;
	m_conn            = NULL;
	m_secIndex        = NULL;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmElement \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

void CLbmElement::CreateElementMemory(int size)
{
	m_conn            = new cgsize_t*      [size];
	m_elementDataSize = new cgsize_t       [size];
	m_elementType     = new ElementType_t  [size];
	m_end             = new cgsize_t       [size];
	m_nbndry          = new int            [size];
	m_npe             = new int            [size];
	m_parentFlag      = new int            [size];
	m_psectionName    = new char*          [size];
	m_start           = new cgsize_t       [size];
	m_secIndex        = new cgsize_t       [size];
}
void  CLbmElement::Lbm_write (std::fstream& outputFile)
{
	outputFile << m_nSections <<' ';
	for(cgsize_t i=1; i <= m_nSections; i++)
		{
			outputFile << m_psectionName[i-1] <<' ';
			outputFile << static_cast<int> ( m_elementType[i-1]) <<' '<< m_start[i-1] <<' '<< m_end[i-1]<<' '
				 << m_elementDataSize[i-1] <<' '<< m_nbndry[i-1] <<' '<< m_parentFlag[i-1] <<' '<< m_npe[i-1] <<' '<< m_secIndex[i-1]<<' ';

			for(cgsize_t j=1; j <= m_elementDataSize[i-1]; j++)
			{
				outputFile << *(*(m_conn + (i-1)) + (j-1))<<' ';
			}
		}
		//write the index of section(s) that contains the material(s)
		outputFile << m_materialIndex.size() <<' ';

		for(MaterialIndex_t::iterator iter = m_materialIndex.begin(); iter !=  m_materialIndex.end() ; iter++)
		{
			outputFile << *iter <<' ';
		}

		//write material(s) information to document
		outputFile << m_sectionsMaterial.size()<<' ';

		for(std::vector<cgsize_t>::const_iterator   iter1 = m_materialIndex.begin(); iter1 !=  m_materialIndex.end() ; iter1++)
		{
			std::map  <cgsize_t,MaterialRecord_t>::const_iterator iter2 = m_sectionsMaterial.find(*iter1);
			outputFile<< iter2->second.second.materialName<<' ';
			outputFile<< iter2->second.second.density<<' '
				<< iter2->second.second.k_viscosity<<' '
				<< iter2->second.second.Lc<<' '
				<< iter2->second.second.Re<<' '
				<< iter2->second.second.tau<<' '
				<< iter2->second.second.vel<<' '
				<< iter2->second.second.gx<<' '
				<< iter2->second.second.gy<<' '
				<< iter2->second.second.gz<<' ';
		}	
}

void  CLbmElement::Lbm_read (std::fstream& inputFile)
{
	inputFile >> m_nSections;
	int elementType = 0;
	CreateElementMemory(m_nSections);
	for(cgsize_t i=1; i <= m_nSections; i++)
	{
		m_psectionName[i-1] = new char [5000];
		inputFile >> m_psectionName[i-1];
		inputFile >> elementType >> m_start[i-1] >> m_end[i-1]
			>> m_elementDataSize[i-1] >> m_nbndry[i-1] >> m_parentFlag[i-1] >> m_npe[i-1] >> m_secIndex[i-1];

	  	m_elementType[i-1] = static_cast<ElementType_t> (elementType);

		*(m_conn + (i-1))  = new cgsize_t [*(m_elementDataSize + (i-1))];

		for(cgsize_t j=1; j <= m_elementDataSize[i-1]; j++)
		{
			inputFile >> *(*(m_conn + (i-1)) + (j-1));
		}
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished with sections \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//read the index of section(s) that contains the material(s)
	cgsize_t index(0);
	size_t   indexCount(0);
	inputFile >> indexCount;
	for(cgsize_t i=0; i < static_cast<cgsize_t>(indexCount); i++)
	{
		if(indexCount)
		{
			inputFile >> index;
			m_materialIndex.push_back(index);
		}
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished reading document sections that contain material\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//read material(s) information to document
	size_t materialCount(0);
	char* materialName = new char [5000];
	CLbmBoco::NodeValueType_t density, k_viscosity, Lc, Re, tau, vel, gx, gy, gz;
		
	inputFile >> materialCount;
	materialCount = 1;
	//FORCE MATERIAL COUNT TO ZERO JUST TO TEST SCMP
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" material count is ="<<materialCount<<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	for(cgsize_t i=0; i < static_cast<cgsize_t>(materialCount); i++)
	{
		if(materialCount)
		{
			inputFile >> materialName;
			inputFile >> density >> k_viscosity >> Lc >> Re >> tau >> vel >> gx >> gy >> gz;
			std::pair <char*, MaterialType_t>       material = MaterialRecord_t(materialName,MaterialType_t(materialName, density, k_viscosity, Lc, Re, tau, 
					vel, gx, gy, gz));
			std::pair <cgsize_t,MaterialRecord_t>  sectionmaterial = std::make_pair(m_materialIndex.at(i),material);
			m_sectionsMaterial.insert(sectionmaterial);
		}
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished reading document sections that contain material\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

CLbmElement::~CLbmElement()
{
	SAFE_DELETE_ARRAY(m_psectionName);
	SAFE_DELETE_ARRAY(m_conn);
	SAFE_DELETE_ARRAY(m_elementType);
	SAFE_DELETE_ARRAY(m_elementDataSize);
	SAFE_DELETE_ARRAY(m_start);
	SAFE_DELETE_ARRAY(m_end);
	SAFE_DELETE_ARRAY(m_nbndry);
	SAFE_DELETE_ARRAY(m_parentFlag);
	SAFE_DELETE_ARRAY(m_npe);
	SAFE_DELETE_ARRAY(m_secIndex);
	std::cerr<<"CLbmElement destructor called\n";
}


//////////////////////////////////////////////////////////////////////////////////
CLbmBoco::CLbmBoco()
{
	m_nBocos         = 0;
	m_bocoType       = NULL;
	m_ptsetType      = NULL;
	m_npnts          = NULL;
	m_normal_index   = NULL;
	m_normalListSize = NULL;
	m_normalDataType = NULL;
	m_ndataset       = NULL;
	m_pnts           = NULL;
	m_bocoXVel       = NULL;
	m_bocoYVel       = NULL;
	m_bocoZVel       = NULL;
	m_bocoPress      = NULL;
	m_bocoIndex      = NULL;
	m_bocoPressureBC = NULL;
	m_bocoVelocityBC = NULL;
	m_bocoExtrapolationBC = NULL;
	m_bocoWall       = NULL;
	m_bocoUserDef    = NULL;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmBoco \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

void CLbmBoco::CreateBocoMemory(int size)
{
	m_bocoType       = new BCType_t       [size];
	m_ptsetType      = new PointSetType_t [size];
	m_npnts          = new cgsize_t       [size];
	m_normal_index   = new cgsize_t       [size];
	m_normalListSize = new cgsize_t       [size];
	m_normalDataType = new DataType_t     [size];
	m_ndataset       = new cgsize_t       [size];
	m_pbocoName      = new char*          [size];
	m_pnts           = new cgsize_t*      [size];
	m_bocoXVel       = new NodeValueType_t[size];
	m_bocoYVel       = new NodeValueType_t[size];
	m_bocoZVel       = new NodeValueType_t[size];
	m_bocoPress      = new NodeValueType_t[size];
	m_bocoIndex      = new cgsize_t       [size];
	m_bocoPressureBC = new BOOL           [size];
	m_bocoVelocityBC = new BOOL           [size];
	m_bocoExtrapolationBC = new BOOL      [size];
	m_bocoWall       = new cgsize_t       [size];
	m_bocoUserDef    = new cgsize_t       [size];
}
void  CLbmBoco::Lbm_write (std::fstream& outputFile)
{
	outputFile <<m_nBocos<<' ';
	for(cgsize_t i=1; i <= m_nBocos; i++)
		{
			outputFile << *(m_pbocoName+(i-1))<<' ';
			outputFile << static_cast<int> (*(m_bocoType+(i-1))) <<' ';
			outputFile << static_cast<int> ( *(m_ptsetType+(i-1))) <<' ';
			outputFile <<*(m_npnts+(i-1)) <<' ';
			outputFile << *(m_normal_index+(i-1))<<' ';
			outputFile << static_cast<int> (*(m_normalDataType+(i-1))) <<' ';
			outputFile << *(m_ndataset+(i-1)) <<' ';
			outputFile<< *(m_bocoXVel + (i-1)) <<' ';
			outputFile << *(m_bocoYVel + (i-1)) <<' ';
			outputFile<<  *(m_bocoZVel + (i-1))<<' ';
			outputFile << *(m_bocoPress + (i-1)) <<' ';
			outputFile<< *(m_bocoIndex + (i-1)) <<' ';
			outputFile << *(m_bocoPressureBC + (i-1)) <<' ';
			outputFile<< *(m_bocoVelocityBC  +(i-1)) <<' ';
			outputFile << *(m_bocoExtrapolationBC +(i-1))<<' ';
			outputFile << *(m_bocoWall+(i-1)) <<' ';
			outputFile << *(m_bocoUserDef+(i-1)) <<' ';

			for(cgsize_t j=1; j <= *(m_npnts+(i-1)); j++)
			{
				outputFile << *(*(m_pnts + (i-1)) + (j-1))<<' ';
			}
		}
}

void  CLbmBoco::Lbm_read (std::fstream& inputFile)
{
	int bocoType       = 0;
	int ptSetType      = 0;
	int normalDataType = 0;
	inputFile >> m_nBocos;
	CreateBocoMemory(m_nBocos);
	for(cgsize_t i=1; i <= m_nBocos; i++)
		{
			m_pbocoName[i-1] = new char [5000];
			inputFile >> m_pbocoName[i-1];
			inputFile >> bocoType ;
			inputFile >> ptSetType;
			inputFile >> m_npnts[i-1];
			inputFile >> m_normal_index[i-1];
			inputFile >> normalDataType;
			inputFile >>  m_ndataset[i-1];
			inputFile >> m_bocoXVel[i-1];
			inputFile >> m_bocoYVel[i-1];
			inputFile >> m_bocoZVel[i-1];
			inputFile >> m_bocoPress[i-1];
			inputFile >> m_bocoIndex[i-1];
			inputFile >> m_bocoPressureBC[i-1];
			inputFile >> m_bocoVelocityBC[i-1];
			inputFile >> m_bocoExtrapolationBC[i-1];
			inputFile >> m_bocoWall[i-1];
			inputFile >> m_bocoUserDef[i-1];

			m_bocoType[i-1]       = static_cast<BCType_t>       (bocoType);
			m_ptsetType[i-1]      = static_cast<PointSetType_t> (ptSetType);
			m_normalDataType[i-1] = static_cast<DataType_t>     (normalDataType);
			*(m_pnts + (i-1)) = new cgsize_t [*(m_npnts + (i-1))];

			for(cgsize_t j=1; j <= m_npnts[i-1]; j++)
			{
				inputFile >> *(*(m_pnts + (i-1)) + (j-1));
			}
		}
}

CLbmBoco::~CLbmBoco()
{
	SAFE_DELETE_ARRAY(m_pbocoName);
	SAFE_DELETE_ARRAY(m_pnts);
	SAFE_DELETE_ARRAY(m_bocoType);
	SAFE_DELETE_ARRAY(m_ptsetType);
	SAFE_DELETE_ARRAY(m_npnts);
	SAFE_DELETE_ARRAY(m_normal_index);
	SAFE_DELETE_ARRAY(m_normalListSize);
	SAFE_DELETE_ARRAY(m_normalDataType);
	SAFE_DELETE_ARRAY(m_bocoXVel);
	SAFE_DELETE_ARRAY(m_bocoYVel);
	SAFE_DELETE_ARRAY(m_bocoZVel);
	SAFE_DELETE_ARRAY(m_bocoPress);
	SAFE_DELETE_ARRAY(m_bocoIndex);
	SAFE_DELETE_ARRAY(m_bocoPressureBC);
	SAFE_DELETE_ARRAY(m_bocoVelocityBC);
	SAFE_DELETE_ARRAY(m_bocoExtrapolationBC);
	SAFE_DELETE_ARRAY(m_bocoWall);
	SAFE_DELETE_ARRAY(m_bocoUserDef);
	std::cerr<<"CLbmBoco destructor called\n";
}

//////////////////////////////////////////////////////////////////////////
CLbmSol::CLbmSol() 
{
	char solName[5000] = "FlowSolution";
	m_psolName = solName;
	//m_nfield = 11;
	m_nfield = 9;
	m_nsol = 1;
	m_solIndex = 0;
	m_NVertex = 0;
	m_iteration = 0;
	m_loc=Vertex;
 	m_datatype = RealDouble;
	m_pfieldName = new char*            [m_nfield];
	m_fieldIndex = new cgsize_t         [m_nfield];

	char xvelocity[100]    = "X-Velocity";
	char yvelocity[100]    = "Y-Velocity";
	char zvelocity[100]    = "Z-Velocity";
	char velMagnitude[100] = "Velocity-Magnitude";
	char density[100]      = "Density";
	char uoverallx[100]    = "Uoverall_X";
	char uoverally[100]    = "Uoverall_Y";
	char uoverallz[100]    = "Uoverall_Z";
	char pressure[100]     = "Pressure";
	/*char uprimex[100]      = "Uprime_X";
	char uprimey[100]      = "Uprime_Y";
	char uprimez[100]      = "Uprime_Z";*/



	m_pfieldName[0]  = xvelocity;
	m_pfieldName[1]  = yvelocity;
	m_pfieldName[2]  = zvelocity;
	m_pfieldName[3]  = velMagnitude;
	m_pfieldName[4]  = density;
	m_pfieldName[5]  = uoverallx;
	m_pfieldName[6]  = uoverally;
	m_pfieldName[7]  = uoverallz;
	m_pfieldName[7]  = pressure;
	/*m_pfieldName[8]  = uprimex;
	m_pfieldName[9]  = uprimey;
	m_pfieldName[10] = uprimez;*/
	
	for(cgsize_t i=1; i <= m_nfield; i++)
	{
		m_fieldIndex[i-1] = i;
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating CLbmSol \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void CLbmSol::CreateSolMemory(int size)
{
	cgsize_t multi(1);

	if(m_phaseModel == CLbmModel::MULTIPHASE)
		multi=2;
	else;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" multi = "<<multi <<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" # of vertices (CLbmSol) " <<m_NVertex<<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	for(cgsize_t i=1; i<=m_nfield; i++)
	{
		m_sol.push_back(new FieldType_t);
		m_solOld.push_back(new FieldType_t);
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" Field in domain " <<m_sol.size()<<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

		for(cgsize_t j=1; j<=multi; j++)
		{
			m_sol.at(i-1)->push_back(new MaterialFieldType_t);
			m_solOld.at(i-1)->push_back(new MaterialFieldType_t);
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" each field container size " <<m_sol.at(i-1)->size()<<"\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		
			for(cgsize_t k=1; k <=m_NVertex; k++)
			{
				m_sol.at(i-1)->at(j-1)->push_back(0.0);
				m_solOld.at(i-1)->at(j-1)->push_back(0.0);
			}
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" material container size " <<m_sol.at(i-1)->at(j-1)->size()<<"\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			
		}
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished creating memory for solution with size " <<m_sol.size()<<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

void  CLbmSol::Lbm_write (std::fstream& outputFile)
{
	cgsize_t multi(1);

	outputFile <<static_cast<int>(m_datatype) <<' '<< static_cast<int>(m_loc) <<' '<< static_cast<int>(m_phaseModel) 
		<<' '<< m_nsol <<' '<< m_solIndex <<' '<< m_NVertex <<' '<< m_nfield <<' '<< m_iteration <<' ';

	outputFile <<m_convergenceData.x_residual <<' '<< m_convergenceData.y_residual <<' '<< m_convergenceData.z_residual<<' ';

	if(m_phaseModel == CLbmModel::MULTIPHASE)
		multi=2;
	else;
	
	for(cgsize_t i=0; i<m_nfield; i++)
		{
			for(cgsize_t j=0; j< multi; j++)
			{
				for(cgsize_t k=0; k< m_NVertex; k++)
				{
					outputFile << m_sol.at(i)->at(j)->at(k)<<' ';
					outputFile << m_solOld.at(i)->at(j)->at(k)<<' ';
					/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
					if(k == 79999)
					{
					std::cout<<" inside Lbm_Write \n";
					std::cout<<" field = "<<i<<" node# = "<<k<<" value = "<< m_sol.at(i)->at(j)->at(k) <<"\n";
					
					std::cin.get();
					}
					else;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
				}
			}
		}
}

void  CLbmSol::Lbm_read (std::fstream& inputFile)
{
	cgsize_t multi(1);
	int datatype = 0;
	int location = 0;
	int phasemodel = 0;
	inputFile >> datatype >> location >> phasemodel >> m_nsol >> m_solIndex >> m_NVertex >> m_nfield >> m_iteration ;

	inputFile >> m_convergenceData.x_residual >> m_convergenceData.y_residual >> m_convergenceData.z_residual;

		m_datatype   = static_cast<DataType_t> (datatype);
		m_loc        = static_cast<GridLocation_t> (location);
		m_phaseModel = static_cast<CLbmModel::PHASEMODEL> (phasemodel);
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" m_phaseModel = "<<m_phaseModel <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		if(m_phaseModel == CLbmModel::MULTIPHASE)
		{
			multi = 2;
		}
		//FORCE FIELD TO 9 ELEMENTS
		m_nfield = 9;
		CreateSolMemory(m_nfield);


		for(cgsize_t i=0; i< m_nfield; i++)
		{
			for(cgsize_t j=0; j< multi; j++)
			{
				/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::cout<<" size for fields container = "<<m_sol.size() <<"for field = "<<j<<"\n";
				std::cin.get();
				////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
				for(cgsize_t k=0; k< m_NVertex; k++)
				{
					//m_sol.at(i-1)->at(j-1)->at(k-1) = 0.0;
					/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
					std::cout<<" size for nodes container = "<<m_sol.at(i-1)->at(j-1)->size() <<"\n";
					std::cin.get();
					////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
					inputFile >> m_sol.at(i)->at(j)->at(k);
					inputFile >> m_solOld.at(i)->at(j)->at(k);
					/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
					if(k == 79999)
					{
					std::cout<<" inside Lbm_Read \n";
					std::cout<<" field = "<<i<<" node# = "<<k<<" value = "<<m_sol.at(i)->at(j)->at(k) <<"\n";
					
					std::cin.get();
					}
					else;
					////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

				}
			}
		}
}

CLbmSol::~CLbmSol()
{
	SAFE_DELETE_ARRAY(m_pfieldName);
	SAFE_DELETE_ARRAY(m_fieldIndex);
	std::cerr<<"CLbmSol destructor called\n";
}

/////////////////////////////////////////////////////////////////////////////
CMultiPhasePatch::CMultiPhasePatch()
{
	x_min = 0.0;
	y_min = 0.0;
	z_min = 0.0;
	x_max = 0.0;
	y_max = 0.0;
	z_max = 0.0;
}
void  CMultiPhasePatch::Lbm_write (std::fstream& outputFile)
{
	outputFile << x_min <<' '<< y_min <<' '<< z_min <<' '<< x_max <<' '<< y_max <<' '<< z_max <<' ';
}

void  CMultiPhasePatch::Lbm_read (std::fstream& inputFile)
{
	inputFile >> x_min >> y_min >>  z_min >> x_max >> y_max >> z_max;
}

CMultiPhasePatch::~CMultiPhasePatch()
{
	std::cerr<<"CMultiPhasePatch destructor called\n";
}

Convergence::Convergence()
{
	x_residual       = 1.0;
	y_residual       = 1.0;
	z_residual       = 1.0;
	density_residual = 1.0;
	temp_residual    = 1.0;
}

