#pragma once
#include <cstdlib>
#include <vector>
#include <map>
#include <fstream>
#include <string>
//#include "CLbmGlobal.h"

////////////////////////////////////
#define SAFE_DELETE(ptr) if (ptr) { delete ptr; ptr = NULL; }
#define SAFE_DELETE_ARRAY(ptr) if(ptr) {delete [] ptr ;ptr = NULL ;}

typedef int  cgsize_t;
typedef bool BOOL;
////////////////////////////////////////////////////////////////////
enum ZoneType_t
{
	Stuctured    =2,
	Unstructured =3
} ;

enum DataType_t
{
	Integer    =2,
	RealSingle  =3,
	RealDouble  =4,
	Character   =5, 
	LongInteger =6
};

enum ElementType_t
{ 
	NODE     =2,
  BAR_2    =3,
	BAR_3    =4,
	TRI_3    =5,
	TRI_6    =6,
	QUAD_4   =7,
	QUAD_8   =8,
	QUAD_9   =9,
	TETRA_4  =10,
	TETRA_10 =11,
	PYRA_5   =12,
	PYRA_14  =13,
	PENTA_6  =14,
	PENTA_15 =15,
	PENTA_18 =16,
	HEXA_8   =17,
	HEXA_20  =18, 
	HEXA_27  =19, 
	MIXED    =20, 
	PYRA_13  =21, 
	NGON_n   =22, 
	NFACE_n  =23
};

enum BCType_t
{
	//BCAxisymmetricWedge =2,
	//BCDegenerateLine =3,
	//BCDegeneratePoint =4,
	BCDirichlet =5,
	//BCExtrapolate =6,
	//BCFarfield =7,
	BCGeneral =8,
	//BCInflow =9,
	//BCInflowSubsonic =10,
	//BCInflowSupersonic =11,
	//BCNeumann =12,
	BCOutflow =13,
	//BCOutflowSubsonic =14,
	//BCOutflowSupersonic =15,
	//BCSymmetryPlane =16,
	//BCSymmetryPolar =17,
	//BCTunnelInflow =18,
	//BCTunnelOutflow =19,
	BCWall =20,
	//BCWallInviscid =21,
	//BCWallViscous =22,
	//BCWallViscousHeatFlux =23,
	//BCWallViscousIsothermal =24, 
	//FamilySpecified =25
};

enum PointSetType_t
{
	PointList =2,
	PointListDonor =3,
	PointRange =4,
	PointRangeDonor =5,
	ElementRange =6,
	ElementList =7,
	CellListDonor =8
};

enum GridLocation_t
{
	Vertex =2,
	CellCenter =3,
	FaceCenter =4,
	IFaceCenter =5,
	JFaceCenter =6, 
	KFaceCenter =7, 
	EdgeCenter =8
};

/////////////////////////////////////////

class Domain;
class Convergence
{
public:
	typedef double NodeValueType_t;
	Convergence();
	NodeValueType_t x_residual;
	NodeValueType_t y_residual;
	NodeValueType_t z_residual;
	NodeValueType_t density_residual;
	NodeValueType_t temp_residual;
};

class CLbmBase 
{
public:
	CLbmBase();
	virtual ~CLbmBase();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);
	
	cgsize_t  index_file;
	cgsize_t  index_base;
	cgsize_t  nBases;
	char*         m_pbaseName;
	cgsize_t  cell_dim;
	cgsize_t  phys_dim;
};

class CLbmModel
{
public:
	CLbmModel();
	virtual ~CLbmModel();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);

	enum SOLVER      {D2Q7, D2Q9, D3Q15, D3Q19};
	enum PHASEMODEL  {OFF,  MULTIPHASE};
	enum RELAXTIME   {SRT, TRT, MRT};
	SOLVER       m_solver;
	PHASEMODEL   m_phaseModel;
	RELAXTIME    m_relaxTime;
	BOOL         m_energyEqn;
	BOOL         m_gravity;
	char       m_initStatus;
};

class CLbmZone
{
public:
	CLbmZone();
	virtual ~CLbmZone();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);
	
	cgsize_t   m_nZones;
	ZoneType_t m_zoneType;
	char*      m_pzoneName; 
	cgsize_t*  m_size;
	cgsize_t   m_index_dim;
	cgsize_t   m_index_zone;
};

class CLbmGrid
{
public:
	CLbmGrid();
	virtual ~CLbmGrid();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);
	
	int         m_nGrids;
	char*       m_pgridName;
	int         m_nCoords;
	char**      m_pcoordName;
	DataType_t* m_dataType;
	double**    m_coord;
	cgsize_t    m_NVertex;
	cgsize_t*   m_coordIndex;
	cgsize_t    m_gridIndex;
	void CreateGridMemory(int size);
};
/////////////////////////////////////////////////////////////////
class CLbmProperties
	{
	public:
		typedef double NodeValueType_t;
		CLbmProperties(char* _materialName = NULL, NodeValueType_t _density = 0.0,	NodeValueType_t _k_viscosity = 0.0,
			NodeValueType_t _Lc = 0.0, NodeValueType_t _Re = 0.0, NodeValueType_t _tau = 0.0, NodeValueType_t _vel= 0.0,
			NodeValueType_t _gx = 0.0, NodeValueType_t _gy = 0.0, NodeValueType_t _gz = 0.0 );
		std::string materialName;
		NodeValueType_t density;
		NodeValueType_t k_viscosity;
		NodeValueType_t Lc;
		NodeValueType_t Re;
		NodeValueType_t tau;
		NodeValueType_t vel;
		NodeValueType_t gx;
		NodeValueType_t gy;
		NodeValueType_t gz;
	};

////////////////////////////////////////////////////////////////
class CLbmElement
{
public:
	typedef   CLbmProperties MaterialType_t;
	typedef std::vector<cgsize_t>               MaterialIndex_t;
	typedef std::map  <char*, MaterialType_t>   MaterialTypes_t;
	typedef std::pair <char*, MaterialType_t>   MaterialRecord_t;
	typedef std::map  <cgsize_t,MaterialRecord_t> SectionsMaterial_t;
	typedef std::pair <cgsize_t,MaterialRecord_t> SectionMaterialRecord_t;
	
	CLbmElement();
	virtual ~CLbmElement();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);
	
	int                m_nSections;
	char**             m_psectionName;
	ElementType_t*     m_elementType;
	cgsize_t*          m_start;
	cgsize_t*          m_end;
	cgsize_t*          m_elementDataSize;
	cgsize_t**         m_conn;
	int*               m_nbndry;
	int*               m_parentFlag;
	cgsize_t*          m_npe;
	cgsize_t*          m_secIndex;
	MaterialIndex_t    m_materialIndex;
	SectionsMaterial_t m_sectionsMaterial;
	void CreateElementMemory(int size);
};

class CLbmBoco
{
public:
	typedef double NodeValueType_t;
	
	CLbmBoco();
	virtual ~CLbmBoco();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);
	
	cgsize_t         m_nBocos;
	char**           m_pbocoName;
	BCType_t*        m_bocoType;
	PointSetType_t*  m_ptsetType;
	cgsize_t*        m_npnts;
	cgsize_t*        m_normal_index;
	cgsize_t*        m_normalListSize;
	DataType_t*      m_normalDataType;
	cgsize_t*        m_ndataset;
	cgsize_t**       m_pnts;
	NodeValueType_t* m_bocoXVel;
	NodeValueType_t* m_bocoYVel;
	NodeValueType_t* m_bocoZVel;
	NodeValueType_t* m_bocoPress;
	cgsize_t*        m_bocoIndex;
	BOOL*            m_bocoPressureBC;
	BOOL*            m_bocoVelocityBC;
	BOOL*            m_bocoExtrapolationBC;
	cgsize_t*        m_bocoWall;
	cgsize_t*        m_bocoUserDef;
	void CreateBocoMemory(int size);
};

class CLbmSol
{
public:
	typedef double NodeValueType_t;
	typedef std::vector<NodeValueType_t>          MaterialFieldType_t;
	typedef std::vector<MaterialFieldType_t*>      FieldType_t;
	typedef std::vector<FieldType_t*>              SolutionType_t;


	CLbmSol();
	virtual ~CLbmSol();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);
	char* 		m_psolName;
	char**		m_pfieldName;
	cgsize_t	m_nsol;
	cgsize_t	m_nfield;
	cgsize_t	m_solIndex;
	cgsize_t	m_NVertex;
	cgsize_t*	m_fieldIndex;
	GridLocation_t	m_loc;
	DataType_t	m_datatype;
	SolutionType_t  m_sol;
	SolutionType_t  m_solOld;
	Convergence	m_convergenceData;
	cgsize_t	m_iteration;
	CLbmModel::PHASEMODEL   m_phaseModel;
	void CreateSolMemory(int size);
};

struct CMultiPhasePatch 
{
	typedef double NodeValueType_t;

	CMultiPhasePatch();
	virtual ~CMultiPhasePatch();
	void  Lbm_write (std::fstream& outputFile);
	void  Lbm_read  (std::fstream& inputFile);

	NodeValueType_t x_min;
	NodeValueType_t y_min;
	NodeValueType_t z_min;
	NodeValueType_t x_max;
	NodeValueType_t y_max;
	NodeValueType_t z_max;
};

class CLbmCase 
{
public:
	CLbmCase(void);
	virtual ~CLbmCase(void);

	void ExportCaseInfo(char* transferName);
	void ImportCaseInfo(char* transferName);
	void Lbm_write_tecplot(char* tecplotFileName);
		
	//struct members
	CLbmBase      	m_base;
	CLbmModel     	m_model;
	CLbmZone     	m_zone;
	CLbmGrid      	m_grid;
	CLbmElement	m_element;
	CLbmBoco	m_boco;
	CLbmSol		m_sol;
	CMultiPhasePatch m_patchedNodes;
  	Domain*       	m_pDomain;//put in private later
	cgsize_t	m_maxIteration; //put in private later

	std::fstream m_importfs;
	std::fstream m_exportfs;
	std::fstream m_exporttecplot;
private:
	//tecplot utility functions
	void setFileHeader(std::fstream& outputFile, char* tecplotFileName);
	void setFileHeaderMulti(std::fstream& outputFile, char* tecplotFileName);
	void setZoneHeader(std::fstream& outputFile);
	void setCoordData(std::fstream& outputFile);
	void setVarData(std::fstream& outputFile);
	void setVarDataMulti(std::fstream& outputFile);
	void setConnData(std::fstream& outputFile);
};


