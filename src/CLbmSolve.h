#pragma once

#ifndef  CLBMSOLVE_H
#define  CLBMSOLVE_H

#include <cstdlib>
#include "CLbmGlobal.h"
#include "CLbmCase.h"
#include "CLbmBlock.h"
#include "CLbmBlock3D.h"
#include "CVariable.h"
#include <vector>
#include <fstream>

//LEVEL ONE:invoked in OnRun()
//==============================================================================================//
class PdfDomain;


class Domain //base class of will implement D2Q9, D3Q15, ... related objects
{
public:
	typedef std::vector<PdfDomain*>       MaterialType_t;
	typedef std::vector<CVariable*>       VariableType_t;
	typedef std::vector<cgsize_t>         PatchNodeType_t;
	typedef std::vector<cgsize_t>         PillarNodeType_t;
	typedef std::vector<PillarNodeType_t*>    MicroPillarsType_t;

	//Getters
	MaterialType_t*    GetDomainMaterial();
	VariableType_t*    GetDomainVariables();
	VariableType_t*    GetDomainTempVariables();
	PatchNodeType_t*   GetDomainSolidObjectNodes();
	MicroPillarsType_t* GetDomainMicroPillars();

	//initializing methods
	void mapBounds(CLbmCase* pCase);
	void mapBoundType(CLbmCase* pCase);
	void mapPeriodicNeighbor(CLbmCase* pCase);
	virtual void mapStreamingStart(CLbmCase* pCase)  = 0;
	virtual void mapCardinalDirection(cgsize_t node) = 0;
	void copyToOthers(cgsize_t node, Node::LBMBOUND bound);
	void copyToOthers(cgsize_t node, cgsize_t mat_index);
	void copyToOthers(cgsize_t node, cgsize_t mat_index, Node::NodeValueType_t density, Node::NodeValueType_t velocity);
	void initializeVars(CLbmCase*  pCase,  VariableType_t* variables);
	void initializeTemps(CLbmCase* pCase,  VariableType_t* tempVariables);
	void initialize(CLbmCase* pCase, Domain& domainVariables);
	void setPatchedSolidObjects(CLbmCase* pCase);
	void setPatchedMicroPillars(CLbmCase* pCase);
	
	
	//process methods
	void collide(CLbmCase* pCase, Domain& domainVariables);
	void post_collide(CLbmCase* pCase, Domain& domainVariables);
	void stream (CLbmCase* pCase, Domain& domainVariables);
	void post_stream (CLbmCase* pCase, Domain& domainVariables);

	void updateVariables(CLbmCase* pCase, Domain& domainVariables);
	void updateVariablesSglComp(CLbmCase* pCase, Domain& domainVariables);

	void writeSolution(CLbmCase* pCase);
	void Update_variables(CLbmCase* pCase);
	void Lbm_write(std::fstream& outputFile);
	void Lbm_read (std::fstream& inputFile);

	virtual void convergence(Domain& domainVariables, CLbmCase* pCase) = 0;
	virtual Convergence*    GetConvergence()          = 0;
	//parallel implementation methods
	virtual void SetNumThreads(cgsize_t  num)         = 0;
	virtual cgsize_t GetNumThreads() const            = 0;
	//dtor
	virtual ~Domain(){};
private:
	//Getters utility
	virtual MaterialType_t*  DomainMaterial()         = 0;
	virtual VariableType_t*  DomainVariables()        = 0;
	virtual VariableType_t*  DomainTempVariables()    = 0;
	virtual PatchNodeType_t* DomainSolidObjectNodes() = 0;
	virtual MicroPillarsType_t* DomainMicroPillars() = 0;
	//iniitalizing utility
	void mapSolidNodes(CLbmCase* pCase, cgsize_t node, cgsize_t bry);
	void mapExtrapNodes(CLbmCase* pCase, cgsize_t node, cgsize_t bry);
	void mapInwardORBuriedNodes(cgsize_t node);
	void copyPeriodicNeighbor(CLbmCase* pCase, cgsize_t node, cgsize_t neighbor);
	void mapPeriodicCorners(CLbmCase* pCase, cgsize_t mat_index, cgsize_t node);
	void mapPeriodicEdges(CLbmCase* pCase, cgsize_t mat_index, cgsize_t node);
	void mapPeriodicSurfaces(CLbmCase* pCase, cgsize_t mat_index, cgsize_t node);
	cgsize_t fetchPeriodicNeighbor(CLbmCase* pCase, Node::LBMBOUND bound, cgsize_t mat_index, cgsize_t periodicNode);
	cgsize_t fetchPeriodicNeighbor(CLbmCase* pCase, Node::LBMBOUND bound, cgsize_t mat_index, cgsize_t periodicNode, Node::NodeValueType_t coord1, cgsize_t index);
	cgsize_t fetchPeriodicNeighbor(CLbmCase* pCase, Node::LBMBOUND bound, cgsize_t mat_index, cgsize_t periodicNode, 
		Node::NodeValueType_t coord1, Node::NodeValueType_t coord2, cgsize_t index1, cgsize_t index2);

	//solid object initialization utility
	void updateSolidNodes(cgsize_t node);
	cgsize_t fetchPrismNode(CLbmCase* pCase, cgsize_t node);
	cgsize_t fetchSphereNode(CLbmCase* pCase, cgsize_t node);
	cgsize_t fetchCylinderNode(CLbmCase* pCase, cgsize_t node);

	//micropillar utility
	void updatePillarNodesWetting(cgsize_t node, PillarNodeType_t* pillar);
	void scanPillar(Node::NodeValueType_t coord1, Node::NodeValueType_t coord2, PillarNodeType_t* pillar);
	bool nodeInRange(cgsize_t node, Node::NodeValueType_t coord1Min, Node::NodeValueType_t coord2Min, Node::NodeValueType_t coord3Min, 
		    	 Node::NodeValueType_t coord1Max, Node::NodeValueType_t coord2Max, Node::NodeValueType_t coord3Max);
	

	//process utility
	//virtual void computeCompositeVars(CLbmCase* pCase, Domain& domainVariables)  = 0;
	virtual void computeCompositeVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)=0;
	//virtual void computeOverallVars(CLbmCase* pCase, Domain& domainVariables)    = 0;
	virtual void computeOverallVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)=0;
	virtual void computeOverallVarsSgl(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)=0;
	//virtual void computeIndividualVars(CLbmCase* pCase, Domain& domainVariables) = 0;
	virtual void computeIndividualVars(CLbmCase* pCase, Domain& domainVariables, MaterialType_t::size_type mat_index, cgsize_t node)=0;

	void computeInteractionForces(CLbmCase* pCase, Domain& domainVariables);
	void computeInteractionForcesSCMP(CLbmCase* pCase, Domain& domainVariables);
	void computeSurfaceForces(CLbmCase* pCase, Domain& domainVariables);

	//////////////// utility for computeInteractionForces
	void computePsi(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G);
	void computePhi(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub);
	Node::NodeValueType_t computePhi(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G);
	Node::NodeValueType_t wall_VirtualDensity(Domain& domainVariables, cgsize_t node, cgsize_t sub);
	Node::NodeValueType_t pressure_EOS_PR(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t density);
	Node::NodeValueType_t pressure_EOS_CS(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t density);
	Node::NodeValueType_t pressure_EOS_vdW(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t density);
	//void computeSummation(CLbmCase* pCase, Domain& domainVariables);
	void computeSummation(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G);
	void computeSummationIntra(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G, cgsize_t pseudoP);
	void computeSummationInter(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub1, cgsize_t sub2, Node::NodeValueType_t G, cgsize_t pseudoP);
	cgsize_t GetInteractionNeighbor(CLbmCase* pCase, cgsize_t f, cgsize_t node, cgsize_t sub);
	void D3Q15_FCohesionAdhesion(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t neighborNode, cgsize_t sub1, cgsize_t sub2, Node::NodeValueType_t G, cgsize_t pseudoP);
	void D3Q19_FCohesionAdhesion(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t neighborNode, cgsize_t sub1, cgsize_t sub2, Node::NodeValueType_t G, cgsize_t pseudoP);
	
	void computeSurfaceSummation(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub, Node::NodeValueType_t Gads);
	//void computeAllTerms(CLbmCase* pCase, Domain& domainVariables);
	void computeAllTerms(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	void computeSurfaceAllTerms(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	cgsize_t switchPdf(CLbmCase* pCase, cgsize_t f);
};

//Derived class that will implement D2Q9 related objects
class Domain_D2Q9: public Domain
{
public:
	Domain_D2Q9(CLbmCase* pCase);
	Domain_D2Q9(CLbmCase* pCase, cgsize_t nvertex);
	virtual void convergence(Domain& domainVariables, CLbmCase* pCase);
	virtual MaterialType_t*  DomainMaterial();
	virtual VariableType_t*  DomainVariables();
	virtual VariableType_t*  DomainTempVariables();
	virtual PatchNodeType_t* DomainSolidObjectNodes();
	virtual MicroPillarsType_t* DomainMicroPillars();

	//initializing methods
	//virtual void mapBounds(CLbmCase* pCase){};
	virtual void mapStreamingStart(CLbmCase* pCase){};
	virtual void mapCardinalDirection(cgsize_t node){};
	//virtual void computeCompositeVars(CLbmCase* pCase, Domain& domainVariables);
	virtual void computeCompositeVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	//virtual void computeOverallVars(CLbmCase* pCase, Domain& domainVariables);
	virtual void computeOverallVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	virtual void computeOverallVarsSgl(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	//virtual void computeIndividualVars(CLbmCase* pCase, Domain& domainVariables);
	virtual void computeIndividualVars(CLbmCase* pCase, Domain& domainVariables, MaterialType_t::size_type mat_index, cgsize_t node);
	virtual Convergence*    GetConvergence();
	virtual void SetNumThreads(cgsize_t  num);
	virtual cgsize_t GetNumThreads() const;
	virtual ~Domain_D2Q9();
private:
	MaterialType_t  m_domainMaterial;
	VariableType_t  m_domainVariables;
	VariableType_t  m_domainTempVariables;
	PatchNodeType_t m_domainSolidObjectNodes;
	MicroPillarsType_t m_domainMicroPillars;
	Convergence     m_convergenceData;
	cgsize_t        m_numThreads;
};

//LEVEL TWO:invoked from Domain
//==============================================================================================//
class LbmDomain;

class PdfDomain //base class of classes that will implement incompressible, compressible, solute, temperature Pdf objects
								// for any lattice type e.g D2Q9, D3Q15, ...
{
public:

	typedef std::vector<LbmDomain*> LatticeType_t;
	LatticeType_t* GetLatticePdf();
	cgsize_t       GetMaterialIndex();
	//virtual void updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index) = 0;
	virtual void updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t node, PdfDomain::LatticeType_t::size_type f)=0;
	virtual void momentToVelocitySpace(CLbmCase* pCase, Domain& domainVariables)=0;
	void initialize(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	void collide(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	void post_collide(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	void stream (CLbmCase* pCase, Domain& domainVariables);
	void post_stream (CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	void bound  (CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t bdry, BCType_t boundType, cgsize_t* pPnts , cgsize_t  numPnts);
	void Lbm_write(std::fstream& outputFile);
	void Lbm_read (std::fstream& inputFile);
	//density at boundaries (if velocity/pressure boundary conditions)
	Node::NodeValueType_t density(Domain& domainVariables, cgsize_t dNqN, cgsize_t node, Node::NodeValueType_t c, cgsize_t vel_index,
					cgsize_t fI=0,  cgsize_t fII=0,  cgsize_t fIII=0,	cgsize_t ffI=0, cgsize_t ffII=0, cgsize_t ffIII=0,
					cgsize_t fIV=0, cgsize_t fV=0, cgsize_t ffIV=0, cgsize_t ffV=0 );

	Node::NodeValueType_t velocity(Domain& domainVariables, cgsize_t dNqN, cgsize_t node, Node::NodeValueType_t c1, Node::NodeValueType_t c2,
					cgsize_t fI=0,  cgsize_t fII=0,  cgsize_t fIII=0,	cgsize_t ffI=0, cgsize_t ffII=0, cgsize_t ffIII=0,
					cgsize_t fIV=0, cgsize_t fV=0, cgsize_t ffIV=0, cgsize_t ffV=0 );
	
	virtual ~PdfDomain(){};

private:
	virtual LatticeType_t* LatticePdf()    = 0;
	virtual cgsize_t       MaterialIndex() = 0;
	//private method for bound method
	void bound_dirichlet(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t* pPnts, cgsize_t  numPnts, cgsize_t i, BCType_t boundType);
	void bound_open_bc(Domain& domainVariables, cgsize_t mat_index, cgsize_t* pPnts, cgsize_t  numPnts);
	void bound_wall_bc(Domain& domainVariables, CLbmCase* pCase, cgsize_t mat_index, cgsize_t* pPnts, cgsize_t  numPnts, BCType_t boundType, cgsize_t i);
	void bound_general_bc(Domain& domainVariables, CLbmCase* pCase, cgsize_t mat_index, cgsize_t* pPnts, cgsize_t  numPnts, cgsize_t i);
	//private method for stream method
	void stream_utility(Domain& domainVariables, cgsize_t* pPnts , cgsize_t  numPnts);
};

//Derived class that will implement D2Q9 Incompressible related Pdf objects
class D2Q9Incomp_domain: public PdfDomain
{
public:
	D2Q9Incomp_domain(CLbmCase* pCase, cgsize_t matIndex);
	D2Q9Incomp_domain(CLbmCase* pCase, cgsize_t matIndex, cgsize_t nvertex);
	virtual LatticeType_t* LatticePdf();
	virtual cgsize_t       MaterialIndex();
	//virtual void updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	virtual void updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t node, PdfDomain::LatticeType_t::size_type f);
	virtual void momentToVelocitySpace(CLbmCase* pCase, Domain& domainVariables){/*place holder*/};

	~D2Q9Incomp_domain();
private:
	LatticeType_t m_latticeType;
	cgsize_t      m_materialIndex;
};

//LEVEL THREE: invoked from 
//==============================================================================================//
/*----------------------------------------------------------------------------------------------*/
//LbmDomain Abstract class
//Properties:
//1. Base class of the model-specific processes implementation class heirachy
/*--------------------------------------------------------------------------------------------*/
class Collider;
class Streamer;
class Pdf_bounder;
class Pdf_setNeighbor;
//class Check_node;
class Srt_collider;
class Srt_collider_improved;
class Trt_collider;
class Mrt_collider;

//Base class of classes that will implement the specific distribution function f1, f2, f3, ....
class LbmDomain
{
public:
	typedef std::vector<cgsize_t>    PatchNodeType_t;
	//API for 
	//Enumeration for different boundary conditions
	enum   BOUND    {BOUNCEBACK, VELOCITY, SYMMETRY, PRESSURE, PERIODIC, EXTRAPOLATION};
	//Enumeration for different formulation for the different boundary conditions
	enum   FORMULA  {ZHOU_HE, ZHOU_HE_I, MOVING, OTHERS, FULL_BB, HALF_BB, EXTRAPOLATE, NBC, CBC_LV};
	static Srt_collider f_srt;
	static Srt_collider_improved f_srt_i;
	static Trt_collider f_trt;
	static Mrt_collider f_mrt;

	void lbmInitialize(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	void lbmCollide(PdfDomain* allPdfs, Domain& domainVariables, CLbmModel::RELAXTIME input, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	void lbmStream (cgsize_t node);
	void lbmWrite  (std::fstream& outputFile);
	void lbmRead   (std::fstream& inputFile);

	void openBCBounder(Domain& domainVariables, cgsize_t index, cgsize_t* pPnts, cgsize_t numPnts);
	void periodicBCBounder(Domain& domainVariables, cgsize_t index, cgsize_t* pPnts, cgsize_t numPnts);

	void solidObjectNodeBounder(Domain& domainVariables, cgsize_t index, PatchNodeType_t* solidObjectNodes, LbmDomain::BOUND b_type, FORMULA f_type);
	void nodeBounder  (Domain& domainVariables, cgsize_t mat_index, cgsize_t* pPnts, cgsize_t numPnts, 
									 BCType_t boundType, BOUND b_type, FORMULA f_type);
	PdfBlock*    pdf();
	PdfMaker*    Impl();

	virtual Pdf_bounder& GetVelBounder()      = 0;
	virtual Pdf_bounder& GetFullBBBounder()   = 0;
	virtual Pdf_bounder& GetPeriodicBounder() = 0;
	virtual Pdf_bounder& GetOpenBounder()     = 0;
	virtual Pdf_setNeighbor& GetNeighbor()    = 0;
	//virtual Check_node&  GetCheckNode()       = 0;
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN) = 0;
	
private:
	//private utility functions 
	virtual PdfBlock* pdfDomain()					= 0;
	virtual PdfMaker* implDomain()					= 0;
	virtual PdfBlock* maker (PdfMaker* pimpl, CLbmCase* pCase)	= 0;
	virtual PdfBlock* maker (PdfMaker* pimpl, CLbmCase* pCase, cgsize_t nvertex) = 0;
	
};
/*--------------------------------------------------------------------------------------------*/
//D2q9_incompDomain Abstract class
//Properties:
//1. Derived from LbmDomain and base to D2q901_incompDomain...
/*--------------------------------------------------------------------------------------------*/
class D2q9_incompDomain:   public LbmDomain
{
public:
		virtual PdfBlock* maker (PdfMaker* pimpl, CLbmCase* pCase);
		virtual PdfBlock* maker (PdfMaker* pimpl, CLbmCase* pCase, cgsize_t nvertex);
};
//////////////////////////////////////////////////////////////////////////////
//BASE FUNCTORS
class Lbm_Initializer
{
public:
	virtual void operator() (PdfBlock* pblock, Domain& domainVariables, cgsize_t mat_index);
};
class Lbm_writer
{
public:
	virtual void operator() (PdfBlock* pblock, std::fstream& outputFile);
};

class Lbm_reader
{
public:
	virtual void operator() (PdfBlock* pblock, std::fstream& inputFile);
};
class Srt_collider
{
public:
	virtual void operator() (PdfBlock* pblock, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
};

class Srt_collider_improved
{
public:
	virtual void operator() (PdfBlock* pblock, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index);
};

class Trt_collider
{
public:
	virtual void operator() (PdfBlock* pblock, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
};

class Mrt_collider
{
public:
	virtual void operator() (PdfBlock* pblock, PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
};

class Pdf_streamer
{
public:
	virtual void operator() (PdfBlock* pblock, cgsize_t node);
};

class Pdf_setNeighbor
{
public:
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index) = 0;
	//void set_Neighbor(PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index, cgsize_t f_opp);
};

class D2Q900_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

class D2Q901_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q902_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q903_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q904_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q905_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q906_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q907_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D2Q908_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
/*==============================================================================*/
class Pdf_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type) = 0;
	Node::NodeValueType_t halfWay_BB(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node);
	Node::NodeValueType_t moving_BB(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node);
	Node::NodeValueType_t extrapolation_BB(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node);
	Node::NodeValueType_t open_NBC(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node);
	Node::NodeValueType_t open_CBC(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, Node::NodeValueType_t vel);
	virtual Node::NodeValueType_t ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
			Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound) = 0;
};

/*==============================================================================*/
class D2Q9_bounder: public Pdf_bounder
{
public:
	virtual Node::NodeValueType_t ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
						   Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound);
	Node::NodeValueType_t ZhouHe_diag(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t fI, cgsize_t fII, cgsize_t node,
		Node::NodeValueType_t c1, Node::NodeValueType_t c2, Node::NodeValueType_t c3,  cgsize_t xvel, cgsize_t yvel, LbmDomain::BOUND b_type, Node::LBMBOUND bound);
};
//DERIVED FUNCTORS FOR D2Q9
/*=============================================================================*/
class D2Q900_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q900_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q900_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q900_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
/*-----------------------------------------------------------------------------*/

class D2Q901_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);

};
class D2Q901_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q901_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q901_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D2Q902_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D2Q902_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q902_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q902_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D2Q903_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D2Q903_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q903_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};

class D2Q903_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D2Q904_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D2Q904_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q904_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q904_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D2Q905_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D2Q905_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q905_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q905_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D2Q906_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D2Q906_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q906_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q906_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D2Q907_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D2Q907_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q907_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q907_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D2Q908_velBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q908_BBack: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D2Q908_periodic: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type){};
};
class D2Q908_openBound: public D2Q9_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};


/*-----------------------------------------------------------------------------*/
////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------------------------------*/
//D2Q900_incompDomain concrete class
//Properties:
//1. Derived from D2Q9Incomp 
//2. Implementation of process specific to D2q900Incomp objects
/*--------------------------------------------------------------------------------------------*/
class D2Q900_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q900_incompDomain(CLbmCase* pCase);
	D2Q900_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode()      {return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock* pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD2Q9_00 impl_00;
	//addresses of pdf_00 domain 
	PdfBlock* pdf_00Domain;
	D2Q900_velBound  v_boundary;
	D2Q900_BBack wall_fbb;
	D2Q900_periodic  io_periodic;
	D2Q900_openBound open_bdry_f;
	D2Q900_setNeighbor  neighbor;
	
};
/*--------------------------------------------------------------------------------------------*/
//D2Q901_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q901_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q901_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q901_incompDomain(CLbmCase* pCase);
	D2Q901_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD2Q9_01 impl_01;
	// addresses of pdf_01 domain 
	PdfBlock*        pdf_01Domain;
	D2Q901_velBound  v_boundary;
	D2Q901_BBack wall_fbb;
	D2Q901_periodic  io_periodic;
	D2Q901_openBound open_bdry_f;
	D2Q901_setNeighbor  neighbor;
	//D2Q901_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D2Q902_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q902_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q902_incompDomain(CLbmCase* pCase);
	D2Q902_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_02 impl_02;
	//addresses of pdf_02 domain 
	PdfBlock* pdf_02Domain;
	D2Q902_velBound  v_boundary;
	D2Q902_BBack wall_fbb;
	D2Q902_periodic  io_periodic;
	D2Q902_openBound open_bdry_f;
	D2Q902_setNeighbor  neighbor;
	//D2Q902_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D2Q903_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q903_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q903_incompDomain(CLbmCase* pCase);
	D2Q903_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_03 impl_03;
	//addresses of pdf_03 domain 
	PdfBlock* pdf_03Domain;
	D2Q903_velBound  v_boundary;
	D2Q903_BBack wall_fbb;
	D2Q903_periodic  io_periodic;
	D2Q903_openBound open_bdry_f;
	D2Q903_setNeighbor  neighbor;
	//D2Q903_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D2Q904_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q904_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q904_incompDomain(CLbmCase* pCase);
	D2Q904_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_04 impl_04;
	//addresses of pdf_04 domain 
	PdfBlock* pdf_04Domain;
	D2Q904_velBound  v_boundary;
	D2Q904_BBack wall_fbb;
	D2Q904_periodic  io_periodic;
	D2Q904_openBound open_bdry_f;
	D2Q904_setNeighbor  neighbor;
	//D2Q904_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D2Q905_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q905_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q905_incompDomain(CLbmCase* pCase);
	D2Q905_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_05 impl_05;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_05Domain;
	D2Q905_velBound  v_boundary;
	D2Q905_BBack wall_fbb;
	D2Q905_periodic  io_periodic;
	D2Q905_openBound open_bdry_f;
	D2Q905_setNeighbor  neighbor;
//	D2Q905_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D2Q906_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q906_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q906_incompDomain(CLbmCase* pCase);
	D2Q906_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_06 impl_06;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_06Domain;
	D2Q906_velBound  v_boundary;
	D2Q906_BBack wall_fbb;
	D2Q906_periodic  io_periodic;
	D2Q906_openBound open_bdry_f;
	D2Q906_setNeighbor  neighbor;
	//D2Q906_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D2Q907_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q907_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q907_incompDomain(CLbmCase* pCase);
	D2Q907_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_07 impl_07;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_07Domain;
	D2Q907_velBound  v_boundary;
	D2Q907_BBack wall_fbb;
	D2Q907_periodic  io_periodic;
	D2Q907_openBound open_bdry_f;
	D2Q907_setNeighbor  neighbor;
	//D2Q907_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D2Q908_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D2Q908_incompDomain:   public D2q9_incompDomain
{
public:
	D2Q908_incompDomain(CLbmCase* pCase);
	D2Q908_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor&  GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder(){return io_periodic;};
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode(){return check_node;};
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD2Q9_08 impl_08;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_08Domain;
	D2Q908_velBound  v_boundary;
	D2Q908_BBack wall_fbb;
	D2Q908_periodic  io_periodic;
	D2Q908_openBound open_bdry_f;
	D2Q908_setNeighbor  neighbor;
	//D2Q908_checkNode check_node;
};

#endif

