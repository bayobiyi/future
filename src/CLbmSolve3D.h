#pragma once

#ifndef  CLBMSOLVE3D_H
#define  CLBMSOLVE3D_H

#include <cstdlib>
#include "CLbmSolve.h"

//LEVEL ONE:invoked in OnRun()
//==============================================================================================//
//D3QN
class Domain_D3QN: public Domain
{
	public:
	virtual void convergence(Domain& domainVariables, CLbmCase* pCase);
	//virtual void computeCompositeVars(CLbmCase* pCase, Domain& domainVariables);
	virtual void computeCompositeVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	//virtual void computeOverallVars(CLbmCase* pCase, Domain& domainVariables);
	virtual void computeOverallVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	virtual void computeOverallVarsSgl(CLbmCase* pCase, Domain& domainVariables, cgsize_t node);
	//virtual void computeIndividualVars(CLbmCase* pCase, Domain& domainVariables);
	virtual void computeIndividualVars(CLbmCase* pCase, Domain& domainVariables, MaterialType_t::size_type mat_index, cgsize_t node);
	virtual      Convergence*    GetConvergence();
	virtual void SetNumThreads(cgsize_t  num);
	virtual cgsize_t GetNumThreads() const;
	virtual ~Domain_D3QN();
private:
	//utilitty methods
	void setNewToOld(Domain& domainVariables, cgsize_t node);
	void setToZero(Domain& domainVariables, CLbmCase* pCase, cgsize_t node);
	Convergence    m_convergenceData;
	cgsize_t       m_numThreads;
};

//D3Q15
class Domain_D3Q15: public Domain_D3QN
{
public:
	//initializing methods
	virtual void mapStreamingStart(CLbmCase* pCase);
	virtual void mapCardinalDirection(cgsize_t node);
	
	Domain_D3Q15(CLbmCase* pCase);
	Domain_D3Q15(CLbmCase* pCase, cgsize_t  nVertex);
	virtual ~Domain_D3Q15();

	virtual MaterialType_t*  DomainMaterial();
	virtual VariableType_t*  DomainVariables();
	virtual VariableType_t*  DomainTempVariables();
	virtual PatchNodeType_t* DomainSolidObjectNodes();
	virtual MicroPillarsType_t* DomainMicroPillars();

private:
	//initialization utility
	virtual bool mapCardinalSurface(cgsize_t node);
	virtual bool mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t node);
	virtual bool mapCardinalEdge(cgsize_t node);
	virtual bool mapEdge(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t fVI, 
		cgsize_t f_orthoI,cgsize_t f_orthoII, cgsize_t node);
	virtual bool mapCardinalCorner(cgsize_t node);
	virtual bool mapCorner(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									   cgsize_t fV, cgsize_t fVI, cgsize_t fVII, cgsize_t node);

	//field
	MaterialType_t  m_domainMaterial;
	VariableType_t  m_domainVariables;
	VariableType_t  m_domainTempVariables;
	PatchNodeType_t m_domainSolidObjectNodes;
	MicroPillarsType_t m_domainMicroPillars;
};
//D3Q19
class Domain_D3Q19: public Domain_D3QN
{
public:
	//initializing methods
	virtual void mapStreamingStart(CLbmCase* pCase);
	virtual void mapCardinalDirection(cgsize_t node);
	
	Domain_D3Q19(CLbmCase* pCase);
	Domain_D3Q19(CLbmCase* pCase, cgsize_t  nVertex);
	virtual ~Domain_D3Q19();

	virtual MaterialType_t*  DomainMaterial();
	virtual VariableType_t*  DomainVariables();
	virtual VariableType_t*  DomainTempVariables();
	virtual PatchNodeType_t* DomainSolidObjectNodes();
	virtual MicroPillarsType_t* DomainMicroPillars();

private:
	//initialization utility
	bool mapCardinalSurface(cgsize_t node);
	bool mapSurface (cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t node);
	bool mapCardinalEdge(cgsize_t node);
	bool mapEdge (cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t node);
	bool mapCardinalCorner(cgsize_t node);
	bool mapCorner (cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t ffI, cgsize_t ffII, cgsize_t ffIII, cgsize_t node);

	//fields
	MaterialType_t  m_domainMaterial;
	VariableType_t  m_domainVariables;
	VariableType_t  m_domainTempVariables;
	PatchNodeType_t m_domainSolidObjectNodes;
	MicroPillarsType_t m_domainMicroPillars;
};

//LEVEL TWO:invoked from Domain
//==============================================================================================//
//D3QNIncomp_domain
class D3QNIncomp_domain: public PdfDomain
{
public:
	//virtual void updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index);
	virtual void updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t node, PdfDomain::LatticeType_t::size_type f);
	
	void mapBoundAndNeighbor( CLbmCase* pCase, LatticeType_t* lattice);
	
	~D3QNIncomp_domain();
private:
	void setNewToOld(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	void resetToZero(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
};
//D3Q15
class D3Q15Incomp_domain: public D3QNIncomp_domain
{
public:
	D3Q15Incomp_domain(CLbmCase* pCase, cgsize_t  matIndex);
	D3Q15Incomp_domain(CLbmCase* pCase, cgsize_t matIndex, cgsize_t  nvertex);
	virtual LatticeType_t* LatticePdf();
	virtual cgsize_t       MaterialIndex();
	virtual void momentToVelocitySpace(CLbmCase* pCase, Domain& domainVariables);

	~D3Q15Incomp_domain();
private:
	Node::NodeValueType_t convertColumn0(cgsize_t node);
	Node::NodeValueType_t convertColumn1(cgsize_t node);
	Node::NodeValueType_t convertColumn2(cgsize_t node);
	Node::NodeValueType_t convertColumn3(cgsize_t node);
	Node::NodeValueType_t convertColumn4(cgsize_t node);
	Node::NodeValueType_t convertColumn5(cgsize_t node);
	Node::NodeValueType_t convertColumn6(cgsize_t node);
	Node::NodeValueType_t convertColumn7(cgsize_t node);
	Node::NodeValueType_t convertColumn8(cgsize_t node);
	Node::NodeValueType_t convertColumn9(cgsize_t node);
	Node::NodeValueType_t convertColumn10(cgsize_t node);
	Node::NodeValueType_t convertColumn11(cgsize_t node);
	Node::NodeValueType_t convertColumn12(cgsize_t node);
	Node::NodeValueType_t convertColumn13(cgsize_t node);
	Node::NodeValueType_t convertColumn14(cgsize_t node);
	LatticeType_t m_latticeType;
	cgsize_t      m_materialIndex;
	void mapStreamStart(CLbmCase* pCase);
};

//D3Q19
class D3Q19Incomp_domain: public D3QNIncomp_domain
{
public:
	D3Q19Incomp_domain(CLbmCase* pCase, cgsize_t  matIndex);
	D3Q19Incomp_domain(CLbmCase* pCase, cgsize_t matIndex, cgsize_t  nvertex);
	virtual LatticeType_t* LatticePdf();
	virtual cgsize_t       MaterialIndex();
	virtual void momentToVelocitySpace(CLbmCase* pCase, Domain& domainVariables){/*place holder*/};
	
	~D3Q19Incomp_domain();
private:
	LatticeType_t m_latticeType;
	cgsize_t      m_materialIndex;
	void mapStreamStart(CLbmCase* pCase);
};
//LEVEL THREE: invoked from 
//==============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//D3q15_incompDomain Abstract class
//Properties:
//1. Derived from LbmDomain and base to D3q1500_incompDomain...
/*--------------------------------------------------------------------------------------------*/
class D3qN_incompDomain:   public LbmDomain
{
public:
		virtual PdfBlock* maker (PdfMaker* pimpl, CLbmCase* pCase);
		virtual PdfBlock* maker (PdfMaker* pimpl, CLbmCase* pCase, cgsize_t nvertex);
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////
//DERIVED FUNCTORS FOR D3Q15
/*=============================================================================*/
class D3Q1500_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1501_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1502_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1503_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1504_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1505_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1506_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1507_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1508_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1509_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q15010_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q15011_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q15012_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q15013_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q15014_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

class D3Q15_bounder: public Pdf_bounder
{
public:
	virtual Node::NodeValueType_t ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
																		Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound);
	virtual Node::NodeValueType_t ZhouHe_diag(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t fI, cgsize_t fII, cgsize_t ffI, cgsize_t ffII, 
		cgsize_t node, Node::NodeValueType_t c1, Node::NodeValueType_t c2, Node::NodeValueType_t c3, Node::NodeValueType_t c4, Node::NodeValueType_t c5,
		cgsize_t xvel, cgsize_t yvel, cgsize_t zvel, LbmDomain::BOUND b_type, Node::LBMBOUND bound);

};

class D3Q1500_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1500_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

class D3Q1500_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1500_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D3Q1501_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1501_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1501_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1501_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D3Q1502_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1502_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1502_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1502_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D3Q1503_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);

};

class D3Q1503_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1503_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1503_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1504_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1504_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1504_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1504_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1505_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1505_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1505_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1505_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1506_velBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1506_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1506_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1506_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1507_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1507_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1507_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1507_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1508_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1508_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1508_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1508_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1509_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);

};
class D3Q1509_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1509_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1509_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q15010_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15010_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15010_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15010_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q15011_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q15011_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15011_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15011_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q15012_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15012_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15012_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15012_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q15013_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q15013_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15013_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15013_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q15014_velBound: public D3Q15_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q15014_BBack: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15014_Periodic: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q15014_openBound: public D3Q15_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

//////////////////////////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------------------------------*/
//D3Q1500_incompDomain concrete class
//Properties:
//1. Derived from D3Q15Incomp 
//2. Implementation of process specific to D3q1500Incomp objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1500_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1500_incompDomain(CLbmCase* pCase);
	D3Q1500_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock* pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q15_00 impl_00;
	//addresses of pdf_00 domain 
	PdfBlock* pdf_00Domain;
	D3Q1500_velBound  v_boundary;
	D3Q1500_BBack wall_fbb;
	D3Q1500_Periodic  io_periodic;
	D3Q1500_openBound open_bdry_f;
	D3Q1500_setNeighbor  neighbor;
	//D3Q1500_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1501_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1501_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1501_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1501_incompDomain(CLbmCase* pCase);
	D3Q1501_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q15_01 impl_01;
	// addresses of pdf_01 domain 
	PdfBlock* pdf_01Domain;
	D3Q1501_velBound  v_boundary;
	D3Q1501_BBack wall_fbb;
	D3Q1501_Periodic  io_periodic;
	D3Q1501_openBound open_bdry_f;
	D3Q1501_setNeighbor  neighbor;
	//D3Q1501_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1502_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1502_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1502_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1502_incompDomain(CLbmCase* pCase);
	D3Q1502_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_02 impl_02;
	//addresses of pdf_02 domain 
	PdfBlock* pdf_02Domain;
	D3Q1502_velBound  v_boundary;
	D3Q1502_BBack wall_fbb;
	D3Q1502_Periodic  io_periodic;
	D3Q1502_openBound open_bdry_f;
	D3Q1502_setNeighbor  neighbor;
	//D3Q1502_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1503_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1503_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1503_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1503_incompDomain(CLbmCase* pCase);
	D3Q1503_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_03 impl_03;
	//addresses of pdf_03 domain 
	PdfBlock* pdf_03Domain;
	D3Q1503_velBound  v_boundary;
	D3Q1503_BBack wall_fbb;
	D3Q1503_Periodic  io_periodic;
	D3Q1503_openBound open_bdry_f;
	D3Q1503_setNeighbor  neighbor;
	//D3Q1503_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1504_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1504_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1504_incompDomain:   public D3qN_incompDomain
{
public:
	 D3Q1504_incompDomain(CLbmCase* pCase);
	 D3Q1504_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_04 impl_04;
	//addresses of pdf_04 domain 
	PdfBlock* pdf_04Domain;
	D3Q1504_velBound  v_boundary;
	D3Q1504_BBack wall_fbb;
	D3Q1504_Periodic  io_periodic;
	D3Q1504_openBound open_bdry_f;
	D3Q1504_setNeighbor  neighbor;
	//D3Q1504_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1505_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1505_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1505_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1505_incompDomain(CLbmCase* pCase);
	D3Q1505_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_05 impl_05;
	//addresses of pdf_05 domain 
	PdfBlock* pdf_05Domain;
	D3Q1505_velBound  v_boundary;
	D3Q1505_BBack wall_fbb;
	D3Q1505_Periodic  io_periodic;
	D3Q1505_openBound open_bdry_f;
	D3Q1505_setNeighbor  neighbor;
	//D3Q1505_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1506_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1506_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1506_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1506_incompDomain(CLbmCase* pCase);
	D3Q1506_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_06 impl_06;
	//addresses of pdf_06 domain 
	PdfBlock* pdf_06Domain;
	D3Q1506_velBound  v_boundary;
	D3Q1506_BBack wall_fbb;
	D3Q1506_Periodic  io_periodic;
	D3Q1506_openBound open_bdry_f;
	D3Q1506_setNeighbor  neighbor;
	//D3Q1506_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1507_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1507_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1507_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1507_incompDomain(CLbmCase* pCase);
	D3Q1507_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_07 impl_07;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_07Domain;
	D3Q1507_velBound  v_boundary;
	D3Q1507_BBack wall_fbb;
	D3Q1507_Periodic  io_periodic;
	D3Q1507_openBound open_bdry_f;
	D3Q1507_setNeighbor  neighbor;
	//D3Q1507_checkNode check_node;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1508_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1508_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1508_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1508_incompDomain(CLbmCase* pCase);
	D3Q1508_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_08 impl_08;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_08Domain;
	D3Q1508_velBound  v_boundary;
	D3Q1508_BBack wall_fbb;
	D3Q1508_Periodic  io_periodic;
	D3Q1508_openBound open_bdry_f;
	D3Q1508_setNeighbor  neighbor;
	//D3Q1508_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1509_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1509_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1509_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1509_incompDomain(CLbmCase* pCase);
	D3Q1509_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_09 impl_09;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_09Domain;
	D3Q1509_velBound  v_boundary;
	D3Q1509_BBack wall_fbb;
	D3Q1509_Periodic  io_periodic;
	D3Q1509_openBound open_bdry_f;
	D3Q1509_setNeighbor  neighbor;
	//D3Q1509_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q15010_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15010_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q15010_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q15010_incompDomain(CLbmCase* pCase);
	D3Q15010_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_10 impl_010;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_010Domain;
	D3Q15010_velBound  v_boundary;
	D3Q15010_BBack wall_fbb;
	D3Q15010_Periodic  io_periodic;
	D3Q15010_openBound open_bdry_f;
	D3Q15010_setNeighbor  neighbor;
	//D3Q15010_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q15011_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15011_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q15011_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q15011_incompDomain(CLbmCase* pCase);
	D3Q15011_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_11 impl_011;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_011Domain;
	D3Q15011_velBound  v_boundary;
	D3Q15011_BBack wall_fbb;
	D3Q15011_Periodic  io_periodic;
	D3Q15011_openBound open_bdry_f;
	D3Q15011_setNeighbor  neighbor;
	//D3Q15011_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q15012_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15012_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q15012_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q15012_incompDomain(CLbmCase* pCase);
	D3Q15012_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_12 impl_012;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_012Domain;
	D3Q15012_velBound  v_boundary;
	D3Q15012_BBack wall_fbb;
	D3Q15012_Periodic  io_periodic;
	D3Q15012_openBound open_bdry_f;
	D3Q15012_setNeighbor  neighbor;
	//D3Q15012_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q15013_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15013_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q15013_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q15013_incompDomain(CLbmCase* pCase);
	D3Q15013_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_13 impl_013;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_013Domain;
	D3Q15013_velBound  v_boundary;
	D3Q15013_BBack wall_fbb;
	D3Q15013_Periodic  io_periodic;
	D3Q15013_openBound open_bdry_f;
	D3Q15013_setNeighbor  neighbor;
	//D3Q15013_checkNode check_node;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q15014_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15014_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q15014_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q15014_incompDomain(CLbmCase* pCase);
	D3Q15014_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q15_14 impl_014;
	//addresses of pdf_01 domain 
	PdfBlock* pdf_014Domain;
	D3Q15014_velBound  v_boundary;
	D3Q15014_BBack wall_fbb;
	D3Q15014_Periodic  io_periodic;
	D3Q15014_openBound open_bdry_f;
	D3Q15014_setNeighbor  neighbor;
	//D3Q15014_checkNode check_node;
};
///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
//DERIVED FUNCTORS FOR D3Q19
/*=============================================================================*/
class D3Q1900_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1901_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1902_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1903_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1904_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1905_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1906_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1907_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1908_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1909_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1910_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1911_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1912_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1913_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};
class D3Q1914_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

class D3Q1915_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

class D3Q1916_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

class D3Q1917_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

class D3Q1918_setNeighbor: public Pdf_setNeighbor
{
	virtual void operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index);
};

//bounder
class D3Q19_bounder: public Pdf_bounder
{
public:
	virtual Node::NodeValueType_t ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
																		Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound);
	virtual Node::NodeValueType_t ZhouHe_diag(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t fI, cgsize_t fII, cgsize_t ffI, cgsize_t ffII, 
		cgsize_t node, Node::NodeValueType_t c1, Node::NodeValueType_t c2, Node::NodeValueType_t c3, Node::NodeValueType_t c4, Node::NodeValueType_t c5,
		cgsize_t xvel, cgsize_t yvel, cgsize_t zvel, LbmDomain::BOUND b_type, Node::LBMBOUND bound);

};

class D3Q1900_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1900_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

class D3Q1900_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1900_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D3Q1901_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1901_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1901_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1901_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D3Q1902_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1902_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1902_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1902_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/

class D3Q1903_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);

};

class D3Q1903_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1903_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1903_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1904_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1904_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1904_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1904_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1905_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1905_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1905_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1905_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1906_Dirichlet: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1906_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1906_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1906_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1907_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1907_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1907_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1907_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1908_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1908_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1908_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1908_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1909_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);

};
class D3Q1909_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1909_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1909_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1910_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1910_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1910_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1910_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1911_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1911_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1911_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1911_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1912_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1912_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1912_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1912_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1913_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1913_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1913_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1913_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1914_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1914_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1914_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1914_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1915_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1915_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1915_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1915_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1916_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1916_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1916_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1916_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1917_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1917_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1917_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1917_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};

/*-----------------------------------------------------------------------------*/
class D3Q1918_Dirichlet: public D3Q19_bounder
{
public:
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
	
};
class D3Q1918_BBack: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1918_Periodic: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};
class D3Q1918_openBound: public D3Q19_bounder
{
	virtual void operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type);
};


//////////////////////////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------------------------------*/
//D3Q1900_incompDomain concrete class
//Properties:
//1. Derived from D3QNIncomp 
//2. Implementation of process specific to D3q1900Incomp objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1900_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1900_incompDomain(CLbmCase* pCase);
	D3Q1900_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN){};
private:
	virtual PdfBlock* pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_00 impl_00;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_00Domain;
	D3Q1900_Dirichlet    dirich_boundary;
	D3Q1900_BBack        wall_bb;
	D3Q1900_Periodic     io_periodic;
	D3Q1900_openBound open_bdry_f;
	D3Q1900_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1901_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1901_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1901_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1901_incompDomain(CLbmCase* pCase);
	D3Q1901_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_01 impl_01;
	//addresses of pdf_01 domain 
	PdfBlock*            pdf_01Domain;
	D3Q1901_Dirichlet    dirich_boundary;
	D3Q1901_BBack        wall_bb;
	D3Q1901_Periodic     io_periodic;
	D3Q1901_openBound open_bdry_f;
	D3Q1901_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1902_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1902_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1902_incompDomain(CLbmCase* pCase);
	D3Q1902_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_02 impl_02;
	//addresses of pdf_02 domain 
	PdfBlock*            pdf_02Domain;
	D3Q1902_Dirichlet    dirich_boundary;
	D3Q1902_BBack        wall_bb;
	D3Q1902_Periodic     io_periodic;
	D3Q1902_openBound open_bdry_f;
	D3Q1902_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1903_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1903_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1903_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1903_incompDomain(CLbmCase* pCase);
	D3Q1903_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_03 impl_03;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_03Domain;
	D3Q1903_Dirichlet    dirich_boundary;
	D3Q1903_BBack        wall_bb;
	D3Q1903_Periodic     io_periodic;
	D3Q1903_openBound open_bdry_f;
	D3Q1903_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1904_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1904_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1904_incompDomain:   public D3qN_incompDomain
{
public:
	 D3Q1904_incompDomain(CLbmCase* pCase);
	 D3Q1904_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_04 impl_04;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_04Domain;
	D3Q1904_Dirichlet    dirich_boundary;
	D3Q1904_BBack        wall_bb;
	D3Q1904_Periodic     io_periodic;
	D3Q1904_openBound open_bdry_f;
	D3Q1904_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1905_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1905_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1905_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1905_incompDomain(CLbmCase* pCase);
	D3Q1905_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_05 impl_05;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_05Domain;
	D3Q1905_Dirichlet    dirich_boundary;
	D3Q1905_BBack        wall_bb;
	D3Q1905_Periodic     io_periodic;
	D3Q1905_openBound    open_bdry_f;
	D3Q1905_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1906_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1906_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1906_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1906_incompDomain(CLbmCase* pCase);
	D3Q1906_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_06 impl_06;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_06Domain;
	D3Q1906_Dirichlet    dirich_boundary;
	D3Q1906_BBack        wall_bb;
	D3Q1906_Periodic     io_periodic;
	D3Q1906_openBound    open_bdry_f;
	D3Q1906_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1907_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1907_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1907_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1907_incompDomain(CLbmCase* pCase);
	D3Q1907_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d or 3d pdf block
	MakeD3Q19_07 impl_07;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_07Domain;
	D3Q1907_Dirichlet    dirich_boundary;
	D3Q1907_BBack        wall_bb;
	D3Q1907_Periodic     io_periodic;
	D3Q1907_openBound    open_bdry_f;
	D3Q1907_setNeighbor  neighbor;
};
/*--------------------------------------------------------------------------------------------*/
//D3Q1908_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1908_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1908_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1908_incompDomain(CLbmCase* pCase);
	D3Q1908_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_08 impl_08;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_08Domain;
	D3Q1908_Dirichlet    dirich_boundary;
	D3Q1908_BBack        wall_bb;
	D3Q1908_Periodic     io_periodic;
	D3Q1908_openBound    open_bdry_f;
	D3Q1908_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1909_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1909_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1909_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1909_incompDomain(CLbmCase* pCase);
	D3Q1909_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_09 impl_09;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_09Domain;
	D3Q1909_Dirichlet    dirich_boundary;
	D3Q1909_BBack        wall_bb;
	D3Q1909_Periodic     io_periodic;
	D3Q1909_openBound    open_bdry_f;
	D3Q1909_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1910_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1910_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1910_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1910_incompDomain(CLbmCase* pCase);
	D3Q1910_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_10 impl_10;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_10Domain;
	D3Q1910_Dirichlet    dirich_boundary;
	D3Q1910_BBack        wall_bb;
	D3Q1910_Periodic     io_periodic;
	D3Q1910_openBound    open_bdry_f;
	D3Q1910_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1911_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1911_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1911_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1911_incompDomain(CLbmCase* pCase);
	D3Q1911_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d 2d pdf block
	MakeD3Q19_11 impl_11;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_11Domain;
	D3Q1911_Dirichlet    dirich_boundary;
	D3Q1911_BBack        wall_bb;
	D3Q1911_Periodic     io_periodic;
	D3Q1911_openBound    open_bdry_f;
	D3Q1911_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1912_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1912_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1912_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1912_incompDomain(CLbmCase* pCase);
	D3Q1912_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	//virtual Check_node&  GetCheckNode();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_12 impl_12;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_12Domain;
	D3Q1912_Dirichlet    dirich_boundary;
	D3Q1912_BBack        wall_bb;
	D3Q1912_Periodic     io_periodic;
	D3Q1912_openBound    open_bdry_f;
	D3Q1912_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1913_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1913_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1913_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1913_incompDomain(CLbmCase* pCase);
	D3Q1913_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_13 impl_13;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_13Domain;
	D3Q1913_Dirichlet    dirich_boundary;
	D3Q1913_BBack        wall_bb;
	D3Q1913_Periodic     io_periodic;
	D3Q1913_openBound    open_bdry_f;
	D3Q1913_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1914_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1914_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1914_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1914_incompDomain(CLbmCase* pCase);
	D3Q1914_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_14 impl_14;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_14Domain;
	D3Q1914_Dirichlet    dirich_boundary;
	D3Q1914_BBack        wall_bb;
	D3Q1914_Periodic     io_periodic;
	D3Q1914_openBound    open_bdry_f;
	D3Q1914_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1915_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1915_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1915_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1915_incompDomain(CLbmCase* pCase);
	D3Q1915_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_15 impl_15;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_15Domain;
	D3Q1915_Dirichlet    dirich_boundary;
	D3Q1915_BBack        wall_bb;
	D3Q1915_Periodic     io_periodic;
	D3Q1915_openBound    open_bdry_f;
	D3Q1915_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1916_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1916_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1916_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1916_incompDomain(CLbmCase* pCase);
	D3Q1916_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_16 impl_16;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_16Domain;
	D3Q1916_Dirichlet    dirich_boundary;
	D3Q1916_BBack        wall_bb;
	D3Q1916_Periodic     io_periodic;
	D3Q1916_openBound    open_bdry_f;
	D3Q1916_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1917_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1917_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1917_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1917_incompDomain(CLbmCase* pCase);
	D3Q1917_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_17 impl_17;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_17Domain;
	D3Q1917_Dirichlet    dirich_boundary;
	D3Q1917_BBack        wall_bb;
	D3Q1917_Periodic     io_periodic;
	D3Q1917_openBound    open_bdry_f;
	D3Q1917_setNeighbor  neighbor;
};

/*--------------------------------------------------------------------------------------------*/
//D3Q1918_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1918_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
class D3Q1918_incompDomain:   public D3qN_incompDomain
{
public:
	D3Q1918_incompDomain(CLbmCase* pCase);
	D3Q1918_incompDomain(CLbmCase* pCase, cgsize_t nvertex);
	virtual Pdf_bounder& GetVelBounder();
	virtual Pdf_bounder& GetFullBBBounder();
	virtual Pdf_setNeighbor& GetNeighbor();
	virtual Pdf_bounder& GetPeriodicBounder();
	virtual Pdf_bounder& GetOpenBounder();
	virtual void BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN);
private:
	virtual PdfBlock*  pdfDomain();
	virtual PdfMaker* implDomain();
	//Enviroments that makes incomp, comp, ... 2d pdf block
	MakeD3Q19_18 impl_18;
	//addresses of pdf_00 domain 
	PdfBlock*            pdf_18Domain;
	D3Q1918_Dirichlet    dirich_boundary;
	D3Q1918_BBack        wall_bb;
	D3Q1918_Periodic     io_periodic;
	D3Q1918_openBound    open_bdry_f;
	D3Q1918_setNeighbor  neighbor;
};

#endif
