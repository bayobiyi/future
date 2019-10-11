#pragma once

#ifndef CLBMBLOCK3D_H
#define CLBMBLOCK3D_H

#include <cstdlib>
#include "CLbmBlock.h"

//=============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//Lbm_D3Q15 Abstract class
//Properties:
//1. Derived from PdfBlock and base to Incomp_D3Q15...
/*--------------------------------------------------------------------------------------------*/
//utility struct
/*----------------------------------------------*/
struct Velocity3D
{
	Node::NodeValueType_t x;
	Node::NodeValueType_t y;
	Node::NodeValueType_t z;
};

struct Force3D
{
	Node::NodeValueType_t x;
	Node::NodeValueType_t y;
	Node::NodeValueType_t z;
};

struct D3Q15_SourceMoments
{
	Node::NodeValueType_t source00;
	Node::NodeValueType_t source01;
	Node::NodeValueType_t source02;
	Node::NodeValueType_t source03;
	Node::NodeValueType_t source04;
	Node::NodeValueType_t source05;
	Node::NodeValueType_t source06;
	Node::NodeValueType_t source07;
	Node::NodeValueType_t source08;
	Node::NodeValueType_t source09;
	Node::NodeValueType_t source10;
	Node::NodeValueType_t source11;
	Node::NodeValueType_t source12;
	Node::NodeValueType_t source13;
	Node::NodeValueType_t source14;
};

struct D3Q15_EquilMoments
{
	Node::NodeValueType_t momentEq00;
	Node::NodeValueType_t momentEq01;
	Node::NodeValueType_t momentEq02;
	Node::NodeValueType_t momentEq03;
	Node::NodeValueType_t momentEq04;
	Node::NodeValueType_t momentEq05;
	Node::NodeValueType_t momentEq06;
	Node::NodeValueType_t momentEq07;
	Node::NodeValueType_t momentEq08;
	Node::NodeValueType_t momentEq09;
	Node::NodeValueType_t momentEq10;
	Node::NodeValueType_t momentEq11;
	Node::NodeValueType_t momentEq12;
	Node::NodeValueType_t momentEq13;
	Node::NodeValueType_t momentEq14;
};

struct D3Q15_Minv
{
	Node::NodeValueType_t Minv0;
	Node::NodeValueType_t Minv1;
	Node::NodeValueType_t Minv2;
	Node::NodeValueType_t Minv3;
	Node::NodeValueType_t Minv4;
	Node::NodeValueType_t Minv5;
	Node::NodeValueType_t Minv6;
	Node::NodeValueType_t Minv7;
	Node::NodeValueType_t Minv8;
	Node::NodeValueType_t Minv9;
	Node::NodeValueType_t Minv10;
	Node::NodeValueType_t Minv11;
	Node::NodeValueType_t Minv12;
	Node::NodeValueType_t Minv13;
	Node::NodeValueType_t Minv14;
};

struct D3Q15_M
{
	cgsize_t M0;
	cgsize_t M1;
	cgsize_t M2;
	cgsize_t M3;
	cgsize_t M4;
	cgsize_t M5;
	cgsize_t M6;
	cgsize_t M7;
	cgsize_t M8;
	cgsize_t M9;
	cgsize_t M10;
	cgsize_t M11;
	cgsize_t M12;
	cgsize_t M13;
	cgsize_t M14;
};

/*----------------------------------------------*/

class Lbm_D3QN: public PdfBlock
{
public:
	//member function
	void setLinkage           (CLbmCase* pCase);
	void setNodeCoord         (CLbmCase* pCase);
	Node::NodeValueType_t     GetWeight();
	
	Force3D     GetTotalForce(Domain& domainVariables, cgsize_t mat_index, cgsize_t node, Node::NodeValueType_t local_rho);
	Velocity3D  GetComponentVel(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	Velocity3D  GetOverallVel(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	Velocity3D  ComputeVelocity(Force3D totalForce, Velocity3D velocity, Node::NodeValueType_t localRho, cgsize_t M);
	Velocity3D  ComputeVelocity(Force3D totalForce, Velocity3D velocity, Node::NodeValueType_t localRho, Node::NodeValueType_t tau, cgsize_t M);
	Velocity3D  ComputeImproveVel(Force3D totalForce, Velocity3D velocity, Node::NodeValueType_t kinVisco, Node::NodeValueType_t psi, Node::NodeValueType_t sigma);
private:
	virtual Node::NodeValueType_t weightFactor() = 0;
	virtual void setHexa_8conn(cgsize_t i)       = 0;
};

/*--------------------------------------------------------------------------------------------*/
//Incomp_D3Q15 Abstract class
//Properties:
//1. Derived from  Lbm_D3Q15 and base to D3q1500_incomp...
/*--------------------------------------------------------------------------------------------*/


class Incomp_D3Q15: public  Lbm_D3QN
{
public:
	//member function
	virtual void srt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void srt_collide_improved(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index);
	virtual void trt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);

	virtual void mrt_collide(PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	
private:
	//utility functions
	Node::NodeValueType_t Forcing_Part1(Force3D totalForce, Node::NodeValueType_t cs_square);
	Node::NodeValueType_t Forcing_Part2(Force3D totalForce, Velocity3D vel, Node::NodeValueType_t cs_square);
	Node::NodeValueType_t ForcingTerm(Force3D totalForce, Velocity3D vel, Node::NodeValueType_t B_e, Node::NodeValueType_t cs_square);

	Node::NodeValueType_t EquilPDF_Part1(Velocity3D equilVel, Node::NodeValueType_t cs_square);
	Node::NodeValueType_t EquilPDF_Part2(Velocity3D equilVel, Node::NodeValueType_t cs_square);
	Node::NodeValueType_t EquilPDFTerm(Velocity3D equilVel, Node::NodeValueType_t cs_square, Node::NodeValueType_t local_rho);

	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node) = 0;
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node) = 0;
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node) = 0;
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node) = 0;	
};

/*--------------------------------------------------------------------------------------------*/
//Incomp_D3Q19 Abstract class
//Properties:
//1. Derived from  Lbm_D3Q19 and base to D3q1900_incomp...
/*--------------------------------------------------------------------------------------------*/
class Incomp_D3Q19: public  Lbm_D3QN
{
public:
	//member function
	virtual void srt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void srt_collide_improved(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index);
	virtual void trt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void mrt_collide(PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node) {return 0.0;};//remove later
	
private:
	//utility functions
};


/*--------------------------------------------------------------------------------------------*/
//D3q1500_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1500_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1500_incomp (CLbmCase* pCase);
	D3q1500_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor:using default ctor
	~D3q1500_incomp();

	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);

	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node){/*place holder*/};

	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);
	
private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	
	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);

};

/*--------------------------------------------------------------------------------------------*/
//D3q1501_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1501_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1501_incomp (CLbmCase* pCase);
	D3q1501_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1501_incomp();
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);
	
private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	CLbmCase*    m_pCase;
	
	NodeType_t   moment;
	NodeType_t   moment_eq;

	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1502_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1502_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1502_incomp (CLbmCase* pCase);
	D3q1502_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1502_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1503_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1503_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1503_incomp (CLbmCase* pCase);
	D3q1503_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1503_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);
	
private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1504_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1504_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1504_incomp (CLbmCase* pCase);
	D3q1504_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1504_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1505_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1505_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1505_incomp (CLbmCase* pCase);
	D3q1505_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1505_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1506_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1506_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	 D3q1506_incomp (CLbmCase* pCase);
	 D3q1506_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~ D3q1506_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1507_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1507_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1507_incomp (CLbmCase* pCase);
	D3q1507_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1507_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);

	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1508_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1508_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1508_incomp (CLbmCase* pCase);
	D3q1508_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1508_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1508_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1509_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q1509_incomp (CLbmCase* pCase);
	D3q1509_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1509_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q15010_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q15010_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q15010_incomp (CLbmCase* pCase);
	D3q15010_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q15010_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;	

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q15011_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q15011_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q15011_incomp (CLbmCase* pCase);
	D3q15011_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q15011_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q15012_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q15012_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q15012_incomp (CLbmCase* pCase);
	D3q15012_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q15012_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q15013_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q15013_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q15013_incomp (CLbmCase* pCase);
	D3q15013_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q15013_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q15014_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q15 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q15014_incomp:   public Incomp_D3Q15
{
public:
	//non-default ctor
	D3q15014_incomp (CLbmCase* pCase);
	D3q15014_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q15014_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);
	//functions that help with MRT collision
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node);
	virtual Node::NodeValueType_t compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node);
	virtual Node::NodeValueType_t computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node);
	virtual Node::NodeValueType_t computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*    m_pCase;
	BoundStreamNodes_t m_boundNodes;
	BoundStreamNodes_t m_streamNodes;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};


//=============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15 derived abstract class
//Properties:
//1. Derived from PdfMaker
//2. Base to MakeD3Q15_00...
/*--------------------------------------------------------------------------------------------*/
class MakeD3QN: public PdfMaker
{
public:
	//API for creating blocks
	virtual PdfBlock* incomp(CLbmCase* pCase)   = 0;
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex) = 0;
	//virtual PdfBlock*   comp(CLbmCase* pCase)   = 0;
		
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_00 derived concrete class
//Properties:
//1. implementation/leaf level class;
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_00: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
	
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_01 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_01: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_02 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_02: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_03 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_03: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_04 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_04: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_05 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_05: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_06 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_06: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_07 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_07: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_08 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_08: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_09 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_09: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_10 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_10: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_11 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_11: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_12 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_12: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_13 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_13: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_14 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q15_14: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};


//////////////////////////////////////////////////////////////////////////////////////////////////
/*--------------------------------------------------------------------------------------------*/
//D3q1900_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1900_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1900_incomp (CLbmCase* pCase);
	D3q1900_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor:using default ctor
	~D3q1900_incomp();

	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);

	//functions to read fields
	virtual NodeType_t* pdfunction();
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node){/*place holder*/};
	
private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);

};

/*--------------------------------------------------------------------------------------------*/
//D3q1901_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1901_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1901_incomp (CLbmCase* pCase);
	D3q1901_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1901_incomp();
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);


private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;
	
	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1902_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1902_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1902_incomp (CLbmCase* pCase);
	D3q1902_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1902_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1903_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1903_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1903_incomp (CLbmCase* pCase);
	D3q1903_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1903_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;
	
	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1904_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1904_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1904_incomp (CLbmCase* pCase);
	D3q1904_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1904_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1905_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1905_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1905_incomp (CLbmCase* pCase);
	D3q1905_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1905_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1906_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1906_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	 D3q1906_incomp (CLbmCase* pCase);
	 D3q1906_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~ D3q1906_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1907_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1907_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1907_incomp (CLbmCase* pCase);
	D3q1907_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1907_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);

	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1908_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1908_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1908_incomp (CLbmCase* pCase);
	D3q1908_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1908_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;
	
	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1908_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1909_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1909_incomp (CLbmCase* pCase);
	D3q1909_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1909_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1910_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1910_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1910_incomp (CLbmCase* pCase);
	D3q1910_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1910_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1911_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1911_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1911_incomp (CLbmCase* pCase);
	D3q1911_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1911_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1912_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1912_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1912_incomp (CLbmCase* pCase);
	D3q1912_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1912_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;
	
	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1913_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1913_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1913_incomp (CLbmCase* pCase);
	D3q1913_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1913_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);


private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1914_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1914_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1914_incomp (CLbmCase* pCase);
	D3q1914_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1914_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1915_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1915_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1915_incomp (CLbmCase* pCase);
	D3q1915_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1915_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1916_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1916_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1916_incomp (CLbmCase* pCase);
	D3q1916_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1916_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1917_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1917_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1917_incomp (CLbmCase* pCase);
	D3q1917_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1917_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

/*--------------------------------------------------------------------------------------------*/
//D3q1918_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D3q1918_incomp:   public Incomp_D3Q19
{
public:
	//non-default ctor
	D3q1918_incomp (CLbmCase* pCase);
	D3q1918_incomp (CLbmCase* pCase, cgsize_t vertex);
	//copy ctor: using default
	~D3q1918_incomp();
	
	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment(){return 0;};
	virtual NodeType_t* pdfmoment_eq(){return 0;};
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();
	virtual Node::NodeValueType_t c_z();
	virtual void setNodeBType(cgsize_t  node);

private:
	NodeType_t   pdf;
	NodeType_t   pdf_eq;
	NodeType_t   moment;
	NodeType_t   moment_eq;
	CLbmCase*    m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t cz;
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setHexa_8conn(cgsize_t i);
	void bounder1(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
	void bounder2(cgsize_t node, Node::NodeValueType_t u_wall, 
									Node::NodeValueType_t v_wall, Domain& domainVariables);
};

////////////////////////////////////////////////////////////////////////////////////////////////

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_00 derived concrete class
//Properties:
//1. implementation/leaf level class;
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_00: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
	
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_01 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_01: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_02 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_02: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_03 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_03: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_04 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_04: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_05 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_05: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_06 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_06: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_07 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_07: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_08 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_08: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_09 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_09: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_10 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_10: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_11 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_11: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_12 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_12: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_13 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_13: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_14 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_14: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_15 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_15: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_16 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_16: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_17 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_17: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_18 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD3Q19_18: public MakeD3QN
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t  nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};
#endif
