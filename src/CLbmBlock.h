#pragma once

#ifndef CLBMBLOCK_H
#define CLBMBLOCK_H

#include <cstdlib>
#include "CLbmCase.h"
#include <vector>
#include <map>
#include <omp.h>

/*--------------------------------------------------------------------------------------------*/
//Node concrete class
//Properties:
//1. Stand alone for now
/*--------------------------------------------------------------------------------------------*/
class PdfDomain;
class CVariable;
class CLbmCase;

class Node
{
public:
	typedef double   NodeValueType_t;
	//ENUMERATIONS
	/////////////////////////////////////////////////////////
	//Enumeration for the different boundaries
	enum   LBMBOUND     {
		EMPTY, WEST, EAST, SOUTH, NORTH, BOTTOM, TOP,
		WB, WT, WS, WN, EB, ET, ES, EN,
		SB, NB, ST, NT, 
		WSB, ESB, WST, EST, WNB, ENB, WNT, ENT
	};
	//Enumeration for boundary node type
	enum   LBMBOUNDPDF {NONE, BURIED, INWARD};
	//class Bad_Node {}; //exception class
	////////////////////////////////////////////////////////
	//ctor
	Node( cgsize_t m_nodeNum =0, NodeValueType_t m_nodeVal = 0.0L);
	//Convert ctor
	Node(NodeValueType_t val);
	//Copy ctor
	Node(const Node& initB);
	//overloaded =  Operator
	Node& operator= (const Node& rhs);
	//overloaded += operator
	Node& operator+= (const Node& rhs);
	//overloaded *= operator
	Node& operator*= (const Node& rhs);
	//overloaded *= operator
	Node& operator/= (const Node& rhs);
	//overloaded == operator
	bool operator== (const Node& rhs);
	//functions that write node fields
	void node_write(std::fstream& outputFile);
	void node_read (std::fstream& inputFile);
	//function that set default state
	static void set_default(cgsize_t nodeNum, NodeValueType_t m_nodeVal);
	cgsize_t          m_nodeNum;
	cgsize_t          m_neighborNum;
	cgsize_t          m_periodicNeighborNum;
	cgsize_t          m_connectTo;
	NodeValueType_t   m_nodeVal;
	NodeValueType_t   m_oldNodeVal;
	NodeValueType_t   m_postColsn;
	static Node       default_node;
	NodeValueType_t   m_xcoord;
	NodeValueType_t   m_ycoord;
	NodeValueType_t   m_zcoord;
	
	bool              m_linkUpdated;
	bool              m_valueUpdated;
	bool              m_streamStart;
	bool              m_nodeStreamed;
	bool              m_periodicUpdated;
	bool              m_wettingFlag;
	bool              m_isSolid;
	bool              m_isExtrapSolid;
	bool              m_isExtrapNotSolid;
	bool              m_isFullBB;
	bool              m_boundFlagUpdated;
	LBMBOUND          m_bound;
	LBMBOUNDPDF       m_boundType;
};

//HELPER FUNCTION FOR NODE
//=============================================================================================//
Node operator+ (Node& lhs, Node& rhs);
Node operator* (Node& lhs, Node& rhs);
Node operator/ (Node& lhs, Node& rhs);




/*--------------------------------------------------------------------------------------------*/
//PdfBlock Abstract class
//Properties:
//1. Base class of the probability density function block heirachy
/*--------------------------------------------------------------------------------------------*/
class PdfBlock
{
public:
	typedef std::vector<Node*>    NodeType_t;
	typedef std::vector<cgsize_t> BoundStreamNodes_t;
	
	//virtual ~PdfBlock(){};

	//overloaded assignment operator
	virtual PdfBlock& operator= (PdfBlock& rhs);
	//overloaded + operation
	virtual PdfBlock&  operator+ (PdfBlock& rhs);
	//overloaded += operation
	virtual PdfBlock&  operator+= (PdfBlock& rhs);
	//overloaded * operation
	virtual PdfBlock&  operator* (PdfBlock& rhs);
	//overloaded * operation
	virtual PdfBlock&  operator* (const Node::NodeValueType_t scalar);

	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs, PdfBlock& rhs);

	//function to reset node value flag
	void ResetNodeValueFlag();
	void ResetValue();

	//public function to read fields
	NodeType_t*           GetPdfunction();
	NodeType_t*           GetPdfunction_eq();
	NodeType_t*           GetPdfmoment();
	NodeType_t*           GetPdfmoment_eq();
	CLbmCase*             GetCaseInfo();
	Node::NodeValueType_t GetCx();
	Node::NodeValueType_t GetCy();
	Node::NodeValueType_t GetCz();
	Node::NodeValueType_t GetWk();
	BoundStreamNodes_t*   GetBoundNodes();
	BoundStreamNodes_t*   GetStreamNodes();

	void setBlock             (CLbmCase* pCase);
	void setBlock             (cgsize_t nvertex);
	void setStreamBound       (CLbmCase* pCase, PdfBlock* oppBlock);
	virtual void setNodeBType (cgsize_t  node)  = 0;
	//virtual void setPdfSwitch (cgsize_t  node)  = 0;
	//member function: sets initial value of Pdfs
	void initialize_Pdf(Domain& domainVariables, cgsize_t mat_index);
	//member function: implements collision
	virtual void srt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f) = 0;
	virtual void srt_collide_improved(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index) = 0;
	virtual void trt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f) = 0;
	virtual void mrt_collide(PdfDomain* allPdfs, Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f) = 0;
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node) = 0;

	//member function: implements streaming
	void stream_domain(cgsize_t node);
	//member functions that writes and read nodes
	void lbm_write(std::fstream& outputFile);
	void lbm_read (std::fstream& inputFile);

private:
	//functions to read fields
	virtual NodeType_t*           pdfunction()     = 0;
	virtual NodeType_t*           pdfunction_eq()  = 0;
	virtual NodeType_t*           pdfmoment()      = 0;
	virtual NodeType_t*           pdfmoment_eq()   = 0;
	virtual CLbmCase*             caseInfo()       = 0;
	virtual Node::NodeValueType_t c_x()            = 0;
	virtual Node::NodeValueType_t c_y()            = 0;
	virtual Node::NodeValueType_t c_z()            = 0;
	
	void StreamScanner  (cgsize_t* pPnts, cgsize_t numPnts, PdfBlock* oppBlock); 
	virtual void setLinkage        (CLbmCase* pCase) = 0;
	virtual void setNodeCoord      (CLbmCase* pCase) = 0;
	virtual Node::NodeValueType_t  GetWeight()       = 0;
};
//HELPER FOR PDFBLOCK
//=============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//Lbm_D2Q9 Abstract class
//Properties:
//1. Derived from PdfBlock and base toIncomp_D2Q9...
/*--------------------------------------------------------------------------------------------*/
class Lbm_D2Q9: public PdfBlock
{
public:
	//member function
	virtual void setLinkage        (CLbmCase* pCase);
	virtual void setNodeCoord      (CLbmCase* pCase);
	virtual Node::NodeValueType_t  GetWeight();

	//move to individual pdfs later
	//virtual void setPdfSwitch (cgsize_t  node){};
private:
	virtual Node::NodeValueType_t weightFactor() = 0;
	virtual void setQuad_4conn(cgsize_t i)       = 0;
};

/*--------------------------------------------------------------------------------------------*/
//Incomp_D2Q9 Abstract class
//Properties:
//1. Derived from  Lbm_D2Q9 and base to D2q900_incomp...
/*--------------------------------------------------------------------------------------------*/
class Incomp_D2Q9: public  Lbm_D2Q9
{
public:
	//member function
	virtual void srt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void srt_collide_improved(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index);
	virtual void trt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void mrt_collide(PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node) {return 0.0;};//remove later

private:
	//flip direction
	void flipper(Domain& domainVariables, cgsize_t f, cgsize_t f_opp, cgsize_t node);
	
};

/*--------------------------------------------------------------------------------------------*/
//Comp_D2Q9 Abstract class
//Properties:
//1. Derived from  Lbm_D2Q9 and base to D2q901_incomp, D2q901_comp ...
/*--------------------------------------------------------------------------------------------*/
class Comp_D2Q9: public  Lbm_D2Q9
{
public:
	//member function
	virtual void srt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void srt_collide_improved(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index);
	virtual void trt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual void mrt_collide(PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f);
	virtual Node::NodeValueType_t compute_moment(PdfDomain* allPdfs, cgsize_t node) {return 0.0;};//remove later
};

/*--------------------------------------------------------------------------------------------*/
//D2q900_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q900_incomp:   public Incomp_D2Q9
{	
public:
	//non-default ctor
	D2q900_incomp (CLbmCase* pCase);
	D2q900_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor:using default ctor
	~D2q900_incomp();

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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q901_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q901_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q901_incomp (CLbmCase* pCase);
	D2q901_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q901_incomp();
	
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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q902_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q902_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q902_incomp (CLbmCase* pCase);
	D2q902_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q902_incomp();
	
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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q903_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q903_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q903_incomp (CLbmCase* pCase);
	D2q903_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q903_incomp();
	
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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q904_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q904_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q904_incomp (CLbmCase* pCase);
	D2q904_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q904_incomp();
	
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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q905_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q905_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q905_incomp (CLbmCase* pCase);
	D2q905_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q905_incomp();

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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q906_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q906_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q906_incomp (CLbmCase* pCase);
	D2q906_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q906_incomp();
	
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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);

};

/*--------------------------------------------------------------------------------------------*/
//D2q907_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q907_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q907_incomp (CLbmCase* pCase);
	D2q907_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q907_incomp();

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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q908_incomp concrete class
//Properties:
//1. Derived from Incomp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q908_incomp:   public Incomp_D2Q9
{
public:
	//non-default ctor
	D2q908_incomp (CLbmCase* pCase);
	D2q908_incomp (CLbmCase* pCase, cgsize_t nvertex);
	//copy ctor: using default
	~D2q908_incomp();

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
	virtual Node::NodeValueType_t c_z(){return 0.0;};
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
	Node::NodeValueType_t weight;

	//private utility functions
	virtual void setQuad_4conn(cgsize_t i);
};

/*--------------------------------------------------------------------------------------------*/
//D2q900_comp concrete class
//Properties:
//1. Derived from D2q9_00. Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q900_comp:   public Comp_D2Q9
{
public:
	//non-default ctor
	D2q900_comp (CLbmCase* pCase);
	//copy ctor: using default
	~D2q900_comp();
		//overloaded assignment operator
	D2q900_comp& operator= (D2q900_comp& rhs);

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

private:
	NodeType_t pdf;
	NodeType_t pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*  m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t weight;
};

/*--------------------------------------------------------------------------------------------*/
//D2q901_comp concrete class
//Properties:
//1. Derived from Comp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q901_comp:   public Comp_D2Q9
{
public:
	//non-default ctor
	D2q901_comp (CLbmCase* pCase);
		//copy ctor: using default
	~D2q901_comp();
	//overloaded assignment operator
	D2q901_comp& operator= (D2q901_comp& rhs);

	//overloaded == operator
	friend bool      operator== (PdfBlock& lhs,  PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction();
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();

private:
	NodeType_t pdf;
	NodeType_t pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*   m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t weight;
};

/*--------------------------------------------------------------------------------------------*/
//D2q902_comp concrete class
//Properties:
//1. Derived from Comp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q902_comp:   public Comp_D2Q9
{
public:
	//non-default ctor
	D2q902_comp (CLbmCase* pCase);
	//copy ctor: using default
	~D2q902_comp();
		//overloaded assignment operator
	D2q902_comp& operator= (D2q902_comp& rhs);

	//overloaded == operator
	friend bool      operator== (const PdfBlock& lhs, const PdfBlock& rhs);
	
	//functions to read fields
	virtual NodeType_t* pdfunction() ;
	virtual NodeType_t* pdfunction_eq();
	virtual NodeType_t* pdfmoment();
	virtual NodeType_t* pdfmoment_eq();
	virtual CLbmCase*   caseInfo();
	virtual Node::NodeValueType_t weightFactor();
	virtual Node::NodeValueType_t c_x();
	virtual Node::NodeValueType_t c_y();

private:
	NodeType_t pdf;
	NodeType_t pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*   m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t weight;
};

/*--------------------------------------------------------------------------------------------*/
//D2q903_comp concrete class
//Properties:
//1. Derived from Comp_D2Q9 Implementation level class
/*--------------------------------------------------------------------------------------------*/
class D2q903_comp:   public Comp_D2Q9
{
public:
	//non-default ctor
	D2q903_comp (CLbmCase* pCase);
	//copy ctor: using default
	~D2q903_comp();
		//overloaded assignment operator
	D2q903_comp& operator= (D2q903_comp& rhs);

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

private:
	NodeType_t pdf;
	NodeType_t pdf_eq;

	NodeType_t   moment;
	NodeType_t   moment_eq;

	CLbmCase*  m_pCase;

	//lattice parameter
	Node::NodeValueType_t cx;
	Node::NodeValueType_t cy;
	Node::NodeValueType_t weight;
};

//=============================================================================================//
//Base class for creating  lbmBlock Objects
class PdfMaker
{
public:
	//API for creating blocks
	virtual PdfBlock* incomp(CLbmCase* pCase)  = 0;
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex) = 0;
	//virtual PdfBlock*   comp(CLbmCase* pCase)  = 0;
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9 derived abstract class
//Properties:
//1. Derived from PdfMaker
//2. Base to MakeD2Q9_00...
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9: public PdfMaker
{
public:
	//API for creating blocks
	virtual PdfBlock* incomp(CLbmCase* pCase)   = 0;
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex)  = 0;
	//virtual PdfBlock*   comp(CLbmCase* pCase)   = 0;
		
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_00 derived concrete class
//Properties:
//1. implementation/leaf level class;
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_00: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
	
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_01 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_01: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_02 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_02: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_03 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_03: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_04 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_04: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_05 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_05: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_06 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_06: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_07 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_07: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_08 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
class MakeD2Q9_08: public MakeD2Q9
{
public:
	virtual PdfBlock* incomp(CLbmCase* pCase);
	virtual PdfBlock* incomp(CLbmCase* pCase, cgsize_t nvertex);
	//virtual PdfBlock*   comp(CLbmCase* pCase);
};

/*--------------------------------------------------------------------------------------------*/
//HELPERS FOR MAPPING OF NODE LOCATION
/*--------------------------------------------------------------------------------------------*/
bool mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, PdfDomain* allPdfs, cgsize_t node);

bool mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									 cgsize_t fV, cgsize_t fVI, PdfDomain* allPdfs, cgsize_t node);

bool mapClassAEdge(cgsize_t fI, cgsize_t fII, cgsize_t fdiagI, cgsize_t fdiagII, PdfDomain* allPdfs, cgsize_t node);

bool mapClassAEdge(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									 cgsize_t fV, cgsize_t fVI, cgsize_t f_orthoI, cgsize_t f_orthoII, PdfDomain* allPdfs, cgsize_t node);

bool mapClassACorner(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									 cgsize_t fV, cgsize_t fVI, cgsize_t fVII, PdfDomain* allPdfs, cgsize_t node);
					

#endif

