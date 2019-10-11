#pragma once
#include <cstdlib>
#include "CLbmCase.h"
#include "CLbmBlock.h"
#include <vector>
#include <fstream>

class CVariable
{
public:
	typedef std::vector<Node*>         NodeType_t;
	typedef std::vector<NodeType_t*>   MaterialVarType_t;
	typedef long double                NodeValueType_t;

	//overloaded assignment operator
	CVariable& operator=   (CVariable& rhs);
	//overloaded + operation
	CVariable& operator+   (CVariable& rhs);
	//overloaded /= operation
	CVariable& operator/  (CVariable& rhs);
	//overloaded * operation
	CVariable& operator*  (CVariable& rhs);
	void Update_variables(CLbmCase* pCase, cgsize_t index);
	void Lbm_write(std::fstream& outputFile);
	void Lbm_read (std::fstream& inputFile);
	Node::NodeValueType_t residual();
	void        setBlock(CLbmCase* pCase);
	void        setBlock(cgsize_t nvertex);
	void        setnodeValue(CLbmCase* pCase, CVariable* var, cgsize_t index);
	void        FlagSolidNodes(CLbmCase* pCase, MaterialVarType_t* var, cgsize_t index);
	void        setNodeCoord(CLbmCase* pCase);
	void        fixVelBocoValue   (CLbmCase* pCase, CVariable* var, cgsize_t* pPnts, cgsize_t numPnts, cgsize_t index,  cgsize_t i);
	void        fixWallBocoValue  (CLbmCase* pCase, CVariable* var, cgsize_t* pPnts, cgsize_t numPnts , cgsize_t index, cgsize_t i);
	void        fixPressBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t index,   cgsize_t i);
	void        fixOutflowBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t index, cgsize_t i);
	void        fixGeneralBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t index, cgsize_t i);
	void        resetNodeStatus(CLbmCase* pCase, CVariable* var, cgsize_t i);
	//testing
	virtual void        setPatchedNodesWettability(CLbmCase* pCase) = 0;






	MaterialVarType_t* GetVariable();
	CLbmCase*   GetCaseInfo();
	cgsize_t   GetIndex();
	virtual ~CVariable();

private:
	//functions to read fields
	virtual MaterialVarType_t* variable() = 0;
	virtual CLbmCase*   caseInfo() = 0;
	virtual cgsize_t    index()    = 0;

};

class CVelocity : public CVariable
{
public:
	virtual ~CVelocity();
};

class CXVelocity: public CVelocity
{
public:
	CXVelocity(CLbmCase* pCase);
	CXVelocity(CLbmCase* pCase, cgsize_t  nvertex);
	virtual ~CXVelocity();
	

	//functions to read fields
	virtual MaterialVarType_t* variable();
	virtual CLbmCase*   caseInfo();
	virtual cgsize_t    index();


	//testing
	virtual void        setPatchedNodesWettability(CLbmCase* pCase){};

private:
	CLbmCase*  m_pCase;
	MaterialVarType_t m_velocity;
	cgsize_t   m_index;
};

class CYVelocity: public CVelocity
{
public:
	CYVelocity(CLbmCase* pCase);
	CYVelocity(CLbmCase* pCase, cgsize_t nvertex);
	virtual ~CYVelocity();

	//functions to read fields
	virtual MaterialVarType_t* variable();
	virtual CLbmCase*   caseInfo();
	virtual cgsize_t    index();

	//testing
	virtual void        setPatchedNodesWettability(CLbmCase* pCase){};
private:
	CLbmCase*  m_pCase;
	MaterialVarType_t m_velocity;
	cgsize_t   m_index;
};

class CZVelocity: public CVelocity
{
public:
	CZVelocity(CLbmCase* pCase);
	CZVelocity(CLbmCase* pCase, cgsize_t nvertex);
	virtual ~CZVelocity();
	
	//functions to read fields
	virtual MaterialVarType_t* variable();
	virtual CLbmCase*   caseInfo();
	virtual cgsize_t    index();

//testing
	virtual void        setPatchedNodesWettability(CLbmCase* pCase){};

private:
	CLbmCase*  m_pCase;
	MaterialVarType_t m_velocity;
	cgsize_t   m_index;
};

class CVelocityMag: public CVelocity
{
public:
	CVelocityMag(CLbmCase* pCase);
	CVelocityMag(CLbmCase* pCase, cgsize_t nvertex);
	virtual ~CVelocityMag();

	//functions to read fields
	virtual MaterialVarType_t* variable();
	virtual CLbmCase*   caseInfo();
	virtual cgsize_t    index();

//testing
	virtual void        setPatchedNodesWettability(CLbmCase* pCase){};

private:
	CLbmCase*  m_pCase;
	MaterialVarType_t m_velocity;
	cgsize_t   m_index;
};

struct SphereParameters
{
	Node::NodeValueType_t a;
	Node::NodeValueType_t b;
	Node::NodeValueType_t c;
	Node::NodeValueType_t r_2D;
	Node::NodeValueType_t r_3D;
};
struct RectPrismParameters
{
	Node::NodeValueType_t a;
	Node::NodeValueType_t b;
	Node::NodeValueType_t c;
	Node::NodeValueType_t l;
	Node::NodeValueType_t w;
	Node::NodeValueType_t h;
};
struct CylinderParameters
{
	Node::NodeValueType_t a;
	Node::NodeValueType_t b;
	Node::NodeValueType_t c;
	Node::NodeValueType_t r_2D;
	Node::NodeValueType_t h;
};
class CDensity: public CVariable
{
public:
	typedef std::vector<cgsize_t> PatchNodeType_t;

	CDensity(CLbmCase* pCase);
	CDensity(CLbmCase* pCase, cgsize_t  nvertex);
	virtual ~CDensity();

	//functions for initializing densities
	void setMaterialDensity(CLbmCase* pCase);
	void setPatchedNodesDensity(CLbmCase* pCase);
	virtual void setPatchedNodesWettability(CLbmCase* pCase);

	void PatchStripeWettingNode(cgsize_t* pPnts, cgsize_t numPnts, Node::NodeValueType_t patternsStart, Node::NodeValueType_t patternsEnd, cgsize_t numPatterns);
	void PatchWedgeWettingNode(cgsize_t* pPnts, cgsize_t numPnts, Node::NodeValueType_t patternsStart, Node::NodeValueType_t patternsEnd, cgsize_t numPatterns);
	void PatchWettingNode(cgsize_t node, Node::NodeValueType_t xCheck, Node::NodeValueType_t xNode, Node::NodeValueType_t yNode, 
			Node::NodeValueType_t wedgeBase, cgsize_t currentWedge, cgsize_t* pPnts, cgsize_t numPnts);

	void PatchWettingNodeMirror(cgsize_t* pPnts, cgsize_t numPnts, Node::NodeValueType_t currentYNode, Node::NodeValueType_t currentXNode, 
				Node::NodeValueType_t wedgeBase, cgsize_t currentWedge);



	//function to setPatchedNodes
	void setPatchedNodes(CLbmCase* pCase);
	void setPrismNodes(CLbmCase* pCase);
	void setCylinderNodes(CLbmCase* pCase);
	void setSphereNodes(CLbmCase* pCase);
	PatchNodeType_t*  GetPatchedNodes();

	//function to read fields
	virtual MaterialVarType_t* variable();
	virtual CLbmCase*   caseInfo();
	virtual cgsize_t    index();

private:
	Node::NodeValueType_t computeInitialDensitySphere(CLbmCase* pCase, Node::NodeValueType_t rho_liquid, Node::NodeValueType_t rho_gas, cgsize_t node);
	Node::NodeValueType_t computeInitialDensityFlat(CLbmCase* pCase, Node::NodeValueType_t rho_liquid, Node::NodeValueType_t rho_gas, cgsize_t node);
	CylinderParameters obtainCylinderParameters(CLbmCase* pCase);
	SphereParameters obtainSphereParameters(CLbmCase* pCase);
	RectPrismParameters obtainRectPrismParameters(CLbmCase* pCase);
	void patchZdirection(CLbmCase* pCase, Node::NodeValueType_t r_3D, Node::NodeValueType_t a, Node::NodeValueType_t b, Node::NodeValueType_t c);
	Node::NodeValueType_t matchZcenter(CLbmCase* pCase, Node::NodeValueType_t c);
	CLbmCase*  m_pCase;
	MaterialVarType_t m_density;
	cgsize_t   m_index;
	PatchNodeType_t m_patchedNode;
	
};

class CTemp: public CVariable
{
public:
	CTemp(CLbmCase* pCase);
	CTemp(CLbmCase* pCase, cgsize_t index);
	virtual ~CTemp();

	//function to read fields
	virtual MaterialVarType_t* variable();
	virtual CLbmCase*			caseInfo();
	virtual cgsize_t			index();
//testing
	virtual void        setPatchedNodesWettability(CLbmCase* pCase){};

private:
	CLbmCase*  m_pCase;
	MaterialVarType_t m_temp;
	cgsize_t   m_index;
};

