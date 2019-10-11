#include "CLbmBlock.h"
#include "CLbmSolve.h"
#include <vector>
#include <cmath>
#include <iostream>

/*-----------------------------------------------------------------------------------------*/
//Node concrete class
//Properties:
//1. Stand alone for now
/*-----------------------------------------------------------------------------------------*/
//ctor
Node::Node( cgsize_t nodeNum, NodeValueType_t val)
{
	m_nodeNum      = nodeNum   ? nodeNum   : default_node.m_nodeNum;
	m_nodeVal      = val       ? val       : default_node.m_nodeVal;
	m_linkUpdated  = false;
	m_valueUpdated = false;
	m_periodicUpdated = false;
	m_wettingFlag = false;
	m_streamStart  = false;
	m_isSolid      = false;
	m_isExtrapSolid= false;
	m_isFullBB     = false;
	m_boundFlagUpdated= false;
	m_nodeStreamed        = false;
	m_oldNodeVal          = 0.0;
	m_postColsn           = 0.0;
	m_xcoord              = 0.0;
	m_ycoord              = 0.0;
	m_zcoord              = 0.0;
	m_connectTo           = 0;
	m_neighborNum         = 0;
	m_periodicNeighborNum = 0;
	m_bound        = Node::EMPTY;
	m_boundType    = Node::NONE;
};

//copy constructor
Node::Node(const Node& initB)
{
	m_nodeNum     = initB.m_nodeNum;
	m_nodeVal     = initB.m_nodeVal;
	m_oldNodeVal  = initB.m_oldNodeVal;
}

//convert ctor
Node::Node(NodeValueType_t val)
{
	m_nodeVal =  val;

}
//overloaded = Operator
Node& Node::operator= (const Node& rhs) 
{
	if(this != &rhs) //beware of self-assignament: rhs = rhs;
	{
		m_oldNodeVal = m_nodeVal;
		m_nodeVal    = rhs.m_nodeVal;
	}
	return *this;
}
//overloaded += Operator
//inline 
Node& Node::operator+= (const Node& rhs) 
{
	m_nodeVal += rhs.m_nodeVal;
	return *this;
}

//overloaded *= Operator
//inline 
Node& Node::operator*= (const Node& rhs) 
{
	m_nodeVal *= rhs.m_nodeVal;
	return *this;
}

//overloaded /= Operator
Node& Node::operator/= (const Node& rhs) 
{
	m_nodeVal /= rhs.m_nodeVal;
	return *this;
}

//overloaded == Operator
//inline 
bool Node::operator== (const Node& rhs)
{
	return m_nodeNum == rhs.m_nodeNum && m_nodeVal == rhs.m_nodeVal;
}

//Default node state
Node Node::default_node(0,0.0);

//function that set the default state
void 
Node::set_default(cgsize_t nodeNum, NodeValueType_t val)
{
	Node::default_node = Node(nodeNum,val);
}
void 
Node::node_write(std::fstream& outputFile)
{
	outputFile<< m_nodeNum <<' '<<	m_nodeVal <<' '<< m_oldNodeVal <<' '<< m_linkUpdated <<' '<<	m_valueUpdated 
						<<' '<< m_streamStart <<' '<<	m_xcoord <<' '<<	m_ycoord   <<' '<<	m_zcoord <<' '<< m_connectTo 
						<<' '<< m_neighborNum <<' '<< m_periodicNeighborNum <<' '<<	m_oldNodeVal <<' '<<	m_nodeStreamed   
						<<' '<<	m_isSolid <<' '<< m_isExtrapSolid<<' '<< m_isExtrapNotSolid <<' '<<	m_isFullBB <<' '<<	m_boundFlagUpdated
						<<' '<< static_cast<int>(m_bound) <<' '<<	static_cast<int>(m_boundType)<<' ';
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"  node           = "<< m_nodeNum <<"\n";
	std::cout<<"  m_nodeVal      = "<<m_nodeVal<<"\n";
	std::cout<<"  m_oldNodeVal   = "<< m_oldNodeVal <<"\n";
	std::cout<<"  m_isSolid node = "<< m_isSolid <<"\n";
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void 
Node::node_read (std::fstream& inputFile)
{
	int bound       = 0;
	int boundtype   = 0;
	inputFile >>m_nodeNum >> m_nodeVal >> m_oldNodeVal >> m_linkUpdated >>	m_valueUpdated 
						>> m_streamStart >>	m_xcoord >> m_ycoord   >>	m_zcoord >>	m_connectTo 
						>> m_neighborNum >> m_periodicNeighborNum  >>	m_oldNodeVal >> m_nodeStreamed   
						>>	m_isSolid >>	m_isExtrapSolid >> m_isExtrapNotSolid >>	m_isFullBB >> m_boundFlagUpdated 
						>> bound >> boundtype	;
	m_bound     = static_cast<LBMBOUND>(bound);
	m_boundType = static_cast<LBMBOUNDPDF>(boundtype);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"  node           = "<< m_nodeNum <<"\n";
	std::cout<<"  m_nodeVal      = "<<m_nodeVal<<"\n";
	std::cout<<"  m_oldNodeVal   = "<< m_oldNodeVal <<"\n";
	std::cout<<"  m_isSolid node = "<< m_isSolid <<"\n";
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//if(m_nodeNum == 8000)
			//std::cout<<" pdf for node = "<<m_nodeNum<<" after uploading is "<< m_nodeVal<<"\n";
}
//HELPER FUNCTION FOR NODE
//==============================================================================================//
Node operator+ (Node& lhs, Node& rhs)
{
	Node r = lhs;
	return r+= rhs;
}

Node operator* (Node& lhs, Node& rhs)
{
	Node r = lhs;
	return r*= rhs;
}

Node operator/ (Node& lhs, Node& rhs)
{
	Node r = lhs;
	return r/= rhs;
}
/*--------------------------------------------------------------------------------------------*/
//PdfBlock Abstract class
//Properties:
//1. Base class of the block heirachy
/*--------------------------------------------------------------------------------------------*/

//overloaded = operator
//inline 
PdfBlock&  
PdfBlock::operator= (PdfBlock& rhs)
{
	if(this != &rhs )  //beware of self-assignament: rhs = rhs
	{                                                 
		for (NodeType_t::size_type i = 0; i < GetPdfunction()->size(); i++)
		{
			if(GetPdfunction()->at(i)->m_nodeNum != rhs.GetPdfunction()->at(i)->m_nodeNum)
			{
				//output some kind of error or exception
				break;
			}
			*(GetPdfunction()->at(i)) = *(rhs.GetPdfunction()->at(i));
		}
	}
	return *this;
}
//overloaded + operation
PdfBlock&  
PdfBlock::operator+ (PdfBlock& rhs)
{	
	for (PdfBlock::NodeType_t::size_type i =0; i < GetPdfunction()->size(); i++)
	{
		//+= to avoid temporary objects
		*(GetPdfunction()->at(i)) += (*(rhs.GetPdfunction()->at(i)));
	}
	return *this;
}

//overloaded += operation
PdfBlock&  
PdfBlock::operator+= (PdfBlock& rhs)
{	
	for (PdfBlock::NodeType_t::size_type i =0; i < GetPdfunction()->size(); i++)
	{
		//+= to avoid temporary objects
		*(GetPdfunction()->at(i)) += (*(rhs.GetPdfunction()->at(i)));
	}
	return *this;
}

//overloaded * operation
PdfBlock&  
PdfBlock::operator* (PdfBlock& rhs)
{	
	for (PdfBlock::NodeType_t::size_type i = 0; i < GetPdfunction()->size(); i++)
	{
		//*= to avoid temporary objects
		*(GetPdfunction()->at(i)) *= (*(rhs.GetPdfunction()->at(i)));
	}
	return *this;
}

PdfBlock&  
PdfBlock::operator* (const Node::NodeValueType_t scalar)
{	
	for (PdfBlock::NodeType_t::size_type i = 0; i < GetPdfunction()->size(); i++)
	{
		//*= to avoid temporary objects
		*(GetPdfunction()->at(i)) *= scalar;
	}
	return *this;
}

//inline 
PdfBlock::NodeType_t* 
PdfBlock::GetPdfunction()
{
	return pdfunction();
}

//inline 
PdfBlock::NodeType_t* 
PdfBlock::GetPdfunction_eq()
{
	return pdfunction_eq();
}
PdfBlock::NodeType_t* 
PdfBlock::GetPdfmoment()
{
	return pdfmoment();
}

//inline 
PdfBlock::NodeType_t* 
PdfBlock::GetPdfmoment_eq()
{
	return pdfmoment_eq();
}
CLbmCase*           
PdfBlock::GetCaseInfo()
{
	return caseInfo();
}

Node::NodeValueType_t 
PdfBlock::GetCx()
{
	return c_x();
}

Node::NodeValueType_t 
PdfBlock::GetCy()
{
	return c_y();
}

Node::NodeValueType_t 
PdfBlock::GetCz()
{
	return c_z();
}

Node::NodeValueType_t 
PdfBlock::GetWk()
{
	return GetWeight();
}

void 
PdfBlock::ResetNodeValueFlag()
{
	for(NodeType_t::size_type i =1; i <= GetPdfunction()->size(); i++)
	{
		GetPdfunction()->at(i-1)->m_valueUpdated    = false;
	}
}

void 
PdfBlock::ResetValue()
{
	for(NodeType_t::size_type i =1; i <= GetPdfunction()->size(); i++)
	{
		GetPdfunction()->at(i-1)->m_nodeVal    = 0.0;
	}
}

//member function: implements streaming
void
PdfBlock::stream_domain(cgsize_t node)
{
	cgsize_t  fromNode = GetPdfunction()->at(node-1)->m_connectTo;
	if(fromNode != 0)
	{
		GetPdfunction()->at(node-1)->m_nodeVal    = GetPdfunction()->at(fromNode-1)->m_nodeVal;
		GetPdfunction()->at(node-1)->m_oldNodeVal = GetPdfunction()->at(fromNode-1)->m_nodeVal;
		GetPdfunction()->at(node-1)->m_nodeStreamed = true;
		GetPdfunction()->at(node-1)->m_valueUpdated = true;
		stream_domain(fromNode);
	}
	else if(fromNode == 0)
		GetPdfunction()->at(node-1)->m_valueUpdated = false;
}
//member functions that read and write nodes
void 
PdfBlock::initialize_Pdf(Domain& domainVariables, cgsize_t mat_index)
{
	Node::NodeValueType_t nodeValue(0.0), density(0.0);
	cgsize_t rho(4);
	for(NodeType_t::size_type i =1; i <= GetPdfunction()->size(); i++)
	{
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" ****** "<<" node =  "<< i <<" ****** "<<"\n";
			std::cout<<" x-coord  "<< GetPdfunction()->at(i-1)->m_xcoord<<"\n";
			std::cout<<" y-coord  "<< GetPdfunction()->at(i-1)->m_ycoord<<"\n";
			std::cout<<" z-coord  "<< GetPdfunction()->at(i-1)->m_zcoord<<"\n";
			std::cout<<" ConnectTo  "<< GetPdfunction()->at(i-1)->m_connectTo<<"\n";
			std::cout<<" neighbor  "<< GetPdfunction()->at(i-1)->m_neighborNum<<"\n";
			std::cout<<" neighbor  "<< GetPdfunction()->at(i-1)->m_periodicNeighborNum<<"\n";
			std::cout<<" bound  "<< GetPdfunction()->at(i-1)->m_bound<<"\n";
			std::cout<<" isSolid  "<< GetPdfunction()->at(i-1)->m_isSolid<<"\n";
			std::cout<<" isInward/Buried  "<< GetPdfunction()->at(i-1)->m_boundType<<"\n";
			std::cout<<" streamStart  "<< GetPdfunction()->at(i-1)->m_streamStart<<"\n";
			std::cin.get();
		 /////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		density   = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(i-1)->m_nodeVal;

		nodeValue = density * GetWk();
		//set density value to section
		GetPdfunction()->at(i-1)->m_nodeVal         = nodeValue;
		GetPdfunction()->at(i-1)->m_oldNodeVal      = nodeValue;
		GetPdfunction()->at(i-1)->m_valueUpdated    = false;	
		/*if(i == 620000)
		{
			std::cout<<" rho  = "<<density<<"\n";
			std::cout<<" fold = "<<GetPdfunction()->at(i-1)->m_nodeVal<<"\n";
			std::cout<<" fnew = "<<GetPdfunction()->at(i-1)->m_oldNodeVal<<"\n";
		}*/
	}
}

//member functions that read and write nodes
void 
PdfBlock::lbm_write(std::fstream& outputFile)
{
	for(NodeType_t::size_type i =1; i <= GetPdfunction()->size(); i++)
	{
		GetPdfunction()->at(i-1)->node_write(outputFile);
		//GetPdfunction_eq()->at(i-1)->node_write(outputFile);
	}
}
void 
PdfBlock::lbm_read (std::fstream& inputFile)
{
	for(NodeType_t::size_type i =1; i <= GetPdfunction()->size(); i++)
	{
		GetPdfunction()->at(i-1)->node_read(inputFile);
		//GetPdfunction_eq()->at(i-1)->node_read(inputFile);
	}
}
//private utility functions
void 
PdfBlock::setBlock(CLbmCase* pCase)
{
	//Add Node objects to the pdf container
	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
	{
		GetPdfunction()->push_back   (new Node(i));
		GetPdfmoment()->push_back   (new Node(i));
		GetPdfunction_eq()->push_back(new Node(i));
		
	}
	//set node coordinates
	setNodeCoord(pCase);
	setLinkage(pCase);
}
void 
PdfBlock::setBlock(cgsize_t nvertex)
{
	//Add Node objects to the pdf container
	for(cgsize_t i = 1; i <= nvertex; i++)
	{
		GetPdfunction()->push_back   (new Node(i));
		GetPdfmoment()->push_back   (new Node(i));
		GetPdfunction_eq()->push_back(new Node(i));
	}
}

void 
PdfBlock::setStreamBound(CLbmCase* pCase, PdfBlock* oppBlock)
{
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0);
	for(cgsize_t i=1; i <= pCase->m_boco.m_nBocos; i++)
	{
		pPnts    = *(pCase->m_boco.m_pnts  + (i-1));
		numPnts  = *(pCase->m_boco.m_npnts + (i-1));
		StreamScanner(pPnts, numPnts, oppBlock);
	}
}

void
PdfBlock::StreamScanner(cgsize_t* pPnts, cgsize_t numPnts, PdfBlock* oppBlock)
{
	cgsize_t node = *pPnts;
	for(cgsize_t i =0; i < numPnts; i++)
	{
		node = *(pPnts + i);
		if(oppBlock->GetPdfunction()->at(node-1)->m_boundType == Node::INWARD)
		{
			if(GetPdfunction()->at(node-1)->m_streamStart == false)
					GetPdfunction()->at(node-1)->m_streamStart = true;
		}
	}
}
 
//============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//Lbm_D2Q9 Abstract class
//Properties:
//1. Derived from PdfBlock and base to Incomp_D2Q9...
/*--------------------------------------------------------------------------------------------*/

//inline 
Node::NodeValueType_t 
Lbm_D2Q9::GetWeight()
{
	return weightFactor();
}

void 
Lbm_D2Q9::setNodeCoord(CLbmCase* pCase)
{
	cgsize_t x(0), y(1);
	// Add the boundary name(s) to boundary nodes
	for(NodeType_t::size_type i = 1; i <= GetPdfunction()->size(); i++)
	{
		GetPdfunction()->at(i-1)->m_xcoord = *(*(pCase->m_grid.m_coord + x)  + (i-1));
		GetPdfunction()->at(i-1)->m_ycoord = *(*(pCase->m_grid.m_coord + y)  + (i-1));
	}
}

void
Lbm_D2Q9::setLinkage(CLbmCase* pCase)
{
	ElementType_t eltType;
	for(CLbmElement::MaterialIndex_t::iterator iter = pCase->m_element.m_materialIndex.begin();
		  iter != pCase->m_element.m_materialIndex.end(); iter++)
	{
		eltType = *(pCase->m_element.m_elementType + *iter);
		switch (eltType)
		{
		case QUAD_4: 
			setQuad_4conn(*iter);
			break;
		default:
			break;
		}
	}
}
//PDFBLOCK HELPER FUNCTION
/*--------------------------------------------------------------------------------------------*/
//Incomp_D2Q9 Abstract class
//Properties:
//1. Derived from Lbm_D2Q9 and base to D2q900_incomp...
/*--------------------------------------------------------------------------------------------*/
//Member function that implements collision
//flip direction
void 
Incomp_D2Q9::flipper(Domain& domainVariables, cgsize_t f, cgsize_t f_opp, cgsize_t node)
{
	//place holder
}
void 
Incomp_D2Q9::srt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	
}
void 
	Incomp_D2Q9::srt_collide_improved(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index)
{

}
void 
Incomp_D2Q9::trt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	
}
void 
Incomp_D2Q9::mrt_collide(PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	
}

/*--------------------------------------------------------------------------------------------*/
//Comp_D2Q9 Abstract class
//Properties:
//1. Derived from Lbm_D2Q9 and base to D2q901_comp...
/*--------------------------------------------------------------------------------------------*/
//member function: implements collision

void 
Comp_D2Q9::srt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{

}
void 
Comp_D2Q9::srt_collide_improved(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index)
{

}
void 
Comp_D2Q9::trt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	
}
void 
Comp_D2Q9::mrt_collide(PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	
}
//FRIENDS TO CONCRETE CLASSES
//==========================================================================================//
//overloaded == operator
bool 
operator== (PdfBlock& lhs, PdfBlock& rhs)
{
		bool temp = true;
		PdfBlock::NodeType_t::const_iterator iter1 = lhs.GetPdfunction()->begin();
		PdfBlock::NodeType_t::const_iterator iter2 = rhs.GetPdfunction()->begin();
		
		while ( iter1 != lhs.GetPdfunction()->end() && temp == true)
		{
			if(*iter1++ == *iter2++)
			{
				temp = true;
			}
			else
			{
				temp = false;
				break;
			}
		}
		return temp;
}
//============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//D2q900_incomp concrete class
//Properties:
//1. Derived from D2q9_00. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q900_incomp::D2q900_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0)
{
	weight = 4.0/9.0;
	setBlock(pCase);
}
D2q900_incomp::D2q900_incomp (CLbmCase* pCase, cgsize_t nvertex):m_pCase(pCase),
cx(0.0), cy(0.0)
{
	weight = 4.0/9.0;
	setBlock(nvertex);
}
void
D2q900_incomp::setQuad_4conn(cgsize_t i){
	//NULL BODY
}
/*
void //remove later does not apply to this model
D2q900_incomp::setHexa_8conn(cgsize_t i){
	//NULL BODY
}*/
//copy ctor: using default copy ctor

//dtor
D2q900_incomp::~D2q900_incomp()
{

}

//functions that reads PdfBlock fields

//inline 
PdfBlock::NodeType_t* 
D2q900_incomp::pdfunction() 
{
	return &pdf;
}
//inline 
PdfBlock::NodeType_t* 
D2q900_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q900_incomp::pdfmoment() 
{
	return &moment;
}
//inline 
PdfBlock::NodeType_t* 
D2q900_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
//inline 
CLbmCase*   
D2q900_incomp::caseInfo()
{
	return m_pCase;
}

//inline 
Node::NodeValueType_t 
D2q900_incomp::weightFactor()
{
	return weight;
}

//inline 
Node::NodeValueType_t 
D2q900_incomp::c_x()
{
	return cx;
}

//inline 
Node::NodeValueType_t 
D2q900_incomp::c_y()
{
	return cy;
}

/*--------------------------------------------------------------------------------------------*/
//D2q901_incomp concrete class
//Properties:
//1. Derived from D2q9_01. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q901_incomp::D2q901_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(0.0)
{
	//set block info
	weight = 1.0/9.0;
	setBlock(pCase);
}
D2q901_incomp::D2q901_incomp (CLbmCase* pCase, cgsize_t nvertex):m_pCase(pCase),
cx(1.0), cy(0.0)
{
	//set block info
	weight = 1.0/9.0;
	setBlock(nvertex);
}
void
D2q901_incomp::setQuad_4conn(cgsize_t i){
	cgsize_t nodeL, nodeR, tonodeL, tonodeR;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		tonodeR = *(*(m_pCase->m_element.m_conn + i) + j);
		nodeR   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		nodeL   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeL = *(*(m_pCase->m_element.m_conn + i) + j+3);
		
			if(GetPdfunction()->at(nodeL-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(nodeL-1)->m_connectTo = tonodeL;
				GetPdfunction()->at(nodeL-1)->m_linkUpdated   = true;
			}
			if(GetPdfunction()->at(nodeR-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(nodeR-1)->m_connectTo = tonodeR;
				GetPdfunction()->at(nodeR-1)->m_linkUpdated   = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q901_incomp::~D2q901_incomp()
{

}

//function that read fields

PdfBlock::NodeType_t* 
D2q901_incomp::pdfunction() 
{
	return &pdf;
}

PdfBlock::NodeType_t* 
D2q901_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q901_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q901_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D2q901_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D2q901_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D2q901_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D2q901_incomp::c_y()
{
	return cy;
}

void 
D2q901_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::WEST)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}

/*--------------------------------------------------------------------------------------------*/
//D2q902_incomp concrete class
//Properties:
//1. Derived from D2q9_02. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q902_incomp::D2q902_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(1.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(pCase);
}
D2q902_incomp::D2q902_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(1.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(nvertex);
}
void
D2q902_incomp::setQuad_4conn(cgsize_t i){
	//L=leftcorner R=rightcorner
	cgsize_t nodeL, nodeR, tonodeL, tonodeR;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeL   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		nodeR   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeR = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeL = *(*(m_pCase->m_element.m_conn + i) + j);
			if(GetPdfunction()->at(nodeR-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(nodeR-1)->m_connectTo = tonodeR;
				GetPdfunction()->at(nodeR-1)->m_linkUpdated   = true;
			}
			if(GetPdfunction()->at(nodeL-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(nodeL-1)->m_connectTo = tonodeL;
				GetPdfunction()->at(nodeL-1)->m_linkUpdated   = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q902_incomp::~D2q902_incomp()
{

}

//function that read fields

PdfBlock::NodeType_t* 
D2q902_incomp::pdfunction()
{
	return &pdf;
}

PdfBlock::NodeType_t* 
D2q902_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q902_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q902_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D2q902_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D2q902_incomp::weightFactor()
{
	return weight;
}
 
Node::NodeValueType_t 
D2q902_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D2q902_incomp::c_y()
{
	return cy;
}

void 
D2q902_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::SOUTH)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
/*--------------------------------------------------------------------------------------------*/
//D2q903_incomp concrete class
//Properties:
//1. Derived from D2q9_03. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q903_incomp::D2q903_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(0.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(pCase);
}
D2q903_incomp::D2q903_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(0.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(nvertex);
}
void
D2q903_incomp::setQuad_4conn(cgsize_t i){
	cgsize_t nodeL, nodeR, tonodeL, tonodeR;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeL   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeL = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeR = *(*(m_pCase->m_element.m_conn + i) + j+2);
		nodeR   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		
			if(GetPdfunction()->at(nodeL-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(nodeL-1)->m_connectTo = tonodeL;
				GetPdfunction()->at(nodeL-1)->m_linkUpdated   = true;
			}
			if(GetPdfunction()->at(nodeR-1)->m_linkUpdated == false)
			{	
				GetPdfunction()->at(nodeR-1)->m_connectTo = tonodeR;
				GetPdfunction()->at(nodeR-1)->m_linkUpdated   = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q903_incomp::~D2q903_incomp()
{

} 

//function that read fields

PdfBlock::NodeType_t* 
D2q903_incomp::pdfunction() 
{
	return &pdf;
}

PdfBlock::NodeType_t* 
D2q903_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q903_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q903_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D2q903_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D2q903_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D2q903_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D2q903_incomp::c_y()
{
	return cy;
}

void 
D2q903_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::EAST)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
/*--------------------------------------------------------------------------------------------*/
//D2q904_incomp concrete class
//Properties:
//1. Derived from D2q9_04. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q904_incomp::D2q904_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(-1.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(pCase);
}
D2q904_incomp::D2q904_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0), cy(-1.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(nvertex);
}
void
D2q904_incomp::setQuad_4conn(cgsize_t i){
	cgsize_t nodeL, nodeR, tonodeL, tonodeR;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeR   = *(*(m_pCase->m_element.m_conn + i) + j);
		nodeL   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeL = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeR   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		
			if(GetPdfunction()->at(nodeL-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(nodeL-1)->m_connectTo = tonodeL;
				GetPdfunction()->at(nodeL-1)->m_linkUpdated   = true;
			}
			if(GetPdfunction()->at(nodeR-1)->m_linkUpdated == false)
			{	
				GetPdfunction()->at(nodeR-1)->m_connectTo = tonodeR;
				GetPdfunction()->at(nodeR-1)->m_linkUpdated   = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q904_incomp::~D2q904_incomp()
{

} 

//function that read fields

inline PdfBlock::NodeType_t* 
D2q904_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q904_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q904_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q904_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q904_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q904_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q904_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q904_incomp::c_y()
{
	return cy;
}

void 
D2q904_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::NORTH)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}

/*--------------------------------------------------------------------------------------------*/
//D2q905_incomp concrete class
//Properties:
//1. Derived from D2q9_05. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q905_incomp::D2q905_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D2q905_incomp::D2q905_incomp (CLbmCase* pCase, cgsize_t nvertex):m_pCase(pCase),
cx(1.0), cy(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void
D2q905_incomp::setQuad_4conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		tonode   = *(*(m_pCase->m_element.m_conn + i) + j);
		node     = *(*(m_pCase->m_element.m_conn + i) + j+2);

			if(GetPdfunction()->at(node-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(node-1)->m_connectTo   = tonode;
				GetPdfunction()->at(node-1)->m_linkUpdated = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q905_incomp::~D2q905_incomp()
{

} 

//function that read fields

inline PdfBlock::NodeType_t* 
D2q905_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q905_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q905_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q905_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q905_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q905_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q905_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q905_incomp::c_y()
{
	return cy;
}

void 
D2q905_incomp::setNodeBType(cgsize_t  node)
{
	if(
					GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
		)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D2q906_incomp concrete class
//Properties:
//1. Derived from D2q9_06. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q906_incomp::D2q906_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D2q906_incomp::D2q906_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void
D2q906_incomp::setQuad_4conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		tonode   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		node = *(*(m_pCase->m_element.m_conn + i) + j+3);

			if(GetPdfunction()->at(node-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(node-1)->m_connectTo   = tonode;
				GetPdfunction()->at(node-1)->m_linkUpdated = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q906_incomp::~D2q906_incomp()
{

} 

//function that read fields

inline PdfBlock::NodeType_t* 
D2q906_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q906_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q906_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q906_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q906_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q906_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q906_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q906_incomp::c_y()
{
	return cy;
}

void 
D2q906_incomp::setNodeBType(cgsize_t  node)
{
	if(
					GetPdfunction()->at(node-1)->m_bound == Node::EAST
			||	GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::ES
		)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D2q907_incomp concrete class
//Properties:
//1. Derived from D2q9_07. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q907_incomp::D2q907_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D2q907_incomp::D2q907_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void
D2q907_incomp::setQuad_4conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+2);

			if(GetPdfunction()->at(node-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(node-1)->m_connectTo   = tonode;
				GetPdfunction()->at(node-1)->m_linkUpdated = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q907_incomp::~D2q907_incomp()
{

} 

//function that read fields

PdfBlock::NodeType_t* 
D2q907_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D2q907_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q907_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q907_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D2q907_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D2q907_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D2q907_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D2q907_incomp::c_y()
{
	return cy;
}

void 
D2q907_incomp::setNodeBType(cgsize_t  node)
{
	if(
					GetPdfunction()->at(node-1)->m_bound == Node::EAST
			||	GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::EN
		)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D2q908_incomp concrete class
//Properties:
//1. Derived from D2q9_08 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q908_incomp::D2q908_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D2q908_incomp::D2q908_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void
D2q908_incomp::setQuad_4conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+3);

			if(GetPdfunction()->at(node-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(node-1)->m_connectTo   = tonode;
				GetPdfunction()->at(node-1)->m_linkUpdated = true;
			}
	}
} 

//copy ctor: using default copy ctor
//dtor
D2q908_incomp::~D2q908_incomp()
{

} 

//function that read fields

inline PdfBlock::NodeType_t* 
D2q908_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q908_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q908_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q908_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q908_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q908_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q908_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q908_incomp::c_y()
{
	return cy;
}

void 
D2q908_incomp::setNodeBType(cgsize_t  node)
{
	if(
					GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::WN
		)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D2q900_comp concrete class
//Properties:
//1. Derived from D2q9_00. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q900_comp::D2q900_comp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0)
{
	weight = 4.0/9.0;
	setBlock(pCase);
}
//copy ctor: using default copy ctor
//convert ctor
//dtor
D2q900_comp::~D2q900_comp()
{

}
//overloaded = operator
inline D2q900_comp& 
D2q900_comp::operator= (D2q900_comp& rhs)
{
	if(this != &rhs) //beware of self-assignament: rhs = rhs
	{                                                   
		for (NodeType_t::size_type i = 0; i < pdf.size(); i++)
		{
			if(pdf.at(i)->m_nodeNum != rhs.GetPdfunction()->at(i)->m_nodeNum)
			{
				//output some kind of error or exception
				break;
			}
			*(pdf.at(i))    = *(rhs.GetPdfunction()->at(i));
		}
	}
	return *this;
}
//function that read fields

inline PdfBlock::NodeType_t* 
D2q900_comp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q900_comp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q900_comp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q900_comp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q900_comp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q900_comp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q900_comp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q900_comp::c_y()
{
	return cy;
}

/*--------------------------------------------------------------------------------------------*/
//D2q901_comp concrete class
//Properties:
//1. Derived from D2q9_01. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//default ctor
//non-default ctor
D2q901_comp::D2q901_comp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(0.0)
{
	weight = 1.0/9.0;
	setBlock(pCase);
}
//copy ctor: using default copy ctor
//convert ctor
//dtor
D2q901_comp::~D2q901_comp()
{

}
inline D2q901_comp& D2q901_comp::operator= (D2q901_comp& rhs)
{
	if(this != &rhs )  //beware of self-assignament: rhs = rhs
	{                                                   
		for (NodeType_t::size_type i = 0; i < pdf.size(); i++)
		{
			if(pdf.at(i)->m_nodeNum != rhs.GetPdfunction()->at(i)->m_nodeNum)
			{
				//output some kind of error or exception
				break;
			}
			*(pdf.at(i))    = *(rhs.GetPdfunction()->at(i));
		}
	}
	return *this;
}
//function that read fields

inline PdfBlock::NodeType_t* 
D2q901_comp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q901_comp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q901_comp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q901_comp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q901_comp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q901_comp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q901_comp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q901_comp::c_y()
{
	return cy;
}

/*--------------------------------------------------------------------------------------------*/
//D2q902_comp concrete class
//Properties:
//1. Derived from D2q9_02. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q902_comp::D2q902_comp (CLbmCase *pCase): m_pCase(pCase)
	, cx(-1.0), cy(0.0)
{
	weight = 1.0/9.0;
	setBlock(pCase);
}
//copy ctor: using default copy ctor
//convert ctor
//dtor
D2q902_comp::~D2q902_comp()
{

}
inline D2q902_comp& D2q902_comp::operator= (D2q902_comp& rhs)
{
	if(this != &rhs)  //beware of self-assignament: rhs = rhs
	{                                                    
		for (NodeType_t::size_type i = 0; i < pdf.size(); i++)
		{
			if(pdf.at(i)->m_nodeNum != rhs.GetPdfunction()->at(i)->m_nodeNum)
			{
				//output some kind of error or exception
				break;
			}
			*(pdf.at(i))    = *(rhs.GetPdfunction()->at(i));
		}
	}
	return *this;
}
//function that read fields

inline PdfBlock::NodeType_t* 
D2q902_comp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q902_comp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q902_comp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q902_comp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q902_comp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q902_comp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q902_comp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q902_comp::c_y()
{
	return cy;
}

/*--------------------------------------------------------------------------------------------*/
//D2q903_comp concrete class
//Properties:
//1. Derived from D2q9_03. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D2q903_comp::D2q903_comp (CLbmCase *pCase): m_pCase(pCase)
	, cx(0.0), cy(1.0)
{
	weight = 1.0/9.0;
	setBlock(pCase);
}
//copy ctor: using default copy ctor
//dtor
D2q903_comp::~D2q903_comp()
{

}
inline D2q903_comp& D2q903_comp::operator= (D2q903_comp& rhs)
{
	if(this != &rhs)  //beware of self-assignament: rhs = rhs
	{                                                    
		for (NodeType_t::size_type i = 0; i < pdf.size(); i++)
		{
			if(pdf.at(i)->m_nodeNum != rhs.GetPdfunction()->at(i)->m_nodeNum)
			{
				//output some kind of error or exception
				break;
			}
			*(pdf.at(i))    = *(rhs.GetPdfunction()->at(i));
		}
	}
	return *this;
}
//function that read fields
inline PdfBlock::NodeType_t* 
D2q903_comp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D2q903_comp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D2q903_comp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D2q903_comp::pdfmoment_eq() 
{
	return &moment_eq;
}
inline CLbmCase*   
D2q903_comp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D2q903_comp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D2q903_comp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D2q903_comp::c_y()
{
	return cy;
}

//============================================================================================//

/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_00 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_00::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q900_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_00::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q900_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_01 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_01::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q901_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_01::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q901_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_02 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_02::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q902_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_02::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q902_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_03 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_03::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q903_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_03::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q903_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_04 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_04::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q904_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_04::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q904_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_05 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_05::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q905_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_05::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q905_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_06 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_06::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q906_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_06::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q906_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_07 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_07::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q907_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_07::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q907_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD2Q9_08 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD2Q9_08::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D2q908_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD2Q9_08::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D2q908_incomp(pCase, nvertex);
	return pblock;
}

/*--------------------------------------------------------------------------------------------*/
//HELPERS FOR MAPPING OF NODE LOCATION
/*--------------------------------------------------------------------------------------------*/
bool mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, PdfDomain* allPdfs, cgsize_t node)
{
	if(  allPdfs->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
		&& allPdfs->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& allPdfs->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
		&& allPdfs->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0)
		return true;
	else
		return false;
}
bool mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
								cgsize_t fV, cgsize_t fVI, PdfDomain* allPdfs, cgsize_t node)
{
	if(  allPdfs->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& allPdfs->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
		&& allPdfs->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
		&& allPdfs->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
		&& allPdfs->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& allPdfs->GetLatticePdf()->at(fVI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0)
		return true;
	else
		return false;
}

bool mapClassAEdge(cgsize_t fI, cgsize_t fII, cgsize_t fdiagI, cgsize_t fdiagII, PdfDomain* allPdfs, cgsize_t node)
{
	if(  allPdfs->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
		&& allPdfs->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fdiagI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fdiagII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0)
		return true;
	else
		return false;
}

bool mapClassAEdge(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									 cgsize_t fV, cgsize_t fVI, cgsize_t f_orthoI,cgsize_t f_orthoII, PdfDomain* allPdfs, cgsize_t node)
{
	if(  allPdfs->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
		&& allPdfs->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0 
		&& allPdfs->GetLatticePdf()->at(fVI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(f_orthoI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& allPdfs->GetLatticePdf()->at(f_orthoII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0)
		return true;
	else
		return false;
}

bool mapClassACorner(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									   cgsize_t fV, cgsize_t fVI, cgsize_t fVII, PdfDomain* allPdfs, cgsize_t node)
{
	if(  allPdfs->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
		&& allPdfs->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0 
		&& allPdfs->GetLatticePdf()->at(fVI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& allPdfs->GetLatticePdf()->at(fVII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0)
		return true;
	else
		return false;
}
