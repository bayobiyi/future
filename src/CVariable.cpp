#include "CVariable.h"
#include <cmath>
#include <iostream>


CVariable& 
CVariable::operator= (CVariable&  rhs) 
{
	if(this != &rhs )  //beware of self-assignament: rhs = rhs
	{ 
		for(MaterialVarType_t::size_type j = 0; j < GetVariable()->size(); j++)
		{
			for (NodeType_t::size_type i = 0; i < GetVariable()->at(j)->size(); i++)
			{
				if(GetVariable()->at(j)->at(i)->m_nodeNum != rhs.GetVariable()->at(j)->at(i)->m_nodeNum)
				{
					//output some kind of error or exception
					break;
				}
				*(GetVariable()->at(j)->at(i)) = *(rhs.GetVariable()->at(j)->at(i));
			}
		}
	}
	return *this;
}

//overloaded + operation
CVariable& 
CVariable::operator+ (CVariable& rhs)
{	
	for(MaterialVarType_t::size_type j = 0; j < GetVariable()->size(); j++)
	{
		for (PdfBlock::NodeType_t::size_type i = 0; i < GetVariable()->at(j)->size(); i++)
		{
			//+= to avoid temporary objects
			*(GetVariable()->at(j)->at(i)) += (*(rhs.GetVariable()->at(j)->at(i)));
		}
	}
	return *this;
}

CVariable& 
CVariable::operator/ (CVariable& rhs)
{	
	for(MaterialVarType_t::size_type j = 0; j < GetVariable()->size(); j++)
	{
		for (PdfBlock::NodeType_t::size_type i =0; i < GetVariable()->at(j)->size(); i++)
		{
			///= to avoid temporary objects
			*(GetVariable()->at(j)->at(i)) /= (*(rhs.GetVariable()->at(j)->at(i)));
		}
	}
	return *this;
}

//overloaded * operation
CVariable& 
CVariable::operator* (CVariable& rhs)
{
	for(MaterialVarType_t::size_type j = 0; j < GetVariable()->size(); j++)
	{
		for (PdfBlock::NodeType_t::size_type i = 0; i < GetVariable()->at(j)->size(); i++)
		{
			//*= to avoid temporary objects
			*(GetVariable()->at(j)->at(i)) *= (*(rhs.GetVariable()->at(j)->at(i)));
		}
	}
	return *this;
}

void 
CVariable::setBlock(CLbmCase* pCase)
{
	//Add Node objects to the velocity container
	cgsize_t multi(2);
	if(pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		for(cgsize_t j = 0; j < multi; j++)
		{
			variable()->push_back(new NodeType_t);
		}

		for(cgsize_t j = 0; j < multi; j++)
		{
			for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
			{
				variable()->at(j)->push_back(new Node(i));
			}
		}
	}
	else
	{
		variable()->push_back(new NodeType_t);
			for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
			{
				variable()->at(0)->push_back(new Node(i));
			}
	}

	//Flag solid nodes
	for(MaterialVarType_t::size_type k = 0; k < variable()->size(); k++)
	{
		FlagSolidNodes(pCase, variable(), k);
	}
	//Set coordinate flags
	setNodeCoord(pCase);
}

void        
CVariable::setNodeCoord(CLbmCase* pCase)
{
	cgsize_t x(0), y(1), z(2);
	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
		{
			GetVariable()->at(0)->at(i-1)->m_xcoord = *(*(pCase->m_grid.m_coord + x)  + (i-1));
			GetVariable()->at(0)->at(i-1)->m_ycoord = *(*(pCase->m_grid.m_coord + y)  + (i-1));
			GetVariable()->at(0)->at(i-1)->m_zcoord = *(*(pCase->m_grid.m_coord + z)  + (i-1));
		}
}

void 
CVariable::setBlock(cgsize_t nvertex)
{
	//Add Node objects to the velocity container
	if(0) //PLEASE FIX THIS LATER!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	{
		for(cgsize_t j = 0; j < 2; j++)
		{
			for(cgsize_t i = 1; i <= nvertex; i++)
			{
				variable()->at(j)->push_back(new Node(i));
			}
		}
	}
	else
	{
			for(cgsize_t i = 1; i <= nvertex; i++)
			{
				variable()->at(0)->push_back(new Node(i));
			}
	}
}

void 
CVariable::FlagSolidNodes(CLbmCase* pCase, MaterialVarType_t* var, cgsize_t index)
{
	cgsize_t  numPnts(0), node(0);
	cgsize_t interior(2);
	cgsize_t* pPnts = NULL;
	for(cgsize_t j=1; j <= pCase->m_boco.m_nBocos; j++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (j-1)) != interior)
			{
				pPnts    = *(pCase->m_boco.m_pnts  + (j-1));
				numPnts  = *(pCase->m_boco.m_npnts + (j-1));
				for(cgsize_t i =0; i < numPnts; i++)
					{
						node = *(pPnts + i);
						//set node as a solid node if it is a wall node
						if(*(pCase->m_boco.m_bocoType + (j-1)) == BCWall
							&& GetVariable()->at(index)->at(node-1)->m_boundFlagUpdated == false)
							{
								GetVariable()->at(index)->at(node-1)->m_isSolid          = true;
								GetVariable()->at(index)->at(node-1)->m_boundFlagUpdated = true;
							}
					}
		}
	}
}

void 
CVariable::Update_variables(CLbmCase* pCase, cgsize_t index)
{ 
	for(MaterialVarType_t::size_type i = 0; i < GetVariable()->size(); i++)
	{
			for (PdfBlock::NodeType_t::size_type j = 0; j < GetVariable()->at(i)->size(); j++)
			{
				GetVariable()->at(i)->at(j)->m_nodeVal    = pCase->m_sol.m_sol.at(index)->at(i)->at(j);
				GetVariable()->at(i)->at(j)->m_oldNodeVal = pCase->m_sol.m_solOld.at(index)->at(i)->at(j);
			}
	}
}


void 
CVariable::Lbm_write(std::fstream& outputFile)
{
	for(MaterialVarType_t::size_type j = 0; j < GetVariable()->size(); j++)
	{

		for (PdfBlock::NodeType_t::size_type i = 0; i < GetVariable()->at(j)->size(); i++)
		{
			GetVariable()->at(j)->at(i)->node_write(outputFile);
		}
	}
}
void 
CVariable::Lbm_read (std::fstream& inputFile)
{
	for(MaterialVarType_t::size_type j = 0; j < GetVariable()->size(); j++)
	{
		for (PdfBlock::NodeType_t::size_type i = 0; i < GetVariable()->at(j)->size(); i++)
		{
			GetVariable()->at(j)->at(i)->node_read(inputFile);
		}
	}
}

void
CVariable::fixVelBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t* pPnts, cgsize_t numPnts, cgsize_t index, cgsize_t i)
{
	cgsize_t    node(0), xvel(0), yvel(1), zvel(2), density(4);
	NodeValueType_t nodeValue(0.0);
	if     (index == xvel)
		nodeValue = *(pCase->m_boco.m_bocoXVel  + (i-1));
	else if(index == yvel)
		nodeValue = *(pCase->m_boco.m_bocoYVel  + (i-1));
	else if(index == zvel)
		nodeValue = *(pCase->m_boco.m_bocoZVel  + (i-1));
	else if(index == density)
		nodeValue = *(pCase->m_boco.m_bocoPress + (i-1));
	else;

	for(cgsize_t j =0; j < numPnts; j++)
		{
			node = *(pPnts + j);
			for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
				{
					if(var->GetVariable()->at(k)->at(node-1)->m_isSolid == false)
						{
							//fix velocities on boundaries
							var->GetVariable()->at(k)->at(node-1)->m_nodeVal         = nodeValue;
							//write to solution container
							//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)          = nodeValue;
							var->GetVariable()->at(k)->at(node-1)->m_valueUpdated    = true;
						}
				}
		}

}

void
CVariable::fixWallBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t* pPnts, cgsize_t numPnts, cgsize_t index, cgsize_t i)
{
	cgsize_t node(0), moving(2);
	cgsize_t xvel(0), yvel(1), zvel(2), density(4);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" Now fixing wall boundary\n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	for(cgsize_t j =0; j < numPnts; j++)
	{
		node = *(pPnts + j);
		//fix velocities on boundaries
		for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
			{
				if(var->GetVariable()->at(k)->at(node-1)->m_valueUpdated == false)
					{
					if(index != density)
						{
							/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
							std::cout<<" non-density node = "<< node <<"\n";
							std::cin.get();
							////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
							if(*(pCase->m_boco.m_bocoWall + (i-1)) == moving)
								{
									if     (index == xvel)
										{
											var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = *(pCase->m_boco.m_bocoXVel  + (i-1));
											//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)     = *(pCase->m_boco.m_bocoXVel  + (i-1));
										}
									else if(index == yvel)
										{
											var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = *(pCase->m_boco.m_bocoYVel  + (i-1));
											//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)    = *(pCase->m_boco.m_bocoYVel  + (i-1));
										}
									else if(index == zvel)
										{
											var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = *(pCase->m_boco.m_bocoZVel  + (i-1));
											//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)     = *(pCase->m_boco.m_bocoZVel  + (i-1));
										}
								}

							else
								{
									var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = 0.0;
									//write to solutaion arrray
									//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)     = 0.0;

									var->GetVariable()->at(k)->at(node-1)->m_valueUpdated    = true;
								}
							
					}
				}
				else;
		}
	}
}

void
CVariable::fixPressBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t index, cgsize_t i)
{
	cgsize_t    node(0);
	NodeValueType_t nodeValue;
	CVariable::NodeType_t::iterator iter;
	nodeValue = *(pCase->m_boco.m_bocoPress + (i-1));
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" Now fixing pressure for boundary = "<<i<<" and density value =  "<<nodeValue<<"\n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	for(cgsize_t j =1; j <= *(pCase->m_boco.m_npnts + (i-1)); j++)
		{
				node = *(*(pCase->m_boco.m_pnts + (i-1)) + (j-1));
				/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::cout<<" node = "<< node <<"\n";
				std::cin.get();
				////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

				for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
					{
						if(var->GetVariable()->at(k)->at(node-1)->m_isSolid == false)
							{
								/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
								std::cout<<" density for node = "<<node<<" is fixed"<<"\n";
								std::cin.get();
								////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
								//fix pressure on boundaries
								var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = nodeValue;
								//write to solutaion arrray
								//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)     = nodeValue;
								var->GetVariable()->at(k)->at(node-1)->m_valueUpdated    = true;
							}
					}
		}
}

void
CVariable::fixGeneralBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t index, cgsize_t i)
{
	cgsize_t    node(0), rho(4);
	if(index == rho && pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" Now fixing pressure for periodic boundary = "<<i<<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		for(cgsize_t j =1; j <= *(pCase->m_boco.m_npnts + (i-1)); j++)
			{
				node = *(*(pCase->m_boco.m_pnts + (i-1)) + (j-1));

				for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
					{
						if(var->GetVariable()->at(k)->at(node-1)->m_isSolid == false)
							{
								/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
								std::cout<<"Material = "<<k<<" periodic node = "<< node <<"\n";
								std::cin.get();
								////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
								//fix density on boundaries
								var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = var->GetVariable()->at(k)->at(node-1)->m_oldNodeVal;
								/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
								std::cout<<" density = "<<var->GetVariable()->at(k)->at(node-1)->m_nodeVal<<"\n";
								std::cin.get();
								////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
								//write to solutaion arrray
								//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1)     = var->GetVariable()->at(k)->at(node-1)->m_oldNodeVal;
								var->GetVariable()->at(k)->at(node-1)->m_valueUpdated    = true;
							}
					}
			}
	}
}

void
CVariable::fixOutflowBocoValue(CLbmCase* pCase, CVariable* var, cgsize_t index, cgsize_t i)
{
	cgsize_t node(0), xvel(0), yvel(1);
	if(index == xvel || index == yvel)
	{
		for(cgsize_t j =1; j <= *(pCase->m_boco.m_npnts + (i-1)); j++)
			{
				node = *(*(pCase->m_boco.m_pnts + (i-1)) + (j-1));
				for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
					{
						var->GetVariable()->at(k)->at(node-1)->m_nodeVal    = 0.0;
						//write to solutaion arrray
						//pCase->m_sol.m_sol.at(index)->at(k)->at(node-1) = 0.0;
						var->GetVariable()->at(k)->at(node-1)->m_valueUpdated    = true;
					}
			}
	}
}

void
CVariable::setnodeValue(CLbmCase* pCase, CVariable* var, cgsize_t index)
{
	//Add boundary condition value
	cgsize_t density(4);
	CVariable::NodeType_t::iterator iter;
	BCType_t boundType;
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), interior(2);
	for(cgsize_t i =1; i <= pCase->m_boco.m_nBocos; i++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" set boundary node values for boundary = "<< i <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (i-1)) != interior)
			{
				boundType = *(pCase->m_boco.m_bocoType + (i-1));
				pPnts     = *(pCase->m_boco.m_pnts     + (i-1));
				numPnts   = *(pCase->m_boco.m_npnts    + (i-1));
				/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::cout<<" boundary type is = "<< boundType <<"\n";
				std::cin.get();
				////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

				switch (boundType)
					{
						case BCDirichlet:
							{
								if(*(pCase->m_boco.m_bocoPressureBC + (i-1)))
								{
									if(index == density)
									{
										fixPressBocoValue(pCase, var, index, i);
									}
								}
								else if(*(pCase->m_boco.m_bocoExtrapolationBC + (i-1)))
								{
									if(index == density)
									{
										fixPressBocoValue(pCase, var, index, i);
									}
								}
								else
								{
									fixVelBocoValue(pCase, var, pPnts, numPnts, index, i);
								}
							}
						break;
						case BCWall:
							{
								fixWallBocoValue(pCase, var, pPnts, numPnts,index, i);
							}
						break;
						case BCGeneral:
							{
								fixGeneralBocoValue(pCase, var, index, i);
							}
						break;
						case BCOutflow:
							{
								//place holder
							}
						break;
					}
			}
	}

	//reset node status to false
	for(cgsize_t i =1; i <= pCase->m_boco.m_nBocos; i++)
		resetNodeStatus(pCase, var, i);
}

CVariable::MaterialVarType_t* 
CVariable::GetVariable()
{
	return variable();
}

CLbmCase* 
CVariable::GetCaseInfo()
{
	return caseInfo();
}

cgsize_t 
CVariable::GetIndex()
{
	return index();
}

/*Node::NodeValueType_t 
CVariable::residual()
{
	cgsize_t rho(4), primary(0);
	Node::NodeValueType_t error(0), norm(0);
	const Node::NodeValueType_t res = 0.00000001;
	if(GetIndex() == rho)
		return 0.0;
	else
	{
		for(NodeType_t::iterator iter = GetVariable()->at(primary)->begin(); iter != GetVariable()->at(primary)->end(); iter++)
		{
			error += pow(((*iter)->m_nodeVal - (*iter)->m_oldNodeVal), 2);
			norm  += pow((*iter)->m_oldNodeVal, 2);
		}
		return sqrt(error/(norm + res));
	}
}*/

Node::NodeValueType_t 
CVariable::residual()
{
	cgsize_t primary(0);
	Node::NodeValueType_t error(0), norm(0);
	const Node::NodeValueType_t res = 0.00000000000001;
	for(CVariable::NodeType_t::size_type node = 1; node <= GetVariable()->at(primary)->size(); node++)
	{
		error += pow((GetVariable()->at(primary)->at(node-1)->m_nodeVal - GetVariable()->at(primary)->at(node-1)->m_oldNodeVal), 2);
		norm  += pow(GetVariable()->at(primary)->at(node-1)->m_nodeVal, 2);
		//store new in old
		GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	}
	return sqrt(error/(norm + res));
}

inline void
CVariable::resetNodeStatus(CLbmCase* pCase, CVariable* var, cgsize_t i)
{
	cgsize_t node(0);
	for(cgsize_t j =1; j <= *(pCase->m_boco.m_npnts + (i-1)); j++)
		{
			node = *(*(pCase->m_boco.m_pnts + (i-1)) + (j-1));
			for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
				{
					if(var->GetVariable()->at(k)->at(node-1)->m_valueUpdated == true)
						{
							var->GetVariable()->at(k)->at(node-1)->m_valueUpdated = false;
						}
				}
		}
}

CVariable::~CVariable()
{
}

// xvelocity
CVelocity::~CVelocity()
{
}

CXVelocity::CXVelocity(CLbmCase* pCase): m_pCase(pCase), m_index(0)
{
	//Add Node objects to the velocity container and set info
	setBlock(pCase);
	//Add boundary condition value
	setnodeValue(pCase, this, GetIndex());
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished construct cxvelocity\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
CXVelocity::CXVelocity(CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase), m_index(0)
{
	//Add Node objects to the velocity container and set info
	setBlock(nvertex);
}
inline CVariable::MaterialVarType_t* 
CXVelocity::variable() 
{
	return &m_velocity;
}

inline CLbmCase* 
CXVelocity::caseInfo() 
{
	return m_pCase;
}

inline cgsize_t 
CXVelocity::index() 
{
	return m_index;
}

CXVelocity::~CXVelocity()
{

}

CYVelocity::CYVelocity(CLbmCase* pCase): m_pCase(pCase), m_index(1)
{
	setBlock(pCase);
	this->setnodeValue(pCase, this, GetIndex());
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished construct cyvelocity\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
CYVelocity::CYVelocity(CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase), m_index(1)
{
	setBlock(nvertex);
}
inline CVariable::MaterialVarType_t* 
CYVelocity::variable() 
{
	return &m_velocity;
}
inline CLbmCase* 
CYVelocity::caseInfo() 
{
	return m_pCase;
}

inline cgsize_t 
CYVelocity::index() 
{
	return m_index;
}

CYVelocity::~CYVelocity()
{

}

CZVelocity::CZVelocity(CLbmCase* pCase): m_pCase(pCase), m_index(2)
{
	setBlock(pCase);
	setnodeValue(pCase, this, GetIndex());
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished construct czvelocity\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
CZVelocity::CZVelocity(CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase), m_index(2)
{
	setBlock(nvertex);
}
inline CVariable::MaterialVarType_t* 
CZVelocity::variable() 
{
	return &m_velocity;
}
inline CLbmCase* 
CZVelocity::caseInfo() 
{
	return m_pCase;
}

inline cgsize_t 
CZVelocity::index() 
{
	return m_index;
}
CZVelocity::~CZVelocity()
{

}

CVelocityMag::CVelocityMag(CLbmCase* pCase): m_pCase(pCase), m_index(3)
{
	setBlock(pCase);
	cgsize_t    node(0), interior(2);
	NodeValueType_t nodeValue(0.0);
	NodeValueType_t temp(0.0);
	CVariable::NodeType_t::iterator iter;
	for(cgsize_t i =1; i <= m_pCase->m_boco.m_nBocos; i++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (i-1)) != interior)
			{
				temp = pow(*(m_pCase->m_boco.m_bocoXVel + (i-1)),2) + 
							 pow(*(m_pCase->m_boco.m_bocoYVel + (i-1)),2) +
							 pow(*(m_pCase->m_boco.m_bocoZVel + (i-1)),2);

				nodeValue = sqrt(temp);
				for(cgsize_t j =1; j <= *(m_pCase->m_boco.m_npnts + (i-1)); j++)
					{
						node = *(*(m_pCase->m_boco.m_pnts + (i-1)) + (j-1));

						for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
							{
								if(GetVariable()->at(k)->at(node-1)->m_valueUpdated == false
										&& nodeValue != 0.0)
								{
									//set velocity magnitude boco on boundaries
								  GetVariable()->at(k)->at(node-1)->m_nodeVal = nodeValue;
									//write to solution container
									pCase->m_sol.m_sol.at(m_index)->at(k)->at(i-1) = nodeValue;
									GetVariable()->at(k)->at(node-1)->m_valueUpdated = true;
								}

							}
					}
			}
	}
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" finished construct cvelocitymeg\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
CVelocityMag::CVelocityMag(CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase), m_index(3)
{
	setBlock(nvertex);
}
inline CVariable::MaterialVarType_t* 
CVelocityMag::variable() 
{
	return &m_velocity;
}
inline CLbmCase* 
CVelocityMag::caseInfo() 
{
	return m_pCase;
}
CVelocityMag::~CVelocityMag()
{

}

inline cgsize_t 
CVelocityMag::index() 
{
	return m_index;
}

void 
CDensity::setMaterialDensity(CLbmCase* pCase)
{
	//set material density
	NodeValueType_t nodeValue(0.0), wall_nodeValue(0.0), rho_liquid(0.0), rho_gas(0.0), temp(0.0);

	for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
	{
		std::map <cgsize_t,CLbmElement::MaterialRecord_t>::const_iterator iter = m_pCase->m_element.m_sectionsMaterial.find(k);
		nodeValue = iter->second.second.density;
			
		for(cgsize_t i = 1; i <= m_pCase->m_grid.m_NVertex; i++)
		{
			//if droplet
			rho_gas    = 0.017;//nodeValue;
			//if bubble
			//rho_liquid = nodeValue;
			
			temp = computeInitialDensitySphere(pCase, rho_liquid, rho_gas, i);
			//temp = computeInitialDensityFlat(pCase, rho_liquid, rho_gas, i);
	
			//set density value to section
			if(m_density.at(k)->at(i-1)->m_isSolid == false)
			{
				m_density.at(k)->at(i-1)->m_nodeVal    = temp;
				//use m_oldNodeVal to store operating density for node (since it is not updated i.e residual calculation is not performed for density)
				m_density.at(k)->at(i-1)->m_oldNodeVal = temp;
				//write to solutaion arrray
				//pCase->m_sol.m_sol.at(m_index)->at(k)->at(i-1) = nodeValue;
			}
			else
			{
				m_density.at(k)->at(i-1)->m_nodeVal = wall_nodeValue;
				//m_density.at(k)->at(i-1)->m_nodeVal = wall_nodeValue;
			}
		}
	}
}
Node::NodeValueType_t
CDensity::computeInitialDensitySphere(CLbmCase* pCase, Node::NodeValueType_t rho_liquid, Node::NodeValueType_t rho_gas, cgsize_t node )
{
	NodeValueType_t xcoord(0.0), ycoord(0.0), zcoord(0.0);
	NodeValueType_t temp1(0.0), temp2(0.0), temp3(0.0), temp4(0.0);
	NodeValueType_t nodeRadius(0.0), interfaceThickness(5.0);
	xcoord = GetVariable()->at(0)->at(node-1)->m_xcoord;
	ycoord = GetVariable()->at(0)->at(node-1)->m_ycoord;
	zcoord = GetVariable()->at(0)->at(node-1)->m_zcoord;
	SphereParameters sphere = obtainSphereParameters(pCase);
	temp1      = (rho_liquid +rho_gas)/2.0;
	temp2      = (rho_liquid - rho_gas)/2.0;
	nodeRadius = sqrt((xcoord - sphere.a)*(xcoord-sphere.a)+(ycoord - sphere.b)*(ycoord-sphere.b)+(zcoord - sphere.c)*(zcoord-sphere.c));
	temp3      = 2.0*(nodeRadius -sphere.r_3D)/interfaceThickness;
	temp4      = (exp(2*temp3)-1.0) / (exp(2*temp3)+1.0);
	return (temp1 - temp2 * temp4);
}

Node::NodeValueType_t
CDensity::computeInitialDensityFlat(CLbmCase* pCase, Node::NodeValueType_t rho_liquid, Node::NodeValueType_t rho_gas, cgsize_t node )
{
	NodeValueType_t xcoord(0.0);
	NodeValueType_t temp1(0.0), temp2(0.0), temp3(0.0), temp4(0.0), temp5(0.0), temp6(0.0);
	NodeValueType_t interfaceThickness(5.0);
	xcoord = GetVariable()->at(0)->at(node-1)->m_xcoord;
	temp1      = rho_gas;
	temp2      = (rho_liquid - rho_gas)/2.0;
	temp3      = 2.0*(xcoord -25.0)/ interfaceThickness;
	temp4      = 2.0*(xcoord -75.0)/ interfaceThickness;
	temp5      = (exp(2*temp3)-1.0) / (exp(2*temp3)+1.0);
	temp6      = (exp(2*temp4)-1.0) / (exp(2*temp4)+1.0);
	return (temp1 + temp2 * (temp5 - temp6));
}
void 
CDensity::setPatchedNodesDensity(CLbmCase* pCase)
{
	NodeValueType_t nodeValue(0.0), rho_liquid(0.0), rho_gas(0.0), temp(0.0);
	cgsize_t node(0);
	for(MaterialVarType_t::size_type k = 0; k < GetVariable()->size(); k++)
	{
		if      (k == 0)
			nodeValue = 0.3060; //Case->m_patchedNodes.primary
		else if (k == 1)
			nodeValue = 0.0; //pCase->m_patchedNodes.secondary
		else;
		
		for(PatchNodeType_t::size_type j=1; j <= GetPatchedNodes()->size(); j++)
		{
			node = GetPatchedNodes()->at(j-1);
			//if droplet
			rho_liquid    = nodeValue;
			//if bubble
			//rho_gas = nodeValue;
			temp = computeInitialDensitySphere(pCase, rho_liquid, rho_gas, node);
			//temp = computeInitialDensityFlat(pCase, rho_liquid, rho_gas, node);

			if(m_density.at(k)->at(node-1)->m_isSolid == false)
			{
				//set density value to section
				m_density.at(k)->at(node-1)->m_nodeVal         = temp;
				//use m_oldNodeVal to store operating density for node (since it is not updated i.e residual calculation is not performed for density)
				m_density.at(k)->at(node-1)->m_oldNodeVal      = temp;
				//write to solutaion arrray
				//pCase->m_sol.m_sol.at(m_index)->at(k)->at(node-1) =  nodeValue;
			}
		}
	}
}

void
CDensity::setPatchedNodesWettability(CLbmCase* pCase)
{
	Node::NodeValueType_t patternsStart(0.0), patternsEnd(100.0);
	cgsize_t numPatterns(4), numPnts(0);
	BCType_t boundType;
	cgsize_t* pPnts = NULL;

	for(cgsize_t i =1; i <= pCase->m_boco.m_nBocos; i++)
	{
		boundType = *(pCase->m_boco.m_bocoType + (i-1));
		pPnts     = *(pCase->m_boco.m_pnts     + (i-1));
		numPnts   = *(pCase->m_boco.m_npnts    + (i-1));

		if(boundType == BCWall)
		{
			//if patching a stripe
			//PatchStripeWettingNode(pPnts,numPnts, patternsStart, patternsEnd, numPatterns);
			PatchWedgeWettingNode(pPnts, numPnts, patternsStart, patternsEnd, numPatterns);

		}
	}
}

void 
CDensity::PatchWedgeWettingNode(cgsize_t* pPnts, cgsize_t numPnts, Node::NodeValueType_t patternsStart, Node::NodeValueType_t patternsEnd, cgsize_t numPatterns)
{
	Node::NodeValueType_t xCheck(0.0), wedgeBase(0.0), xNode(0.0), yNode(0.0), patternsBaseStart(20.0), slope(22); //slope for 5 wedges was 27.5
	cgsize_t node(0);

	wedgeBase = (patternsEnd - patternsStart) / numPatterns;

	for(cgsize_t wedgeNum = 0; wedgeNum < numPatterns; wedgeNum++)
	{
		for(cgsize_t i=0; i < numPnts; i++)
		{		
			node   =  *(pPnts + i);	
	
			if(GetVariable()->at(0)->at(node-1)->m_wettingFlag == false)
			{	
				//get coordinate of node for comparison
				xNode  =  GetVariable()->at(0)->at(node-1)->m_xcoord;
				yNode  =  GetVariable()->at(0)->at(node-1)->m_ycoord;

				if(yNode >= patternsBaseStart)
				{
					//if wedge gradeint pattern
					xCheck =  (yNode - patternsBaseStart)/slope + (wedgeNum * wedgeBase);
					PatchWettingNode(node, xCheck, xNode,yNode, wedgeBase, wedgeNum, pPnts, numPnts);
				}
				else
				{
					for(MaterialVarType_t::size_type mat_index = 0; mat_index < GetVariable()->size(); mat_index++)
					{
						//GetVariable()->at(mat_index)->at(node-1)->m_wettingFlag = true;
						//m_density.at(0)->at(node-1)->m_nodeVal         = 150;
						//m_density.at(1)->at(node-1)->m_nodeVal         = 0;
					}
				}	
				
			}
		}
	}
}


void 
CDensity::PatchStripeWettingNode(cgsize_t* pPnts, cgsize_t numPnts, Node::NodeValueType_t patternsStart, Node::NodeValueType_t patternsEnd, cgsize_t numPatterns)
{
	Node::NodeValueType_t coordMax(0.0), coordMin(0.0), patternBase(0.0), xNode(0.0), yNode(0.0);
	cgsize_t node(0);
	
	if(numPatterns > 1)
		patternBase = (patternsEnd - patternsStart) / (numPatterns * 2.0);
	else
		patternBase = (patternsEnd - patternsStart) / (numPatterns);

	for(cgsize_t stripeNum = 0; stripeNum < numPatterns; stripeNum++)
	{
		coordMin = patternsStart + (stripeNum * 2.0 * patternBase);
		coordMax = coordMin + patternBase;

		for(cgsize_t i=0; i < numPnts; i++)
		{		
			node   =  *(pPnts + i);	
	
			if(GetVariable()->at(0)->at(node-1)->m_wettingFlag == false)
			{	
				//get coordinate of node for comparison
				xNode  =  GetVariable()->at(0)->at(node-1)->m_xcoord;
				yNode  =  GetVariable()->at(0)->at(node-1)->m_ycoord;

				if((coordMin <= yNode && yNode <= coordMax) || (coordMin <= xNode && xNode <= coordMax))
				{
					for(MaterialVarType_t::size_type mat_index = 0; mat_index < GetVariable()->size(); mat_index++)
					{
						GetVariable()->at(mat_index)->at(node-1)->m_wettingFlag = true;
						//m_density.at(0)->at(node-1)->m_nodeVal         = 150;
						//m_density.at(1)->at(node-1)->m_nodeVal         = 0;
					}
				}
				
			}
		}
	}
}

/*void
CDensity::setPatchedNodesWettability(CLbmCase* pCase)
{
	//std::cout<<"starting wettability patch\n";
	//std::cin.get();
	Node::NodeValueType_t baseStart(190.0), wedgeBase(20.0), wedgeHeight(400.0), xNode(0.0), yNode(0.0), xCheck(0.0), factor(20);
	cgsize_t numWedges(5), node(0);
	BCType_t boundType;
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0);

	for(cgsize_t i =1; i <= pCase->m_boco.m_nBocos; i++)
	{
		boundType = *(pCase->m_boco.m_bocoType + (i-1));
		pPnts     = *(pCase->m_boco.m_pnts     + (i-1));
		numPnts   = *(pCase->m_boco.m_npnts    + (i-1));

		if(boundType == BCWall)
		{
			//std::cout<<"patching wall "<< i <<"\n";
			//std::cin.get();
			for(cgsize_t j=0; j < numWedges; j++)
			{
				//std::cout<<"patching wedge "<< j <<"\n";
				for(cgsize_t k=0; k < numPnts; k++)
				{
					
					node   =  *(pPnts + k);
					xNode  =  GetVariable()->at(0)->at(node-1)->m_xcoord;
					yNode  =  GetVariable()->at(0)->at(node-1)->m_ycoord;
					
					if(GetVariable()->at(0)->at(node-1)->m_wettingFlag == false)
					{
						if(yNode >= baseStart)
						{
							//if wedge gradeint pattern
							xCheck =  (yNode - baseStart)/factor + (j * wedgeBase);
							PatchWettingNode(node, xCheck, xNode,yNode, wedgeBase, j, pPnts, numPnts);
							//PatchStripeWettingNode(node, yNode, numWedges, j);
						}
						else
						{
							for(MaterialVarType_t::size_type mat_index = 0; mat_index < GetVariable()->size(); mat_index++)
							{
								GetVariable()->at(mat_index)->at(node-1)->m_wettingFlag = true;
								m_density.at(0)->at(node-1)->m_nodeVal         = 150;
								m_density.at(1)->at(node-1)->m_nodeVal         = 0;
							}
						}	
					}
				}
			
			}
		}
	}
	//std::cout<<"ending wettability patch\n";
	//std::cin.get();
}*/
void 
CDensity::PatchWettingNode(cgsize_t node, Node::NodeValueType_t xCheck, Node::NodeValueType_t xNode, Node::NodeValueType_t yNode, 
			Node::NodeValueType_t wedgeBase, cgsize_t currentWedge, cgsize_t* pPnts, cgsize_t numPnts)
{
	if(xCheck <= xNode && xNode <= (0.5*wedgeBase + currentWedge*wedgeBase))
	{
		for(MaterialVarType_t::size_type mat_index = 0; mat_index < GetVariable()->size(); mat_index++)
		{
			GetVariable()->at(mat_index)->at(node-1)->m_wettingFlag = true;
			//m_density.at(0)->at(node-1)->m_nodeVal         = 0.3060;
			//m_density.at(1)->at(node-1)->m_nodeVal         = 0;
		}
		PatchWettingNodeMirror(pPnts, numPnts, yNode, xNode, wedgeBase, currentWedge);
	}
}

void
CDensity::PatchWettingNodeMirror(cgsize_t* pPnts, cgsize_t numPnts, Node::NodeValueType_t currentYNode, Node::NodeValueType_t currentXNode, 
				Node::NodeValueType_t wedgeBase, cgsize_t currentWedge)
{
	Node::NodeValueType_t xMirror(0.0), tol(0.00001), xNode(0.0), yNode(0.0);
	cgsize_t node(0);
	for(cgsize_t i=0; i < numPnts; i++)
	{
		node   =  *(pPnts + i);
		xNode  =  GetVariable()->at(0)->at(node-1)->m_xcoord;
		yNode  =  GetVariable()->at(0)->at(node-1)->m_ycoord;
		
		if(abs(yNode - currentYNode) < tol)
		{
			xMirror = wedgeBase * (1.0 + currentWedge) - (currentXNode - wedgeBase*currentWedge);
			
			if(abs(xMirror - xNode) < tol)
			{
				for(MaterialVarType_t::size_type mat_index = 0; mat_index < GetVariable()->size(); mat_index++)
				{
					GetVariable()->at(mat_index)->at(node-1)->m_wettingFlag = true;	
					//m_density.at(0)->at(node-1)->m_nodeVal         = 0.3060;
					//m_density.at(1)->at(node-1)->m_nodeVal         = 0;
				}
			}
		}
	}
}

CDensity::CDensity(CLbmCase* pCase): m_pCase(pCase), m_index(4)
{
	setBlock(pCase);
	setPatchedNodes(pCase);
	setMaterialDensity(pCase);
	if(GetPatchedNodes()->size() - 1)
	{
		setPatchedNodesDensity(pCase);
	}
	else;
}

CDensity::CDensity(CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase), m_index(4)
{
	//Add Node objects to the density container and set info
	setBlock(nvertex);
}

void 
CDensity::setPatchedNodes(CLbmCase* pCase)
{
	//Different Patch Types
	//setPrismNodes(pCase);
	//setCylinderNodes(pCase);
	setSphereNodes(pCase);
}

void
CDensity::setPrismNodes(CLbmCase* pCase)
{
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" xmax =  " << pCase->m_patchedNodes.x_max <<"  ymax = "<<pCase->m_patchedNodes.y_max 
					<<"  zmax = "<<pCase->m_patchedNodes.z_max <<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	RectPrismParameters temp = obtainRectPrismParameters(pCase);

	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
		{
				/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::cout<<" node =  " << i <<"  x = "<<GetVariable()->at(0)->at(i-1)->m_xcoord
					<<"  y = "<<GetVariable()->at(0)->at(i-1)->m_ycoord <<"  z =  "<<GetVariable()->at(0)->at(i-1)->m_zcoord<<"\n";
				std::cin.get();
				////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
				if(GetVariable()->at(0)->at(i-1)->m_xcoord <= pCase->m_patchedNodes.x_max && GetVariable()->at(0)->at(i-1)->m_xcoord >= pCase->m_patchedNodes.x_min)
				{
					if(GetVariable()->at(0)->at(i-1)->m_ycoord <= pCase->m_patchedNodes.y_max && GetVariable()->at(0)->at(i-1)->m_ycoord >= pCase->m_patchedNodes.y_min)
					{
						if(abs(temp.h) == 0)
						{
							if(abs(GetVariable()->at(0)->at(i-1)->m_zcoord - temp.c) == 0)
							{
								m_patchedNode.push_back(GetVariable()->at(0)->at(i-1)->m_nodeNum);
							}
							else;
						}
						else if(abs(temp.h) != 0)
						{
							/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
							std::cout<<" node number =  " << GetVariable()->at(0)->at(i-1)->m_nodeNum <<"\n";
							std::cin.get();
							////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
							if(GetVariable()->at(0)->at(i-1)->m_zcoord <= pCase->m_patchedNodes.z_max 
								&& GetVariable()->at(0)->at(i-1)->m_zcoord >= pCase->m_patchedNodes.z_min)
							{
								m_patchedNode.push_back(GetVariable()->at(0)->at(i-1)->m_nodeNum);
							}
						}
						else;
					}
					else;
				}
				else;
		}
}

void
CDensity::setCylinderNodes(CLbmCase* pCase)
{
	//patch cylinder
	//obtain cylinder parameters first
	CylinderParameters temp = obtainCylinderParameters(pCase);
	Node::NodeValueType_t r_square(0.0);
	//patch plane
	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
	{
		r_square = pow((GetVariable()->at(0)->at(i-1)->m_xcoord - temp.a), 2) + 
								pow((GetVariable()->at(0)->at(i-1)->m_ycoord - temp.b), 2);

		if(	r_square <= pow(temp.r_2D, 2)
				&&
				(GetVariable()->at(0)->at(i-1)->m_zcoord >= pCase->m_patchedNodes.z_min && 
				GetVariable()->at(0)->at(i-1)->m_zcoord <= pCase->m_patchedNodes.z_max)
			)
		{
			m_patchedNode.push_back(GetVariable()->at(0)->at(i-1)->m_nodeNum);
		}
		else;
	}
}

void
CDensity::setSphereNodes(CLbmCase* pCase)
{
	//patch sphere
	//obtain sphere parameters first
	SphereParameters temp = obtainSphereParameters(pCase);
	Node::NodeValueType_t r_square(0.0), center_z(0.0);
	//match the center z position
	center_z = matchZcenter(pCase, temp.c);
	//patch plane
	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
	{
		r_square = pow((GetVariable()->at(0)->at(i-1)->m_xcoord - temp.a), 2) + 
								pow((GetVariable()->at(0)->at(i-1)->m_ycoord - temp.b), 2);

		if(r_square <= pow(temp.r_2D, 2) && temp.r_3D != 0.0)
		{
			if(abs(GetVariable()->at(0)->at(i-1)->m_zcoord-temp.c) == center_z && GetVariable()->at(0)->at(i-1)->m_zcoord <= temp.c)
			{
				m_patchedNode.push_back(GetVariable()->at(0)->at(i-1)->m_nodeNum);
			}
			else;
		}
		else if(r_square <= pow(temp.r_2D, 2) && temp.r_3D == 0.0)
		{
			if(abs(GetVariable()->at(0)->at(i-1)->m_zcoord - temp.c) == 0)
			{
				m_patchedNode.push_back(GetVariable()->at(0)->at(i-1)->m_nodeNum);
			}
			else;
		}
		else;
	}
	if(temp.r_3D != 0.0)
	{
		patchZdirection(pCase,temp.r_3D, temp.a, temp.b, temp.c);
	}
	else;
}

Node::NodeValueType_t 
CDensity::matchZcenter(CLbmCase* pCase, Node::NodeValueType_t c)
{
	Node::NodeValueType_t newDiff(0.0), oldDiff(100000000), center_z(0.0);

	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
	{
		newDiff  = abs(GetVariable()->at(0)->at(i-1)->m_zcoord - c);
		center_z = (newDiff < oldDiff) ? newDiff : oldDiff;
		oldDiff  = center_z;
	}
	return center_z;
}

RectPrismParameters 
CDensity::obtainRectPrismParameters(CLbmCase* pCase)
{
	RectPrismParameters temp;
	temp.l = (pCase->m_patchedNodes.x_max - pCase->m_patchedNodes.x_min);
	temp.w = (pCase->m_patchedNodes.y_max - pCase->m_patchedNodes.y_min);
	temp.h = (pCase->m_patchedNodes.z_max - pCase->m_patchedNodes.z_min);
	temp.a    = pCase->m_patchedNodes.x_min + temp.l/2.0;
	temp.b    = pCase->m_patchedNodes.y_min + temp.w/2.0;
	temp.c    = pCase->m_patchedNodes.z_min + temp.h/2.0;
	return temp;
}
CylinderParameters 
CDensity::obtainCylinderParameters(CLbmCase* pCase)
{
	CylinderParameters temp;
	temp.r_2D = (pCase->m_patchedNodes.x_max - pCase->m_patchedNodes.x_min)/2.0;
	temp.h    = pCase->m_patchedNodes.z_max - pCase->m_patchedNodes.z_min;
	temp.a    = pCase->m_patchedNodes.x_min + temp.r_2D;
	temp.b    = pCase->m_patchedNodes.y_min + temp.r_2D;
	temp.c    = pCase->m_patchedNodes.z_min + temp.h/2;
	return temp;
}
SphereParameters 
CDensity::obtainSphereParameters(CLbmCase* pCase)
{
	SphereParameters temp;
	temp.r_2D = (pCase->m_patchedNodes.x_max - pCase->m_patchedNodes.x_min)/2.0;
	temp.r_3D = (pCase->m_patchedNodes.z_max - pCase->m_patchedNodes.z_min)/2.0;
	temp.a    = pCase->m_patchedNodes.x_min + temp.r_2D;
	temp.b    = pCase->m_patchedNodes.y_min + temp.r_2D;
	temp.c    = pCase->m_patchedNodes.z_min + temp.r_3D;
	return temp;
}

void 
CDensity::patchZdirection(CLbmCase* pCase, Node::NodeValueType_t r_3D, Node::NodeValueType_t a, Node::NodeValueType_t b, Node::NodeValueType_t c)
{
	cgsize_t node(0);
	PatchNodeType_t::size_type end = GetPatchedNodes()->size();

	for(PatchNodeType_t::size_type j=1; j <= end; j++)
	{
		node = GetPatchedNodes()->at(j-1);
		Node::NodeValueType_t r_square = pow((GetVariable()->at(0)->at(node-1)->m_xcoord - a), 2) + 
																			pow((GetVariable()->at(0)->at(node-1)->m_ycoord - b), 2);
		Node::NodeValueType_t Z_u = c + sqrt(pow(r_3D,2) - r_square);
		Node::NodeValueType_t Z_l = c - sqrt(pow(r_3D,2) - r_square);

		for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
			{
				if(  m_density.at(0)->at(i-1)->m_xcoord == m_density.at(0)->at(node-1)->m_xcoord 
					&& m_density.at(0)->at(i-1)->m_ycoord == m_density.at(0)->at(node-1)->m_ycoord)
				{
					if(m_density.at(0)->at(i-1)->m_zcoord <= Z_u && m_density.at(0)->at(i-1)->m_zcoord >= Z_l)
					{
						m_patchedNode.push_back(m_density.at(0)->at(i-1)->m_nodeNum);
					}
					else;
				}
			else;
		}
	}
}

CDensity::PatchNodeType_t*
CDensity::GetPatchedNodes()
{
	return &m_patchedNode;
}


CVariable::MaterialVarType_t* 
CDensity::variable() 
{
	return &m_density;
}

CLbmCase* 
CDensity::caseInfo() 
{
	return m_pCase;
}

cgsize_t 
CDensity::index() 
{
	return m_index;
}

CDensity::~CDensity()
{

}

// Temp for multiphase
CTemp::CTemp(CLbmCase* pCase): m_pCase(pCase)
{
	m_index = 0;
	setBlock(pCase);
}

CTemp::CTemp(CLbmCase* pCase, cgsize_t index): m_pCase(pCase)
{
	//Add Node objects to the psi container and set info
	m_index = index;
	setBlock(pCase);
}

CVariable::MaterialVarType_t* 
CTemp::variable() 
{
	return &m_temp;
}

CLbmCase* 
CTemp::caseInfo() 
{
	return m_pCase;
}

cgsize_t 
CTemp::index() 
{
	return m_index;
}

CTemp::~CTemp()
{

}

