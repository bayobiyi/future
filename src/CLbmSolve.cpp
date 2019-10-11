
#include "CLbmBlock.h"
#include "CLbmSolve.h"
#include <time.h>
#include <vector>
#include <cmath>
#include <iostream>
#include <omp.h>



//LEVEL ONE:initialized in main()
//==============================================================================================//
//Getters
Domain::MaterialType_t*
Domain::GetDomainMaterial()
{
	return DomainMaterial();
}

Domain::VariableType_t*
Domain::GetDomainVariables()
{
	return DomainVariables();
}

Domain::VariableType_t*
Domain::GetDomainTempVariables()
{
	return DomainTempVariables();
}

Domain::PatchNodeType_t*
Domain::GetDomainSolidObjectNodes()
{
	return DomainSolidObjectNodes();
}

Domain::MicroPillarsType_t*
Domain::GetDomainMicroPillars()
{
	return DomainMicroPillars();
}
//initialization methods
void
Domain::mapBounds(CLbmCase* pCase)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" start mapping bounds \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), node(0);
	cgsize_t interior(2);
	for(cgsize_t j=1; j <= pCase->m_boco.m_nBocos; j++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (j-1)) != interior)
			{
				//first obtain nodes on boundary
				pPnts    = *(pCase->m_boco.m_pnts  + (j-1));
				numPnts  = *(pCase->m_boco.m_npnts + (j-1));

				for(cgsize_t i =0; i < numPnts; i++)
				{
						node = *(pPnts + i);
						//map cardinal direction
						mapCardinalDirection(node);
						//map isSolid, halfBB, fullBB
						mapSolidNodes(pCase, node, j);
						//map isNotSolidExtrap
						mapExtrapNodes(pCase, node, j);
						//map inward/buried nodes
				}
			}
	}
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"finished mapping bound\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}


void
Domain::mapBoundType(CLbmCase* pCase)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" start mapping bounds \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), node(0);
	cgsize_t interior(2);
	for(cgsize_t j=1; j <= pCase->m_boco.m_nBocos; j++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (j-1)) != interior)
			{
				//first obtain nodes on boundary
				pPnts    = *(pCase->m_boco.m_pnts  + (j-1));
				numPnts  = *(pCase->m_boco.m_npnts + (j-1));

				for(cgsize_t i =0; i < numPnts; i++)
				{
						node = *(pPnts + i);
						mapInwardORBuriedNodes(node);
						//map start of streaming
				}
			}
	}
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"finished mapping bound type \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

void 
Domain::mapPeriodicCorners(CLbmCase* pCase, cgsize_t mat_index, cgsize_t node)
{
	cgsize_t  one(1), neighbor(0);

	if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::WSB)
	{
		neighbor = fetchPeriodicNeighbor(pCase, Node::ENT, mat_index, node);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::WNB)
	{
		neighbor = fetchPeriodicNeighbor(pCase, Node::EST, mat_index, node);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::ESB)
	{
		neighbor = fetchPeriodicNeighbor(pCase, Node::WNT, mat_index, node);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::ENB)
	{
		neighbor = fetchPeriodicNeighbor(pCase, Node::WST, mat_index, node);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else;
}

void 
Domain::mapPeriodicEdges(CLbmCase* pCase, cgsize_t mat_index, cgsize_t node)
{
	cgsize_t  one(1), neighbor(0);
	Node::NodeValueType_t coord1(0.0);
	cgsize_t xcoord(0), ycoord(1), zcoord(2);

	if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::WB)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::ET, mat_index, node, coord1, ycoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::SB)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::NT, mat_index, node, coord1, xcoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::EB)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::WT, mat_index, node, coord1, ycoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::NB)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::ST, mat_index, node, coord1, xcoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::WS)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::EN, mat_index, node, coord1, zcoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::ES)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::WN, mat_index, node, coord1, zcoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else;
}

void 
Domain::mapPeriodicSurfaces(CLbmCase* pCase, cgsize_t mat_index, cgsize_t node)
{
	cgsize_t  one(1), neighbor(0);
	Node::NodeValueType_t coord1(0.0), coord2(0.0);
	cgsize_t xcoord(0), ycoord(1), zcoord(2);
	if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::WEST)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord;
		coord2 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::EAST, mat_index, node, coord1,coord2, ycoord, zcoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}

	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::SOUTH)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord;
		coord2 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::NORTH, mat_index, node, coord1,coord2, xcoord, zcoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	
	else if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM)
	{
		coord1 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord;
		coord2 = GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord;
		neighbor = fetchPeriodicNeighbor(pCase, Node::TOP, mat_index, node, coord1,coord2, xcoord, ycoord);
		copyPeriodicNeighbor(pCase, node, neighbor);
	}
	else;
}

void
Domain::mapPeriodicNeighbor(CLbmCase* pCase)
{
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), node(0), zero(0);
	cgsize_t interior(2);
	for(cgsize_t i=0; i < pCase->m_boco.m_nBocos; i++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + i) != interior)
		{
			//first obtain nodes on boundary
			pPnts    = *(pCase->m_boco.m_pnts  + i);
			numPnts  = *(pCase->m_boco.m_npnts + i);

			for(cgsize_t j =0; j < numPnts; j++)
			{
				node = *(pPnts + j);
				mapPeriodicCorners(pCase, zero, node);
				mapPeriodicEdges(pCase, zero, node);
				mapPeriodicSurfaces(pCase, zero, node);
			}
		}
		else;
	}
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"finished mapping periodic bound\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

void 
Domain::copyPeriodicNeighbor(CLbmCase* pCase, cgsize_t node, cgsize_t neighbor)
{
	for(MaterialType_t::size_type sub=0 ; sub < GetDomainMaterial()->size(); sub++)
	{
		for(PdfDomain::LatticeType_t::size_type i =1; i < GetDomainMaterial()->at(sub)->GetLatticePdf()->size(); i++)
		{
			GetDomainMaterial()->at(sub)->GetLatticePdf()->at(i)->pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum = neighbor;
		}
	}
}

cgsize_t 
Domain::fetchPeriodicNeighbor(CLbmCase* pCase, Node::LBMBOUND bound, cgsize_t mat_index, cgsize_t periodicNode)
{
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), node(0), one(1), temp(0);
	cgsize_t periodic(0), interior(2);

	for(cgsize_t i=0; i < pCase->m_boco.m_nBocos; i++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + i) != interior && *(pCase->m_boco.m_bocoUserDef + i)  == periodic)
		{
				//first obtain nodes on boundary
			pPnts    = *(pCase->m_boco.m_pnts  + i);
			numPnts  = *(pCase->m_boco.m_npnts + i);

			for(cgsize_t j =0; j < numPnts; j++)
			{
				node = *(pPnts + j);
				if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == bound)
				{
					temp = node;
					copyPeriodicNeighbor(pCase, temp, periodicNode);

				}
				else;
			}
		}
	}			
	return temp;
}
cgsize_t 
Domain::fetchPeriodicNeighbor(CLbmCase* pCase, Node::LBMBOUND bound, cgsize_t mat_index, cgsize_t periodicNode, Node::NodeValueType_t coord1, cgsize_t index)
{
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), node(0), one(1), temp(0);
	cgsize_t periodic(0), interior(2);
	cgsize_t xcoord(0), ycoord(1), zcoord(2);
	Node::NodeValueType_t tolerance(0.0001);
	for(cgsize_t i=0; i < pCase->m_boco.m_nBocos; i++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + i) != interior && *(pCase->m_boco.m_bocoUserDef + i)  == periodic)
		{
			//first obtain nodes on boundary
			pPnts    = *(pCase->m_boco.m_pnts  + i);
			numPnts  = *(pCase->m_boco.m_npnts + i);
			
			for(cgsize_t j =0; j < numPnts; j++)
			{
				node = *(pPnts + j);
				if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == bound)
				{
					if(index == xcoord)
					{
						if(abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord - coord1) <= tolerance)
						{
							temp = node;
							copyPeriodicNeighbor(pCase, temp, periodicNode);
						}
						else;
					}
					else if(index == ycoord)
					{
						if(abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord - coord1) <= tolerance)
						{
							temp = node;
							copyPeriodicNeighbor(pCase, temp, periodicNode);
						}
						else;
					}
					else if(index == zcoord)
					{
						if(abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord - coord1) <= tolerance)
						{
							temp = node;
							copyPeriodicNeighbor(pCase, temp, periodicNode);
						}
						else;
					}
					else;
				}
				else;
			}
		}
	}			
	return temp;
}

cgsize_t 
Domain::fetchPeriodicNeighbor(CLbmCase* pCase, Node::LBMBOUND bound, cgsize_t mat_index, cgsize_t periodicNode, 
Node::NodeValueType_t coord1, Node::NodeValueType_t coord2, cgsize_t index1, cgsize_t index2)
{
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), node(0), one(1), temp(0);
	cgsize_t periodic(0), interior(2);
	cgsize_t xcoord(0), ycoord(1), zcoord(2);
	Node::NodeValueType_t tolerance(0.0001);

	for(cgsize_t i=0; i < pCase->m_boco.m_nBocos; i++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + i) != interior && *(pCase->m_boco.m_bocoUserDef + i)  == periodic)
		{
			//first obtain nodes on boundary
			pPnts    = *(pCase->m_boco.m_pnts  + i);
			numPnts  = *(pCase->m_boco.m_npnts + i);

			for(cgsize_t j =0; j < numPnts; j++)
			{
				node = *(pPnts + j);
				if(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_bound == bound)
				{
					if(index1 == xcoord && index2 == ycoord)
					{
						if((abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord - coord1) <= tolerance)
								&&
						(abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord - coord2) <= tolerance))
						{
							temp = node;
							copyPeriodicNeighbor(pCase, temp, periodicNode);
						}
						else;
					}
					else if(index1 == ycoord && index2 == zcoord)
					{
						if((abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord - coord1) <= tolerance)
								&&
						(abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord - coord2) <= tolerance))
						{
							temp = node;
							copyPeriodicNeighbor(pCase, temp, periodicNode);
						}
						else;
					}
					else if(index1 == xcoord && index2 == zcoord)
					{
						if((abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord - coord1) <= tolerance)
								&&
						(abs(GetDomainMaterial()->at(mat_index)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord - coord2) <= tolerance))
						{
							temp = node;
							copyPeriodicNeighbor(pCase, temp, periodicNode);
						}
						else;
					}
					else;
				}
				else;
			}
		}
	}				
	return temp;
}


void 
Domain::mapSolidNodes(CLbmCase* pCase, cgsize_t node, cgsize_t bry)
{
	cgsize_t Zhou_He_FullBB(0), extrapolate(3);
	
	if((*(pCase->m_boco.m_bocoType + (bry-1))) == BCWall)
	{
		/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" map solid nodes \n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
			{
				//map solid nodes in variables
				copyToOthers(node, i);

				for(PdfDomain::LatticeType_t::size_type j =0; j < GetDomainMaterial()->at(i)->GetLatticePdf()->size(); j++)
				{
					//overule any prior mapping
					GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_isSolid  = true;

					if(*(pCase->m_boco.m_bocoWall + (bry-1)) == Zhou_He_FullBB)
						GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_isFullBB = true;

					else if(*(pCase->m_boco.m_bocoWall + (j-1)) == extrapolate )
						GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_isExtrapSolid = true;

					else;

					GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_boundFlagUpdated = true;
				}
			}
	}
}

void 
Domain::setPatchedSolidObjects(CLbmCase* pCase)
{
	cgsize_t node(0);
	for(cgsize_t i = 1; i <= pCase->m_grid.m_NVertex; i++)
	{
		// if prism object
		node = fetchPrismNode(pCase, i);
		//node = fetchSphereNode(pCase, i);
		//node = fetchCylinderNode(pCase, i);
		if(node)
			GetDomainSolidObjectNodes()->push_back(node);
		else;
		
	}
	
	//update state of solid nodes in domain
	for(PatchNodeType_t::size_type j=1; j <= GetDomainSolidObjectNodes()->size(); j++)
	{
		updateSolidNodes(GetDomainSolidObjectNodes()->at(j-1));
	}
	
}

cgsize_t 
Domain::fetchPrismNode(CLbmCase* pCase, cgsize_t node)
{
	cgsize_t zero(0), one(1), solidNode(0);

	Node::NodeValueType_t x_min(0.0), x_max(80.0), y_min(91.0), y_max(133.0), z_min(19.0), z_max(61.0);

	if(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord <= x_max 
	&& GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord >= x_min)
	{
		if(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord <= y_max 
		&& GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord >= y_min)
		{
			if(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord <= z_max 
			&& GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord >= z_min)
			{
				solidNode = node;
			}
			else;
		}
		else;
	}
	else;
	return solidNode;	
}

cgsize_t 
Domain::fetchSphereNode(CLbmCase* pCase, cgsize_t node)
{
	cgsize_t zero(0), one(1), solidNode(0);

	Node::NodeValueType_t a(40.0), b(122.0), c(40.0);

	Node::NodeValueType_t sphereRadius(10.0), r_square(0.0); 

	r_square = pow((GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord - a), 2) + 
		   pow((GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord - b), 2) +
		   pow((GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord - c), 2);

	if(r_square <= pow(sphereRadius,2))
	{
		solidNode = node;
	}
	else;
	return solidNode;
}

cgsize_t 
Domain::fetchCylinderNode(CLbmCase* pCase, cgsize_t node)
{
	cgsize_t zero(0), one(1), solidNode(0);
	Node::NodeValueType_t r_square(0.0), cyldRadius(21.0);
	Node::NodeValueType_t b(112.0), c(40.0), x_min(0.0), x_max(80.0);
	
	r_square = pow((GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord - b), 2) + 
		   pow((GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord - c), 2);

	if(
		r_square < pow(cyldRadius,2)
		&&
		(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord   <= x_max 
		&& GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord >= x_min)
	  )
	{
			solidNode = node;
	}
	else;
	
	return solidNode;
}

// Micropillars Patch
void 
Domain::setPatchedMicroPillars(CLbmCase* pCase)
{
	Node::NodeValueType_t coord1(0.0), coord2(0.0), coord3(0.0), size1(100.0), size2(300.0), height(3.0), side(4.0);
	Node::NodeValueType_t coord1Min(0.0), coord2Min(0.0), coord3Min(0.0);
	Node::NodeValueType_t coord1Max(0.0), coord2Max(0.0), coord3Max(0.0);
	cgsize_t counter(0), coord1_counter(0), coord2_counter(0), node(0);

	Node::NodeValueType_t coord1_interval =6.0;
	Node::NodeValueType_t coord2_interval =6.0;

	while (coord1Min < size1)
	{
		coord1Max =   (coord1Min + side);
		
		while (coord2Min < size2)
		{
			coord2Max =   (coord2Min + side);

			coord3Min = 0.0;
			coord3Max = height;
			
			counter += 1;
			GetDomainMicroPillars()->push_back(new PillarNodeType_t);

			for(cgsize_t node = 1; node <= pCase->m_grid.m_NVertex; node++)
			{
				if(nodeInRange(node, coord1Min, coord2Min, coord3Min, coord1Max, coord2Max, coord3Max))
				{
					GetDomainMicroPillars()->at(counter-1)->push_back(node);
				}
				
			}
			coord2Min += (side + coord2_interval);
		} 
		coord2Min = 0.0;
		coord2Max = 0.0;
		coord1Min += (side + coord1_interval);  	
	}
	//update state of solid nodes in domain
	for(MicroPillarsType_t::size_type i=1; i <= GetDomainMicroPillars()->size(); i++)
	{
		for(PillarNodeType_t::size_type j=1; j <= GetDomainMicroPillars()->at(i-1)->size(); j++)
		{
			node = GetDomainMicroPillars()->at(i-1)->at(j-1);
			updatePillarNodesWetting(node , GetDomainMicroPillars()->at(i-1));
			updateSolidNodes(node);	
		}
	}
}

bool 
Domain::nodeInRange(cgsize_t node, Node::NodeValueType_t coord1Min, Node::NodeValueType_t coord2Min, Node::NodeValueType_t coord3Min, 
		    Node::NodeValueType_t coord1Max, Node::NodeValueType_t coord2Max, Node::NodeValueType_t coord3Max)
{
	cgsize_t zero(0), one(1);
	bool status(false);
	if(
		(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord   <= coord1Max &&
		 GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord   >= coord1Min  )
		&&
		(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord   <= coord2Max &&
		 GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord   >= coord2Min  )
		&&
		(GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord   <= coord3Max &&
		 GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord   >= coord3Min  )
	)
	{
		status = true;
		/*std::cout<<"node = "<<node<<" is part of a pillar \n";
		std::cout<<"x-coord = "<<GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_xcoord<<"\n";
		std::cout<<"y-coord = "<<GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_ycoord<<"\n";
		std::cout<<"z-coord = "<<GetDomainMaterial()->at(zero)->GetLatticePdf()->at(one)->pdf()->GetPdfunction()->at(node-1)->m_zcoord<<"\n";
		std::cin.get();*/
	}
	return status;
}


void 
Domain::updateSolidNodes(cgsize_t node)
{
	//cgsize_t rho(4);
	for(MaterialType_t::size_type j=0 ; j < GetDomainMaterial()->size(); j++)
	{
		//fix density and velocity on solid nodes
		copyToOthers(node, j, 0.0, 0.0);
		//copy isSolid status to all variables for current node
		copyToOthers(node, j);

		for(PdfDomain::LatticeType_t::size_type k =0; k < GetDomainMaterial()->at(j)->GetLatticePdf()->size(); k++)
		{
			GetDomainMaterial()->at(j)->GetLatticePdf()->at(k)->pdf()->GetPdfunction()->at(node-1)->m_isSolid  = true;
			//GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_nodeVal = 0.4;
		}
			
	}
}

void 
Domain::updatePillarNodesWetting(cgsize_t node, PillarNodeType_t* pillar)
{
	cgsize_t rho(4);
	Node::NodeValueType_t coord1(0.0), coord2(0.0);

	for(MaterialType_t::size_type j=0 ; j < GetDomainMaterial()->size(); j++)
	{
		if(GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_wettingFlag == true)
		{
			//compare all other nodes coord1 and coord2 to that of current node 
			coord1 = GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_xcoord;
			coord2 = GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_ycoord;
			scanPillar(coord1, coord2, pillar);
		}
		/*else
		{
			GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_nodeVal = 0.02;
		}*/			
	}
}

void 
Domain::scanPillar(Node::NodeValueType_t coord1, Node::NodeValueType_t coord2, PillarNodeType_t* pillar)
{
	cgsize_t rho(4), node(0);
	Node::NodeValueType_t tolerance(0.0001);

	for(PillarNodeType_t::size_type j=1; j <= pillar->size(); j++)
	{
		node = pillar->at(j-1);
		for(MaterialType_t::size_type j=0 ; j < GetDomainMaterial()->size(); j++)
		{
			if((abs(GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_xcoord - coord1) <= tolerance) 
			&& (abs(GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_ycoord - coord2) <= tolerance)) 
			{
				GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_wettingFlag = true;
				//GetDomainVariables()->at(rho)->GetVariable()->at(j)->at(node-1)->m_nodeVal = 0.4;
			}		
		}
			
	}
}

void
Domain::copyToOthers(cgsize_t node, cgsize_t mat_index)
{
	//copy isSolid to other variables (I CHANGED INDEX FROM K=1 T0 K=0 ALTHOUGH XVEL HAVS ALREADY BEEN UPDATED
	for(Domain::VariableType_t::size_type k=0; k < GetDomainVariables()->size(); k++)
	{
		GetDomainVariables()->at(k)->GetVariable()->at(mat_index)->at(node-1)->m_isSolid = true;
	}
}
void 
Domain::copyToOthers(cgsize_t node, cgsize_t mat_index, Node::NodeValueType_t density, Node::NodeValueType_t velocity)
{
	cgsize_t rho(4), xvel(0), yvel(1), zvel(2);
	for(Domain::VariableType_t::size_type k=0; k < GetDomainVariables()->size(); k++)
	{
		GetDomainVariables()->at(k)->GetVariable()->at(mat_index)->at(node-1)->m_isSolid = true;
		if(k == rho)
			GetDomainVariables()->at(k)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = density;
		else if (k == xvel || k == yvel || k == zvel )
			GetDomainVariables()->at(k)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = velocity;
	}

}
void 
Domain::mapExtrapNodes(CLbmCase* pCase, cgsize_t node, cgsize_t bry)
{

	if(*(pCase->m_boco.m_bocoExtrapolationBC + (bry-1)))
	{
		/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" map extrapolation nodes \n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
			{
				for(PdfDomain::LatticeType_t::size_type j =0; j < GetDomainMaterial()->at(i)->GetLatticePdf()->size(); j++)
				{
					if(GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_boundFlagUpdated == false)
					{
						GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_isExtrapNotSolid  = true;

						GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_boundFlagUpdated  = true;
					}
					else;
				}
			}
	}
}

void 
Domain::mapInwardORBuriedNodes(cgsize_t node)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" map inward or buried nodes \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
	{
		for(PdfDomain::LatticeType_t::size_type j =1; j < GetDomainMaterial()->at(i)->GetLatticePdf()->size(); j++)
		{
			if(GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0)
			{
				GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->setNodeBType(node);
			}
			else;
		}
	}
}


void
Domain::copyToOthers(cgsize_t node, Node::LBMBOUND bound)
{
	//copy bound to pdfs 
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
	{
		//copy bound to each pdf block
		for(PdfDomain::LatticeType_t::size_type j =0; j < GetDomainMaterial()->at(i)->GetLatticePdf()->size(); j++)
		{
			GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(node-1)->m_bound = bound;
		}
		//copy bound to other variables
		for(Domain::VariableType_t::size_type k=1; k < GetDomainVariables()->size(); k++)
		{
			GetDomainVariables()->at(k)->GetVariable()->at(i)->at(node-1)->m_bound = bound;
		}
	}
}

void 
Domain::initializeVars(CLbmCase*  pCase,  VariableType_t* variables)
{
	cgsize_t u_overallx(5), u_overally(6), u_overallz(7), pressure(8), u_primex(8), u_primey(9), u_primez(10);
		//Initialize domain variables
		variables->push_back(new CXVelocity(pCase));   //xvelocity index = 0
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in cxvelocity\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CYVelocity(pCase));   //yvelocity index = 1
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in cyvelocity\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CZVelocity(pCase));   //zvelocity index = 2
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in czvelocity\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CVelocityMag(pCase)); //velocity magnitude index = 3
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in cvelmag\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CDensity(pCase));	    //density index = 4
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in cdensity\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CTemp(pCase, u_overallx));        //u_overallx  index = 5
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in uoverall_x\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CTemp(pCase, u_overally));        //u_overally  index = 6
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in uoverall_y\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CTemp(pCase, u_overallz));        //u_overallz  index = 7
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in uoverall_z\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		variables->push_back(new CTemp(pCase, pressure));        //presure  index = 8
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" pushed in uoverall_z\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		if( pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
		{
			variables->push_back(new CTemp(pCase, u_primex));        //u_primex  index = 8
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" pushed in uprime_x\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			variables->push_back(new CTemp(pCase, u_primey));        //u_primey  index = 9
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" pushed in uprime_y\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			variables->push_back(new CTemp(pCase, u_primez));        //u_primez  index = 10
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" pushed in uprime_z\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/	
		}
}
void 
Domain::initializeTemps(CLbmCase* pCase,  VariableType_t* tempVariables)
{
	tempVariables->push_back(new CTemp(pCase));        //psi			index = 0
	tempVariables->push_back(new CTemp(pCase));        //ftempx			index = 1
	tempVariables->push_back(new CTemp(pCase));        //ftempy			index = 2
	tempVariables->push_back(new CTemp(pCase));        //ftempz			index = 3
	tempVariables->push_back(new CTemp(pCase));        //fsurftempx 		index = 4
	tempVariables->push_back(new CTemp(pCase));        //fsurftempy 		index = 5
	tempVariables->push_back(new CTemp(pCase));        //fsurftempz 		index = 6
	tempVariables->push_back(new CTemp(pCase));        //phi               		index = 7

	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" pushed all temps\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"domainTempVariables contains = " << tempVariables->size() << " elements\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void 
Domain::Update_variables(CLbmCase* pCase)
{
	for(VariableType_t::size_type i=0 ; i < GetDomainVariables()->size(); i++)
		GetDomainVariables()->at(i)->Update_variables(pCase, i);
}

void
Domain::initialize(CLbmCase* pCase, Domain& domainVariables)
{
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
		GetDomainMaterial()->at(i)->initialize(pCase, domainVariables,i);
}
//process methods

void 
Domain::Lbm_write(std::fstream& outputFile)
{
	/*for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
		GetDomainMaterial()->at(i)->Lbm_write(outputFile);
	for(VariableType_t::size_type i=0 ; i < GetDomainVariables()->size(); i++)
		GetDomainVariables()->at(i)->Lbm_write(outputFile);*/
}

void 
Domain::Lbm_read(std::fstream& inputFile)
{
	/*for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
		GetDomainMaterial()->at(i)->Lbm_read(inputFile);
	for(VariableType_t::size_type i=0 ; i < GetDomainVariables()->size(); i++)
		GetDomainVariables()->at(i)->Lbm_read(inputFile);*/
}

void
Domain::collide(CLbmCase* pCase, Domain& domainVariables)
{
	//double starttime;
	//double endtime;
	//starttime = omp_get_wtime();
	//First update interparticle forces
	if(pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" in multiphase about to compute interaction forces \n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//MCMP----------------------------------------------------------------------------------------------------------------------------
		computeInteractionForces(pCase, domainVariables);
		//computeSurfaceForces(pCase, domainVariables);
		//endtime = omp_get_wtime();
		//printf("Precollsion Work took %f sec. time.\n", endtime-starttime);

		for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i+=2)
		{
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" now perform collision for material = "<< i <<"\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			GetDomainMaterial()->at(i)->collide(pCase, domainVariables, i);
			GetDomainMaterial()->at(i+1)->collide(pCase, domainVariables, i+1);
		}
		//MCMP-----------------------------------------------------------------------------------------------------------------------------
	}
	
	else
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" in single component about to compute interaction forces \n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//SCMP---------------------------------------------------------------------------------------------------------------------------
		computeInteractionForcesSCMP(pCase, domainVariables);
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" finished computing interaction force \n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
		{
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" now perform collision for material = "<< i <<"\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			GetDomainMaterial()->at(i)->collide(pCase, domainVariables, i);
			/*******************do this if MRT*************************/
			GetDomainMaterial()->at(i)->momentToVelocitySpace(pCase, domainVariables);
			/**********************************************************/
		}
		//SCMP----------------------------------------------------------------------------------------------------------------------------
	}
}
void
Domain::post_collide(CLbmCase* pCase, Domain& domainVariables)
{
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now perform post collision for material = "<< i <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetDomainMaterial()->at(i)->post_collide(pCase, domainVariables,i);
	}
}
void
Domain::stream(CLbmCase* pCase, Domain& domainVariables)
{
	//SCMP---------------------------------------------------------------------------------------------------------------
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now perform streaming for material= "<< i <<" \n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetDomainMaterial()->at(i)->stream(pCase, domainVariables);
	}
	//SCMP---------------------------------------------------------------------------------------------------------------

	/*MCMP---------------------------------------------------------------------------------------------------------------
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i+=2)
	{
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now perform streaming for material= "<< i <<" \n";
		std::cin.get();
		///////////////////////////////////////////////////////////////////////////////////////////////////////////
		GetDomainMaterial()->at(i)->stream(pCase, domainVariables);
		GetDomainMaterial()->at(i+1)->stream(pCase, domainVariables);
	}
	//MCMP---------------------------------------------------------------------------------------------------------------*/
}
void
Domain::post_stream(CLbmCase* pCase, Domain& domainVariables)
{
	for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now perform post stream for material= "<< i <<" \n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetDomainMaterial()->at(i)->post_stream(pCase, domainVariables, i);
	}
}

void
Domain::updateVariablesSglComp(CLbmCase* pCase, Domain& domainVariables)
{
	cgsize_t single(0), node(0), node1(0), node2(0), node3(0);
	cgsize_t xvel(0), yvel(1), zvel(2), rho(4);
	cgsize_t end_nodes = pCase->m_grid.m_NVertex;

	PdfDomain::LatticeType_t::size_type f;
	PdfDomain::LatticeType_t::size_type end_pdf = GetDomainMaterial()->at(single)->GetLatticePdf()->size();

#pragma omp parallel private(f)
{//----------------------------start of parallel region----------------------------------------------------------
#pragma omp for private(node)
	for(node = 1; node <= end_nodes; node++)
	{
		//store new to old for fluid variables for each component 
		if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(single)->at(node-1)->m_isSolid == false)
		{
			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(single)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(single)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(single)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(single)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(single)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(single)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(single)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(single)->at(node-1)->m_nodeVal;

			

			//reset to zero before summation 
			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(single)->at(node-1)->m_nodeVal     = 0.0;
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(single)->at(node-1)->m_nodeVal    = 0.0;
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(single)->at(node-1)->m_nodeVal    = 0.0;
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(single)->at(node-1)->m_nodeVal    = 0.0;
		}
		else;
	}
	
	for(f =0; f < end_pdf; f++)
	{
#pragma omp for private(node1)
		for(node1 = 1; node1 <= end_nodes; node1++)
		{
			
			GetDomainMaterial()->at(single)->updateAllSummations(pCase, domainVariables,single, node1, f);
		}
	}
#pragma omp for private(node2)
	for(node2 = 1; node2 <= end_nodes; node2++)
	{
		computeOverallVarsSgl(pCase, domainVariables,   node2);
	}
#pragma omp for private(node3)
	for(node3 = 1; node3 <= end_nodes; node3++)
	{
		computeIndividualVars(pCase, domainVariables, single, node3);
	}
}//---------------------------end of parallel region----------------------------------------------------------------------------

}

void
Domain::updateVariables(CLbmCase* pCase, Domain& domainVariables)
{
	cgsize_t primary(0), secondary(1), node(0), node1(0), node2(0), node3(0);
	cgsize_t xvel(0), yvel(1), zvel(2), rho(4), uoverall_x(5), uoverall_y(6), uoverall_z(7);
	cgsize_t end_nodes = pCase->m_grid.m_NVertex;

	PdfDomain::LatticeType_t::size_type f;
	PdfDomain::LatticeType_t::size_type end_pdf = GetDomainMaterial()->at(primary)->GetLatticePdf()->size();

#pragma omp parallel private(f)
{/*----------------------------start of parallel region----------------------------------------------------------*/
#pragma omp for private(node)
	for(node = 1; node <= end_nodes; node++)
	{
		//store new to old for fluid variables for each component 
		if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_isSolid == false)
		{
			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;

			//store new to old for overall fluid velocities 
			domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

			//reset to zero before summation 
			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal     = 0.0;
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal    = 0.0;
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal    = 0.0;
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal    = 0.0;

			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal     = 0.0;
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal    = 0.0;
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal    = 0.0;
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal    = 0.0;
		}
		else;
	}
	
	for(f =0; f < end_pdf; f++)
	{
#pragma omp for private(node1)
		for(node1 = 1; node1 <= end_nodes; node1++)
		{
			
			GetDomainMaterial()->at(primary)->updateAllSummations(pCase, domainVariables,primary, node1, f);
			GetDomainMaterial()->at(secondary)->updateAllSummations(pCase, domainVariables,secondary, node1, f);
		}
	}

	//if multiphase compute composite velocity
	if(pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
#pragma omp for private(node2)
		for(node2 = 1; node2 <= end_nodes; node2++)
		{
			computeCompositeVars(pCase, domainVariables, node2);
			computeOverallVars(pCase, domainVariables,   node2);
		}
	}
	else;
#pragma omp for private(node3)
	for(node3 = 1; node3 <= end_nodes; node3++)
	{
		computeIndividualVars(pCase, domainVariables, primary, node3);
		computeIndividualVars(pCase, domainVariables, secondary, node3);
	}
}/*---------------------------end of parallel region----------------------------------------------------------------------------*/

}

void
Domain::writeSolution(CLbmCase* pCase)
{
	//write to solution container
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"# of fields in solution container = "<<pCase->m_sol.m_sol.size() <<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"# of materials in field 5  container = "<<pCase->m_sol.m_sol.at(5)->size()<<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"# of nodes in material 1 of field 1  container = "<<pCase->m_sol.m_sol.at(0)->at(0)->size() <<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"size of domainVariables container = "<<GetDomainVariables()->size() <<"\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	//Store all physical variables for export first 
	for(Domain::VariableType_t::size_type i=0; i<GetDomainVariables()->size(); i++)
	{
		for(CVariable::MaterialVarType_t::size_type j=0; j<GetDomainVariables()->at(i)->GetVariable()->size(); j++)
		{
			/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<"field = "<<i<<"contains material = "<<GetDomainVariables()->at(i)->GetVariable()->size();
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			for(CVariable::NodeType_t::size_type k=0; k<GetDomainVariables()->at(i)->GetVariable()->at(j)->size(); k++)
			{
				pCase->m_sol.m_sol.at(i)->at(j)->at(k) = GetDomainVariables()->at(i)->GetVariable()->at(j)->at(k)->m_nodeVal;
				pCase->m_sol.m_solOld.at(i)->at(j)->at(k) = GetDomainVariables()->at(i)->GetVariable()->at(j)->at(k)->m_oldNodeVal;	
				/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::cout<<" node = "<<k<<" node value = "<<pCase->m_sol.m_sol.at(i)->at(j)->at(k)<<"\n";
				std::cin.get();
				////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			}
		}
	}
}

void 
Domain::computeInteractionForcesSCMP(CLbmCase* pCase, Domain& domainVariables)
{
	cgsize_t single(0), xVel(0), fx(1), fy(2), fz(3), node(0), node1(0);
	Node::NodeValueType_t G_sgl(-1.0);
	PdfDomain::LatticeType_t::size_type f;
	PdfDomain::LatticeType_t::size_type end_f = GetDomainMaterial()->at(single)->GetLatticePdf()->size();
	cgsize_t end = pCase->m_grid.m_NVertex;
	/*std::cout<<"Now computing interaction force for SCMP\n";
	std::cout<<"number of nodes = "<<end<<"\n";*/
	
#pragma omp parallel private(f)
{/*-----------------------start of parallel region-----------------------------------------------------------------------------*/
#pragma omp for private(node)
	for(node = 1; node <= end; node++)
	{
		domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(single)->at(node-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(single)->at(node-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(single)->at(node-1)->m_nodeVal = 0.0;

		//compute psi
		computePsi(pCase, domainVariables, node, single, G_sgl);
	}
	for(f =1; f < end_f; f++)
	{
#pragma omp for private(node1)
		for(node1 = 1; node1 <= end; node1++)
		{
			cgsize_t psi(0);
			
			if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(single)->at(node1-1)->m_isSolid == false)
			{
				computeSummationIntra(pCase, domainVariables, f, node1, single, G_sgl, psi);
			}
			else;
		}
		//std::cout<<"==============Finished with pdf  "<<f<<"=====================\n";
	}
}/*----------------------end of parallel region-----------------------------------------------------------------------------*/

}

void 
Domain::computeInteractionForces(CLbmCase* pCase, Domain& domainVariables)
{
	cgsize_t primary(0), secondary(1), xVel(0), fx(1), fy(2), fz(3), node(0), node1(0), node2(0);
	Node::NodeValueType_t G(0.001), G_pri(-1.0), G_sec(0.0);
	PdfDomain::LatticeType_t::size_type f;
	PdfDomain::LatticeType_t::size_type end_f = GetDomainMaterial()->at(primary)->GetLatticePdf()->size();
	cgsize_t end = pCase->m_grid.m_NVertex;
#pragma omp parallel private(f)
{/*-----------------------start of parallel region-----------------------------------------------------------------------------*/
#pragma omp for private(node1)
	for(node1 = 1; node1 <= end; node1++)
	{
		domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(primary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(primary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(primary)->at(node1-1)->m_nodeVal = 0.0;

		domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(secondary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(secondary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(secondary)->at(node1-1)->m_nodeVal = 0.0;

		//computePsi(pCase, domainVariables, node1, primary);
		//computePsi(pCase, domainVariables, node1, secondary);
		//computePhi(pCase, domainVariables, node1, primary);
		//computePhi(pCase, domainVariables, node1, secondary);
	}
	for(f =1; f < end_f; f++)
	{
#pragma omp for private(node)
		for(node = 1; node <= end; node++)
		{
			cgsize_t phi(7), psi(0);// note phi is 7 
			if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(primary)->at(node-1)->m_isSolid == false)
			{
				//computeSummationIntra(pCase, domainVariables, f, node, primary,   G_pri, psi);
				//computeSummationIntra(pCase, domainVariables, f, node, secondary, G_sec, psi);
				//computeSummationInter(pCase, domainVariables, f, node, primary, secondary, G, phi);		
				//computeSummationInter(pCase, domainVariables, f, node, secondary, primary, G, phi);
				computeSummationInter(pCase, domainVariables, f, node, primary, secondary, G, psi);		
				computeSummationInter(pCase, domainVariables, f, node, secondary, primary, G, psi);
				//computeSummation(pCase, domainVariables, f, node, primary, G);
				//computeSummation(pCase, domainVariables, f, node, secondary, G);
			}
			else;
		}
	}
/*#pragma omp for private(node2)
	for(node2 = 1; node2 <= end; node2++)
	{
		if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(primary)->at(node2-1)->m_isSolid == false)
		{
			computeAllTerms(pCase, domainVariables, node2);
		}
		else;
	}*/
}/*----------------------end of parallel region-----------------------------------------------------------------------------*/
	
}


void 
Domain::computeSurfaceForces(CLbmCase* pCase, Domain& domainVariables)
{
	
	cgsize_t primary(0), secondary(1), xVel(0), fsurf_x(4), fsurf_y(5), fsurf_z(6), node(0), node1(0), node2(0);
	Node::NodeValueType_t Gads_pri, Gads_sec;
	PdfDomain::LatticeType_t::size_type f;
	PdfDomain::LatticeType_t::size_type end_f = GetDomainMaterial()->at(primary)->GetLatticePdf()->size();
	cgsize_t end = pCase->m_grid.m_NVertex;
	Gads_pri =     0.0146;
	Gads_sec =    -0.0146;
#pragma omp parallel private(f)
{/*-----------------------start of parallel region-----------------------------------------------------------------------------*/
#pragma omp for private(node1)
	for(node1 = 1; node1 <= end; node1++)
	{
		domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(primary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(primary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(primary)->at(node1-1)->m_nodeVal = 0.0;

		domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(secondary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(secondary)->at(node1-1)->m_nodeVal = 0.0;
		domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(secondary)->at(node1-1)->m_nodeVal = 0.0;
	}
	
	for(f =1; f < end_f; f++)
	{
#pragma omp for private(node)
		for(node = 1; node <= end; node++)
		{
			if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(primary)->at(node-1)->m_isSolid == false)
			{
				computeSurfaceSummation(pCase, domainVariables, f, node, primary,   Gads_pri);
				computeSurfaceSummation(pCase, domainVariables, f, node, secondary, Gads_sec);
			}
			else;
		}
	}
	
#pragma omp for private(node2)
	for(node2 = 1; node2 <= end; node2++)
	{
		if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(primary)->at(node2-1)->m_isSolid == false)
		{
			computeSurfaceAllTerms(pCase, domainVariables, node2);
		}
		else;
	}
	
}/*----------------------end of parallel region-----------------------------------------------------------------------------*/
	
}

void 
Domain::computePsi(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G)
{
	cgsize_t xVel(0), rho(4), press(8), pseudoP(0);
	Node::NodeValueType_t D3Q15_typeII(1.0/72.0), nodeValue(0.0), pressure(0.0), cs_square(1.0/3.0), c0(10.0), temp1(0.0), temp2(0.0);
	
	if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(sub)->at(node-1)->m_isSolid == false)
		nodeValue = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(sub)->at(node-1)->m_nodeVal;
	else
		nodeValue = wall_VirtualDensity(domainVariables, node, sub);

	
	if(sub == 0)
	{
		//Use Yuan-Schaefer (from Kamali et al)
		pressure = pressure_EOS_CS(pCase, domainVariables, node, sub, nodeValue);
		//set value of pressure
		domainVariables.GetDomainVariables()->at(press)->GetVariable()->at(sub)->at(node-1)->m_nodeVal = pressure;
		temp1 = (2.0 / G) * (pressure - (cs_square * nodeValue));
			
		if(temp1 < 0)
		{
			std::cout<<"node = "<<node<<" => Error! subtracting = "<<(pressure - (cs_square * nodeValue))<<" is positive <="<< " Pressure = "<<pressure<<
			" density = "<<nodeValue<<" cs_square * density = "<<(cs_square * nodeValue)<<" at x-coord = "<<
			domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(sub)->at(node-1)->m_xcoord<<"\n";
			std::cin.get();
		}
		nodeValue = sqrt(temp1);		
	}
	else; 
	domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub)->at(node-1)->m_nodeVal = nodeValue;
}

Node::NodeValueType_t 
Domain::computePhi(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G)
{
	cgsize_t xVel(0), rho(4), primary(0), secondary(1);
	Node::NodeValueType_t nodeValue(0.0), a0(0.005), rho_sec(0.0), rho_pri(0.0); 
	Node::NodeValueType_t rho_sec0 =  0.03;
	Node::NodeValueType_t rho_pri0 = -0.08/log(a0);
	//std::cout<<" rho_pri0 = "<<rho_pri0<<"\n";
	//std::cin.get();
	//Set psi as density
	if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(sub)->at(node-1)->m_isSolid == false)
	{
		if(sub == primary)
		{
			rho_sec = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
			nodeValue = (1.0 - exp(-rho_sec / rho_sec0));	
		}
		else if(sub == secondary)
		{
			rho_pri = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
			nodeValue = (a0 - exp(-rho_pri / rho_pri0));
		}
		else;
	}
	return nodeValue;
}


void 
Domain::computePhi(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub)
{
	cgsize_t xVel(0), rho(4), phi(7), primary(0), secondary(1);
	Node::NodeValueType_t nodeValue(0.0), a0(0.005), rho_sec(0.0), rho_pri(0.0); 
	Node::NodeValueType_t rho_sec0 =  0.003;
	Node::NodeValueType_t rho_pri0 = -0.008/log(a0);

	//Set psi as density
	if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(sub)->at(node-1)->m_isSolid == false)
	{
		if(sub == primary)
		{
			rho_sec = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
			nodeValue = (1.0 - exp(-rho_sec / rho_sec0));	
		}
		else
		{
			rho_pri = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
			nodeValue = (a0 - exp(-rho_pri / rho_pri0));
		}
	}
	domainVariables.GetDomainTempVariables()->at(phi)->GetVariable()->at(sub)->at(node-1)->m_nodeVal = nodeValue;
}


Node::NodeValueType_t
Domain::wall_VirtualDensity(Domain& domainVariables, cgsize_t node, cgsize_t sub)
{
	cgsize_t xVel(0), rho(4), solidNeighbor(0), num_solidNeighbor(0);
	Node::NodeValueType_t virtualDensity(0.0), DeltaS_pri(0.0), DeltaS_sec(0.0);
	PdfDomain::LatticeType_t::size_type f;
	PdfDomain::LatticeType_t::size_type end_f = GetDomainMaterial()->at(sub)->GetLatticePdf()->size();
	
	for(f =1; f < end_f; f++)
	{
		solidNeighbor = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;

		if(solidNeighbor != 0 && domainVariables.GetDomainTempVariables()->at(xVel)->GetVariable()->at(sub)->at(solidNeighbor-1)->m_isSolid == false)
		{
			num_solidNeighbor += 1;
			virtualDensity += domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(sub)->at(solidNeighbor-1)->m_nodeVal;
		}
		else;
	}
	virtualDensity /= num_solidNeighbor;
	
	//fabien Jansen or Yu Chen model
	//check to see if there is a wettability patch (implement different deltas for patch)
	if (domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(sub)->at(node-1)->m_wettingFlag == true)
	{
		//change to a method later
		/**********************************************************************/
		if(sub == 0)
		{
			virtualDensity +=  0.0 /*DeltaS_sec*/;
		}
		else
		{
			virtualDensity +=  0.0 /*DeltaS_pri*/;
		}
		/**********************************************************************/
	}
	else
	{
		
		/**********************************************************************/
		if(sub == 0)
		{
			virtualDensity += DeltaS_pri;
		}
		else
		{
			virtualDensity += DeltaS_sec;
		}
		/**********************************************************************/
	}
	
	//update value of density on the wall with the virtual density
	//domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(sub)->at(node-1)->m_nodeVal = 0.0;
	
	return virtualDensity;
}

Node::NodeValueType_t
Domain::pressure_EOS_PR(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t density)
{
	//Peng and Robinson EOS
	Node::NodeValueType_t alpha(0.0), pressure(0.0), accentric(0.344), Temp(0.0), Temp_crit(0.0), a(2.0/49.0), b(2.0/21.0), R(1.0), temp1(0.0);
	
		
	Temp_crit = (0.1702 * a) / (R * b);
	Temp      =  0.75 * (Temp_crit); // specify temperature from coexistence curve here 
	temp1     = (1.0 + (0.37464 + 1.54226 * accentric - 0.26992 * accentric * accentric) * (1.0 - sqrt(Temp/Temp_crit))); 
	alpha     = temp1 * temp1;

	pressure = (density * R * Temp) / (1.0 - b*density)  - 
		   (a * alpha * density * density) / (1.0 + (2.0 * b * density) - (b*b*density*density));
	return pressure;
}

Node::NodeValueType_t
Domain::pressure_EOS_CS(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t density)
{
	//Carnahan and Starling EOS
	Node::NodeValueType_t pressure(0.0), rho_crit(0.0), press_crit(0.0), Temp(0.0), Temp_crit(0.0), a(1.0/*0.09926*/), b(4.0/*0.18727*/), R(1.0/*0.2*/), temp1(0.0), temp2(0.0);
	
	Temp_crit  = (0.3773 * a) / (R * b);
	Temp       =  0.8 * (Temp_crit); // specify temperature from coexistence curve here 
	temp1      = (b * density) / 4.0; 
	temp2      = (1.0 - temp1);

	pressure = (density * R * Temp) * ((1.0 + temp1 + (temp1*temp1) - (temp1*temp1*temp1))/(temp2*temp2*temp2)) - (a * density * density);  
	return pressure;
}

Node::NodeValueType_t
Domain::pressure_EOS_vdW(CLbmCase* pCase, Domain& domainVariables, cgsize_t node, cgsize_t sub, Node::NodeValueType_t density)
{
	//van der Waal
	Node::NodeValueType_t pressure(0.0), Temp(0.0), a(9.0/49.0), b(2.0/21.0), R(1.0), Temp_crit(4.0/7.0);
	Temp = 0.58 * Temp_crit;
	pressure = (density * R * Temp)/(1.0 - b*density) - (a*density*density);  
	return pressure;
}

void 
Domain::computeSummationInter(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub1, cgsize_t sub2, Node::NodeValueType_t G, cgsize_t pseudoP)
{
	cgsize_t xvel(0);
	cgsize_t neighborNode = GetInteractionNeighbor(pCase,f, node, sub1);
	if(neighborNode != 0)
	{
		//if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(sub1)->at(neighborNode-1)->m_isSolid == false)// if I am using Marty and Chen
		//{
			if (pCase->m_model.m_solver == CLbmModel::D3Q15)
			{
				D3Q15_FCohesionAdhesion(pCase, domainVariables, f, node, neighborNode, sub1, sub2, G, pseudoP);
			}
			else if (pCase->m_model.m_solver == CLbmModel::D3Q19)
			{
				D3Q19_FCohesionAdhesion(pCase, domainVariables, f, node, neighborNode, sub1, sub2, G, pseudoP);
			}
			else;
		//}
		//else;
	}
	else;
}

void 
Domain::computeSummationIntra(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G, cgsize_t pseudoP)
{
	cgsize_t xVel(0), rho(4);
	Node::NodeValueType_t Gads(0.0);
	cgsize_t neighborNode = GetInteractionNeighbor(pCase,f, node, sub);

	if(neighborNode != 0)
	{
		if (pCase->m_model.m_solver == CLbmModel::D3Q15)
		{
			if(domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(sub)->at(neighborNode-1)->m_isSolid == false)
			{
				D3Q15_FCohesionAdhesion(pCase, domainVariables, f, node, neighborNode, sub, sub, G, pseudoP);
			}
			else
			{
				if (domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(sub)->at(neighborNode-1)->m_wettingFlag == true)
				{
					Gads = -1.20; //make hydrophilic 
				}
				else
				{
					Gads = -0.84; //make hydrophobic
				}
				D3Q15_FCohesionAdhesion(pCase, domainVariables, f, node, neighborNode, sub, sub, Gads, pseudoP);
			}		
		}
		else if (pCase->m_model.m_solver == CLbmModel::D3Q19)
		{
			//space holder
		}
		else;
	}
	else;
}

void 
Domain::computeSummation(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub, Node::NodeValueType_t G)
{
	cgsize_t psi(0), fx(1), fy(2), fz(3), xvel(0), tempNode(0), tempf(0); 
	Node::NodeValueType_t D3Q15_typeII(1.0/72.0);
	
	if	(GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk() == D3Q15_typeII)
	{
		//G = G / sqrt(3.0);
		G = G/24.0;
	}
	else
	{
		G = G / 3.0;
	}

	cgsize_t neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
	if(!neighborNode)
	{
		tempNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum;

		if(tempNode == 0) //boundary is not periodic (open boundary)
		{
			tempf = switchPdf(pCase, f);
			neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(tempf)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
			//change to a method or subroutine later
			if(!neighborNode)
			{
				//f1_neighbor  = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
				//f2_neighbor  = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
				//neighborNode = (f1_neighbor != 0) ? f1_neighbor : f2_neighbor;
				neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
			}
		}
		else //boundary is periodic
		{
			neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(tempNode-1)->m_neighborNum;
		}
	}
	else;
	
	if(neighborNode != 0)
	{
		if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(sub)->at(neighborNode-1)->m_isSolid == false)
		{
			if (pCase->m_model.m_solver == CLbmModel::D3Q15)
			{
				domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
					(-G * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCx()
					* domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(sub)->at(neighborNode-1)->m_nodeVal);
						
				domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
					(-G * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCy()
					* domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(sub)->at(neighborNode-1)->m_nodeVal);
						
				domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
					(-G * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCz()
					* domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(sub)->at(neighborNode-1)->m_nodeVal);			
			}
			else if (pCase->m_model.m_solver == CLbmModel::D3Q19)
			{
				domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
					(-G * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCx()
					* GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk()
					* domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(sub)->at(neighborNode-1)->m_nodeVal);
						
				domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
					(-G * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCy()
					* GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk()
					* domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(sub)->at(neighborNode-1)->m_nodeVal);
						
				domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
					(-G * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCz()
					* GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk()
					* domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(sub)->at(neighborNode-1)->m_nodeVal);
							
			}
			else;
		}
		else;
	}
	else;				
}


cgsize_t 
Domain::GetInteractionNeighbor(CLbmCase* pCase, cgsize_t f, cgsize_t node, cgsize_t sub)
{
	cgsize_t tempNode(0), tempf(0), neighborNode(0);

	neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
	
	if(!neighborNode)
	{
		tempNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum;
		
		if(tempNode == 0) //boundary is not periodic (open boundary)
		{
			tempf = switchPdf(pCase, f);
			neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(tempf)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
			//change to a method or subroutine later
			if(!neighborNode)
			{
				//f1_neighbor  = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
				//f2_neighbor  = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;
				//neighborNode = (f1_neighbor != 0) ? f1_neighbor : f2_neighbor;
				neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
			}
		}
		else //boundary is periodic
		{
			neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(tempNode-1)->m_neighborNum;
		}
	}
	else;
	return neighborNode;
}	

void 
Domain::D3Q15_FCohesionAdhesion(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t neighborNode, cgsize_t sub1, cgsize_t sub2, Node::NodeValueType_t G, cgsize_t pseudoP)
{
	cgsize_t fx(1), fy(2), fz(3), psi(0), phi(7); 
	Node::NodeValueType_t D3Q15_typeII(1.0/72.0), pseudoP_node(0.0), pseudoP_neighborNode(0.0);
	
	//compute psi and/or phi first
	if(pseudoP == psi)
	{
		pseudoP_node           = domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal;
		pseudoP_neighborNode   = domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub2)->at(neighborNode-1)->m_nodeVal;	
	}
	else if(pseudoP == phi)
	{
		pseudoP_node           = computePhi(pCase, domainVariables, node, sub1, G);
		pseudoP_neighborNode   = computePhi(pCase, domainVariables, neighborNode, sub2, G);
	}
	else;
	if(GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetWk() == D3Q15_typeII)
	{
		//G = G / sqrt(3.0) ;
		G = G / 24.0;	 
	}
	else
	{
		//G =  G;
		G = G / 3.0;
	}

	domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal += 
		(-G * GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetCx()
		* pseudoP_node * pseudoP_neighborNode);
						
	domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal += 
		(-G * GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetCy()
		* pseudoP_node * pseudoP_neighborNode);
						
	domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal += 
		(-G * GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetCz()
		* pseudoP_node * pseudoP_neighborNode);
}

void 
Domain::D3Q19_FCohesionAdhesion(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t neighborNode, cgsize_t sub1, cgsize_t sub2, Node::NodeValueType_t G, cgsize_t pseudoP)
{
	cgsize_t fx(1), fy(2), fz(3); 
	Node::NodeValueType_t D3Q19_typeII(1.0/72.0);

	if(GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetWk() == D3Q19_typeII)
	{
		G = 0;
	}
	else;
	
	domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal += 
		(-G * GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetCx()
		* GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetWk()
		* domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal
		* domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub2)->at(neighborNode-1)->m_nodeVal);
						
	domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal += 
		(-G * GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetCy()
		* GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetWk()
		* domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal
		* domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub2)->at(neighborNode-1)->m_nodeVal);
						
	domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal += 
		(-G * GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetCz()
		* GetDomainMaterial()->at(sub1)->GetLatticePdf()->at(f)->pdf()->GetWk()
		* domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub1)->at(node-1)->m_nodeVal
		* domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(sub2)->at(neighborNode-1)->m_nodeVal);
}

cgsize_t 
Domain::switchPdf(CLbmCase* pCase, cgsize_t f)
{
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), 
					f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14),
					f15(15), f16(16), f17(17), f18(18);
	
	if(f == f1)      return f2;
	else if(f == f2) return f1;
	else if(f == f3) return f4;
	else if(f == f4) return f3;
	else if(f == f5) return f6;
	else if(f == f6) return f5;
	else;
	if(pCase->m_model.m_solver == CLbmModel::D3Q15)
	{
		if(f == f7)       return f14;
		else if(f == f8)  return f13;
		else if(f == f9)	return f12;
		else if(f == f10)	return f11;
		else if(f == f11)	return f10;
		else if(f == f12)	return f9;
		else if(f == f13)	return f8;
		else if(f == f14)	return f7;
		else;
	}
	else if (pCase->m_model.m_solver == CLbmModel::D3Q19)
	{ //Not correct: fix later
		if(f == f7)       return f10;
		else if(f == f8)  return f9;
		else if(f == f9)	return f8;
		else if(f == f10)	return f7;
		else if(f == f11)	return f14;
		else if(f == f12)	return f13;
		else if(f == f13)	return f12;
		else if(f == f14)	return f11;
		else if(f == f15)	return f18;
		else if(f == f16)	return f17;
		else if(f == f17)	return f16;
		else if(f == f18)	return f15;
		else
			return 0;
	}
}
void 
Domain::computeSurfaceSummation(CLbmCase* pCase, Domain& domainVariables, cgsize_t f, cgsize_t node, cgsize_t sub, Node::NodeValueType_t Gads)
{
	cgsize_t fsurf_x(4), fsurf_y(5), fsurf_z(6), xvel(0), tempNode(0);
	Node::NodeValueType_t D3Q15_typeII(1.0/72.0);
	
	if(GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk() == D3Q15_typeII)
	{
		Gads = Gads / sqrt(3.0);
		//Gads = 0.0;
	}
	else;
	
	cgsize_t neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_neighborNum;

	if(neighborNode == 0)
	{
		tempNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum;

		if(tempNode != 0)
		{
			neighborNode = GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(tempNode-1)->m_neighborNum;
		}
		else;
	}
	else;

	if(neighborNode != 0)
	{
		if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(sub)->at(neighborNode-1)->m_isSolid == true)
		{
			if (pCase->m_model.m_solver == CLbmModel::D3Q15)
			{
				//using pan et al
				domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
							(-Gads * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCx());

				domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
							(-Gads * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCy());

				domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
							(-Gads * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCz());	

			}
			else if(pCase->m_model.m_solver == CLbmModel::D3Q19)
			{
				domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
								(-Gads * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCx() * 
								GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk());

				domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
								(-Gads * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCy() * 
								GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk());

				domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(sub)->at(node-1)->m_nodeVal += 
								(-Gads * GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetCz() * 
								GetDomainMaterial()->at(sub)->GetLatticePdf()->at(f)->pdf()->GetWk());
			}
			else;
		}
		else;
	}
	else;
}

void 
Domain::computeAllTerms(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	Node::NodeValueType_t  x_temp(0.0), y_temp(0.0), z_temp(0.0);
	cgsize_t psi(0), fx(1), fy(2), fz(3), xvel(0);
	cgsize_t primary(0), secondary(1);	
		
	if (domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(secondary)->at(node-1)->m_isSolid == false)
	{		
		x_temp = domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
		//obtain fx for primary and secondary materials
		domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal * 
			domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
				
		domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(primary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(primary)->at(node-1)->m_nodeVal * x_temp;
       
		y_temp = domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
		//obtain fy for primary and secondary materials
		domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal *
			domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

		domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(primary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(primary)->at(node-1)->m_nodeVal * y_temp;
					
		z_temp = domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
		//obtain fx for primary and secondary materials
		domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal *
			domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

		domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(primary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(primary)->at(node-1)->m_nodeVal * z_temp;		
	}
	else;
}

void 
Domain::computeSurfaceAllTerms(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	cgsize_t fsurf_x(4), fsurf_y(5), fsurf_z(6), psi(0), xvel(0);
	cgsize_t primary(0), secondary(1);
	
	if (domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_isSolid == false)
	{	
		domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal *= 
					(domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(primary)->at(node-1)->m_nodeVal);
		domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal *= 
					(domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(primary)->at(node-1)->m_nodeVal);
		domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal *= 
					(domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(primary)->at(node-1)->m_nodeVal);

		domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal *= 
					(domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal);
		domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal *= 
					(domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal);
		domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal *= 
					(domainVariables.GetDomainTempVariables()->at(psi)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal);

	}
	else;
}
//=============================================================================================//


Domain_D2Q9::Domain_D2Q9(CLbmCase* pCase)
{
	//Initialize materials in domain
	//if material is incompressible
	if( pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		cgsize_t twoPhase(2);
		for(cgsize_t i=0 ; i < twoPhase; i++)
		{
			m_domainMaterial.push_back(new D2Q9Incomp_domain(pCase, i));
		}

		//Initialize domain variable	//Initialize domain variables
		m_domainVariables.push_back(new CXVelocity(pCase));
		m_domainVariables.push_back(new CYVelocity(pCase));
		m_domainVariables.push_back(new CZVelocity(pCase));
		m_domainVariables.push_back(new CVelocityMag(pCase));
		m_domainVariables.push_back(new CDensity(pCase));
	}
	else
	{

	}
}
Domain_D2Q9::Domain_D2Q9(CLbmCase* pCase, cgsize_t nvertex)
{
	//Initialize materials in domain
	//if material is incompressible
	m_domainMaterial.push_back(new D2Q9Incomp_domain(pCase, nvertex));


	//Initialize domain variables
	m_domainVariables.push_back(new CXVelocity(pCase, nvertex));
	m_domainVariables.push_back(new CYVelocity(pCase, nvertex));
	m_domainVariables.push_back(new CZVelocity(pCase, nvertex));
	m_domainVariables.push_back(new CVelocityMag(pCase, nvertex));
	m_domainVariables.push_back(new CDensity(pCase, nvertex));
}

void 
Domain_D2Q9::computeCompositeVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	//place holder
}

void 
Domain_D2Q9::computeOverallVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	//place holder
}
void 
Domain_D2Q9::computeOverallVarsSgl(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	//place holder
}
void 
Domain_D2Q9::computeIndividualVars(CLbmCase* pCase, Domain& domainVariables, MaterialType_t::size_type mat_index, cgsize_t node)
{

	//place holder
}
void 
Domain_D2Q9::convergence(Domain& domainVariables, CLbmCase* pCase)
{
	cgsize_t xvel(0), yvel(1);
	m_convergenceData.x_residual = domainVariables.GetDomainVariables()->at(xvel)->residual();
	m_convergenceData.y_residual = domainVariables.GetDomainVariables()->at(yvel)->residual();
	m_convergenceData.z_residual = 0.0;
}

Convergence*   
Domain_D2Q9::GetConvergence() 
{
	return &m_convergenceData;
}

Domain_D2Q9::~Domain_D2Q9()
{
	
}

Domain::MaterialType_t*
Domain_D2Q9::DomainMaterial()
{
	return &m_domainMaterial;
}

Domain::VariableType_t*
Domain_D2Q9::DomainVariables()
{
	return &m_domainVariables;
}

Domain::VariableType_t*
Domain_D2Q9::DomainTempVariables()
{
	return &m_domainTempVariables;
}

Domain::PatchNodeType_t*
Domain_D2Q9::DomainSolidObjectNodes()
{
	return &m_domainSolidObjectNodes;
}

Domain::MicroPillarsType_t*
Domain_D2Q9::DomainMicroPillars()
{
	return &m_domainMicroPillars;
}

void 
Domain_D2Q9::SetNumThreads(cgsize_t  num)
{
	//place holder
}
cgsize_t 
Domain_D2Q9::GetNumThreads() const
{
	//place holder
	return 0;
}
//LEVEL TWO:invoked from Domain
//==============================================================================================//

//function to obtain the density at the wall given the velocity at the wall
Node::NodeValueType_t 
PdfDomain::density(Domain& domainVariables, cgsize_t dNqN,cgsize_t node, Node::NodeValueType_t c, cgsize_t vel_index,
					cgsize_t fI,  cgsize_t fII,  cgsize_t fIII,	cgsize_t ffI, cgsize_t ffII, cgsize_t ffIII,
					cgsize_t fIV, cgsize_t fV, cgsize_t ffIV, cgsize_t ffV)
{
	Node::NodeValueType_t wallVel(0.0), temp1(0.0), temp2(0.0), temp3(0.0);
	wallVel = domainVariables.GetDomainVariables()->at(vel_index)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
	temp1 = 1.0 / (1.0 +  c * wallVel );
	temp2 =		GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal   + GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal 
					+ GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal;
	if(fIV != 0 && fV != 0)
	{
		temp2 = temp2 +  GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal  + GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal;
	}

	temp3 = 2*(GetLatticePdf()->at(ffI)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal  +  GetLatticePdf()->at(ffII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal 
					+ GetLatticePdf()->at(ffIII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal);
 if(ffIV != 0 && ffV != 0)
	{
		temp3 = temp3 + 2.0 * ( GetLatticePdf()->at(ffIV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal  + GetLatticePdf()->at(ffV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal);
 }
	return temp1 * (temp2 + temp3);
}


//function to obtain the velocity at the wall given the density
Node::NodeValueType_t 
PdfDomain::velocity(Domain& domainVariables, cgsize_t dNqN, cgsize_t node, Node::NodeValueType_t c1, Node::NodeValueType_t c2,
					cgsize_t fI,  cgsize_t fII,  cgsize_t fIII,	cgsize_t ffI, cgsize_t ffII, cgsize_t ffIII,
					cgsize_t fIV, cgsize_t fV, cgsize_t ffIV, cgsize_t ffV)
{
	Node::NodeValueType_t wallRho(0.0), temp1(0.0), temp2(0.0), temp3(0.0);
	cgsize_t rho(4);
	wallRho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;

	temp2 =		GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal   + GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal 
					+ GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal;
	if(fIV != 0 && fV != 0)
	{
		temp2 = temp2 +  GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal  + GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal;
	}

	temp3 = 2*(GetLatticePdf()->at(ffI)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal  +  GetLatticePdf()->at(ffII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal 
					+ GetLatticePdf()->at(ffIII)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal);
 if(ffIV != 0 && ffV != 0)
	{
		temp3 = temp3 + 2.0 * ( GetLatticePdf()->at(ffIV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal  + GetLatticePdf()->at(ffV)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal);
 }

	temp1 = c2  + c1 * ((temp2 + temp3)/wallRho);
	return temp1;
}

PdfDomain::LatticeType_t* 
PdfDomain::GetLatticePdf()
{
	return LatticePdf();
}
cgsize_t 
PdfDomain::GetMaterialIndex()
{
	return MaterialIndex();
}
void
PdfDomain::initialize(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index)
{
	
	for(LatticeType_t::size_type i =0; i < GetLatticePdf()->size(); i++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"material = "<< mat_index <<" and pdf = "<<i<<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetLatticePdf()->at(i)->lbmInitialize(pCase, domainVariables,mat_index);
	}
}
void
PdfDomain::collide(CLbmCase* pCase, Domain& domainVariables, cgsize_t index)
{	
	//obtain relaxation type
	CLbmModel::RELAXTIME input = pCase->m_model.m_relaxTime;
	//obtain relaxation
	Node::NodeValueType_t tau(0.6);
	std::map <cgsize_t,CLbmElement::MaterialRecord_t>::const_iterator iter = pCase->m_element.m_sectionsMaterial.find(index);
	//tau     = 3.0 * iter->second.second.k_viscosity + 0.5;
	
	PdfDomain::LatticeType_t::size_type j;
	PdfDomain::LatticeType_t::size_type end = GetLatticePdf()->size();
	//cgsize_t end=2;
	//cgsize_t i;
	//cgsize_t TID;
#pragma omp parallel private(j)
{/*-----------start of parallel region----------------------*/
//#pragma omp for private(j) 
        for(j=0; j < end ; j++)
	{	
		//TID = omp_get_thread_num();
		//printf("Thread: %d is executing loop work for iteration= %d \n", TID, i);
		GetLatticePdf()->at(j)->lbmCollide(this, domainVariables,input, tau, index, j);
		//Now transform from moment space to velocity space
	}/*-----------end of omp for---------------------*/
}/*---------end of parallel region-------------------------*/
}
void
PdfDomain::Lbm_write(std::fstream& outputFile)
{
	for(LatticeType_t::size_type i =0; i < GetLatticePdf()->size(); i++)
		GetLatticePdf()->at(i)->lbmWrite(outputFile);
}

void
PdfDomain::Lbm_read(std::fstream& inputFile)
{
	for(LatticeType_t::size_type i =0; i < GetLatticePdf()->size(); i++)
		GetLatticePdf()->at(i)->lbmRead(inputFile);
}
void
PdfDomain::post_collide(CLbmCase* pCase, Domain& domainVariables, cgsize_t index)
{
	BCType_t boundType;
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), interior(2);
	clock_t time;
	time = clock();
	//update wall PDF after collision
	for(cgsize_t j =1; j <= pCase->m_boco.m_nBocos; j++)
	{
		/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" about to post collide for boundary = "<< j <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (j-1)) != interior)
			{
				boundType = *(pCase->m_boco.m_bocoType + (j-1));
				pPnts     = *(pCase->m_boco.m_pnts     + (j-1));
				numPnts   = *(pCase->m_boco.m_npnts    + (j-1));
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//WALL-update distribution functions on the wall
				
				if(boundType == BCWall) //If the boundary is not a wall we will first stream then obtain unknown pdfs
					{
						/*/////////////////////////////////////////////////////////////////////////////////////////////////////////////
						std::cout<<" boundary = "<< j <<" is a wall, therefore eligible\n";
						std::cin.get();
						////////////////////////////////////////////////////////////////////////////////////////////////////////////*/						
						bound(pCase,domainVariables,index, j, boundType, pPnts,numPnts);
					}
				
				else if(boundType == BCDirichlet)
				{
					if(*(pCase->m_boco.m_bocoExtrapolationBC + (j-1)))
						{
							bound(pCase,domainVariables, index, j, boundType, pPnts,numPnts);
						}
				}
		}
	}
	time = (clock() - time);/*CLOCKS_PER_SEC;*/

	if(!domainVariables.GetDomainSolidObjectNodes()->empty()) 
	{
		LbmDomain::FORMULA f_type  = LbmDomain::HALF_BB;
		LbmDomain::BOUND   b_type  = LbmDomain::BOUNCEBACK;

		for(LatticeType_t::size_type f =1; f < GetLatticePdf()->size(); f++)
		{
			//call nodeBounder to set parameter and update each pdf
			GetLatticePdf()->at(f)->solidObjectNodeBounder(domainVariables,index,domainVariables.GetDomainSolidObjectNodes(),b_type,f_type);
		}
	}


	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout<<"time to perform post collision = "<< time <<"\n";
	//std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void
PdfDomain::stream(CLbmCase* pCase, Domain& domainVariables)
{
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), interior(2);
	clock_t time;
	time = clock();
	for(cgsize_t j =1; j <= pCase->m_boco.m_nBocos; j++)
	{
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (j-1)) != interior)
			{
				pPnts     = *(pCase->m_boco.m_pnts     + (j-1));
				numPnts   = *(pCase->m_boco.m_npnts    + (j-1));
        			stream_utility(domainVariables, pPnts , numPnts);
			}
	}
	time = (clock() - time);/*CLOCKS_PER_SEC;*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout<<"time to perform streaming = "<< time <<"\n";
	//std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void
PdfDomain::post_stream(CLbmCase* pCase, Domain& domainVariables, cgsize_t index)
{
	BCType_t boundType;
	cgsize_t* pPnts = NULL;
	cgsize_t  numPnts(0), interior(2);
	clock_t time;
	time = clock();
	for(cgsize_t j =1; j <= pCase->m_boco.m_nBocos; j++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" about to post stream for boundary = "<< j <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//do this if boundary condition is not interior
		if(*(pCase->m_boco.m_bocoUserDef + (j-1)) != interior)
			{
				boundType = *(pCase->m_boco.m_bocoType + (j-1));
				pPnts     = *(pCase->m_boco.m_pnts     + (j-1));
				numPnts   = *(pCase->m_boco.m_npnts    + (j-1));
        //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				//Non-WALL-update distribution functions on the non-wall boundaries
				
				if(boundType != BCWall) //If the boundary is not a wall we will first stream then obtain unknown pdfs
				{
					//if(  (*(pCase->m_boco.m_bocoPressureBC + (j-1)))
						//|| (*(pCase->m_boco.m_bocoVelocityBC + (j-1))) )
						//{
							bound(pCase,domainVariables, index, j, boundType, pPnts,numPnts);
						//}
				}
				else;
		}
	}
	time = (clock() - time);/*CLOCKS_PER_SEC;*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
		//std::cout<<"time to perform post streaming = "<< time <<"\n";
		//std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void
PdfDomain::stream_utility(Domain& domainVariables, cgsize_t* pPnts , cgsize_t  numPnts)
{
	/*Will only stream if node is:
	   - a stream start node
		 - nodeStreamed flag is false (i.e distribution function has not been streamed by another boundary)
	*/
	/*cgsize_t i;
	LatticeType_t::size_type k;
	
	for(i =0; i < numPnts; i++)
		{
			cgsize_t  node(0);
			node = *(pPnts + i);
			//stream all Pdfs on node as we iterate through each node on boundary
			for(k =1; k < GetLatticePdf()->size(); k++)
				{
					//check to see if node streams for the current distribution function
					if(GetLatticePdf()->at(k)->pdf()->GetPdfunction()->at(node-1)->m_streamStart)
						{
							//check to see if node has not been streamed by other boundary
							if(GetLatticePdf()->at(k)->pdf()->GetPdfunction()->at(node-1)->m_nodeStreamed == false)
								GetLatticePdf()->at(k)->lbmStream(node);
							else;
					}
					else;
			}
	}*/

	//////////////modification for parallel/////////////////////////
	//double starttime;
	//double endtime;
	//starttime = omp_get_wtime();
	cgsize_t i; 
	LatticeType_t::size_type k;
	LatticeType_t::size_type end = GetLatticePdf()->size();
	cgsize_t  node(0);
	//cgsize_t end1=3;
	//cgsize_t end2=13;
	//cgsize_t k;
	//cgsize_t i;
	//cgsize_t inner_TID;
#pragma omp parallel private(k)
{/*------------begin parallel region----------------*/
#pragma omp for private(i, node)
	for(k =1; k < end/*GetLatticePdf()->size()*/; k++)
	{/*-------------begin omp for loop---------------*/
		for(i =0; i < numPnts; i++)
		{/*---------begin omp for loop-------------------*/
			node = *(pPnts + i);
			//stream all Pdfs on node as we iterate through each node on boundary
			
					//check to see if node streams for the current distribution function
					if(GetLatticePdf()->at(k)->pdf()->GetPdfunction()->at(node-1)->m_streamStart)
						{
							//check to see if node has not been streamed by other boundary
							if(GetLatticePdf()->at(k)->pdf()->GetPdfunction()->at(node-1)->m_nodeStreamed == false)
								GetLatticePdf()->at(k)->lbmStream(node);
							else;
						}
					else;
		}/*----------end of omp for loop-----------------*/
	}/*---------end omp for loop------------------*/
}/*---------------end parallel region------------------*/
}

//apply boundary conditions for each pdf in the lattice
void
PdfDomain::bound(CLbmCase* pCase, Domain& domainVariables, cgsize_t index, cgsize_t bdry, BCType_t boundType, cgsize_t* pPnts , cgsize_t  numPnts)
{
	//Distribution function update flag needs to be false for an update to be performed
	switch (boundType)
		{
			case BCDirichlet:
				{
					bound_dirichlet(pCase, domainVariables, index, pPnts, numPnts, bdry, boundType);
				}
			break;
			case BCOutflow:
				{
					// OUTFLOW BOUNDARY CONDITION
					bound_open_bc(domainVariables, index,  pPnts, numPnts);
				}
			break;
			case BCWall: 
				{
					bound_wall_bc(domainVariables, pCase, index, pPnts, numPnts,boundType, bdry);
				}
			break;
			case BCGeneral: 
				{
					bound_general_bc(domainVariables, pCase, index, pPnts, numPnts, bdry);
				}
			break;	
			default:
			break;
		}
}

//private methods that are utility methods for bound
void
PdfDomain::bound_dirichlet(CLbmCase* pCase, Domain& domainVariables, cgsize_t index, cgsize_t* pPnts, cgsize_t  numPnts, cgsize_t i,  BCType_t boundType)
{
	LbmDomain::BOUND   b_type = LbmDomain::PRESSURE;
	LbmDomain::FORMULA f_type = LbmDomain::ZHOU_HE;
	if(*(pCase->m_boco.m_bocoPressureBC + (i-1)))
		{
			b_type  = LbmDomain::PRESSURE;	
			//specify which formulation to use for velocity inlet
		  f_type = LbmDomain::ZHOU_HE;
			//LbmDomain::FORMULA f_type = LbmDomain::ZHOU_HE_I;
			/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" Pressure boundary for boundary =  "<< i <<"\n";
			std::cin.get();
			////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		}
	else if(*(pCase->m_boco.m_bocoVelocityBC + (i-1)))
	{
		b_type  = LbmDomain::VELOCITY;
		//specify which formulation to use for velocity inlet
		f_type = LbmDomain::ZHOU_HE;
		//LbmDomain::FORMULA f_type = LbmDomain::ZHOU_HE_I;
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" Velocity boundary for boundary =  "<< i <<"\n";
			std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(*(pCase->m_boco.m_bocoExtrapolationBC + (i-1)))
	{
		b_type  = LbmDomain::EXTRAPOLATION;
		f_type = LbmDomain::EXTRAPOLATE;
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<" Extrapolation boundary for boundary =  "<< i <<"\n";
			std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else;
	//continue with dirichlet boundary condition
	for(LatticeType_t::size_type j =1; j < GetLatticePdf()->size(); j++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" bounding dirichlet boundary of  pdf =  "<< j <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetLatticePdf()->at(j)->nodeBounder(domainVariables,index,pPnts,numPnts,boundType,b_type,f_type);
	}
}

void
PdfDomain::bound_open_bc(Domain& domainVariables, cgsize_t index, cgsize_t* pPnts, cgsize_t  numPnts)
{
	for(LatticeType_t::size_type j =1; j < GetLatticePdf()->size(); j++)
	{
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" bounding open boundary of  pdf =  "<< j <<"\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetLatticePdf()->at(j)->openBCBounder(domainVariables,index, pPnts,numPnts);
	}
}

void
PdfDomain::bound_wall_bc(Domain& domainVariables, CLbmCase* pCase, cgsize_t index, cgsize_t* pPnts, cgsize_t  numPnts, BCType_t boundType, cgsize_t i)
{
	cgsize_t Zhou_He_FullBB(0), Zhou_He_HalfBB(1), moving(2), extrapolation(3);
	LbmDomain::FORMULA f_type  = LbmDomain::HALF_BB;
	LbmDomain::BOUND   b_type  = LbmDomain::BOUNCEBACK;
	//specify which formulation to use for bounce back
	if(*(pCase->m_boco.m_bocoWall + (i-1)) == Zhou_He_FullBB)
	{
		f_type = LbmDomain::FULL_BB;
	}
	else if(*(pCase->m_boco.m_bocoWall + (i-1)) == Zhou_He_HalfBB)
	{
		f_type = LbmDomain::HALF_BB;
	}
	else if(*(pCase->m_boco.m_bocoWall + (i-1)) == moving)
	{
		f_type = LbmDomain::MOVING;
	}
	else if(*(pCase->m_boco.m_bocoWall + (i-1)) == extrapolation)
	{
		f_type = LbmDomain::EXTRAPOLATE;
	}
        else;

	//perform bounce back
	for(LatticeType_t::size_type j =1; j < GetLatticePdf()->size(); j++)
	{
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"About to update nodes on wall for material = "<<index<<" boundary " << i <<" pdf =  "<<j<<"\n";
		std::cin.get();
		///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		GetLatticePdf()->at(j)->nodeBounder(domainVariables,index, pPnts,numPnts,boundType,b_type,f_type);
	}
}

void
PdfDomain::bound_general_bc(Domain& domainVariables, CLbmCase* pCase, cgsize_t index, cgsize_t* pPnts, cgsize_t  numPnts, cgsize_t i)
{
	cgsize_t periodic(0), open(1);
	if(*(pCase->m_boco.m_bocoUserDef + (i-1)) == periodic)
	{
		//apply periodic condition
		for(LatticeType_t::size_type j =1; j < GetLatticePdf()->size(); j++)
		{
			/*////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<"About to update nodes  for material = "<<index<<" boundary " << i <<" pdf =  "<<j<<"\n";
			std::cin.get();
			///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			GetLatticePdf()->at(j)->periodicBCBounder(domainVariables,index, pPnts,numPnts);
		}
	}
	else if(*(pCase->m_boco.m_bocoUserDef + (i-1)) == open)
	{
		for(LatticeType_t::size_type j =1; j < GetLatticePdf()->size(); j++)
		{
			GetLatticePdf()->at(j)->openBCBounder(domainVariables,index, pPnts,numPnts);
		}
	}
}

//Derived class
//==============//
D2Q9Incomp_domain::D2Q9Incomp_domain(CLbmCase* pCase, cgsize_t matIndex)
{
	m_materialIndex = matIndex;
	m_latticeType.push_back(new D2Q900_incompDomain(pCase));
	m_latticeType.push_back(new D2Q901_incompDomain(pCase));
	m_latticeType.push_back(new D2Q902_incompDomain(pCase));
	m_latticeType.push_back(new D2Q903_incompDomain(pCase));
	m_latticeType.push_back(new D2Q904_incompDomain(pCase));
	m_latticeType.push_back(new D2Q905_incompDomain(pCase));
	m_latticeType.push_back(new D2Q906_incompDomain(pCase));
	m_latticeType.push_back(new D2Q907_incompDomain(pCase));
	m_latticeType.push_back(new D2Q908_incompDomain(pCase));
	//map name of boundary nodes
	//for(LatticeType_t::size_type i =0; i < m_latticeType.size(); i++)
		//m_latticeType.at(i)->pdf()->setPdfBound(pCase, this);

	//map start of streaming
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8);
	m_latticeType.at(f1)->pdf()->setStreamBound(pCase, m_latticeType.at(f3)->pdf());
	m_latticeType.at(f2)->pdf()->setStreamBound(pCase, m_latticeType.at(f4)->pdf());
	m_latticeType.at(f3)->pdf()->setStreamBound(pCase, m_latticeType.at(f1)->pdf());
	m_latticeType.at(f4)->pdf()->setStreamBound(pCase, m_latticeType.at(f2)->pdf());
	m_latticeType.at(f5)->pdf()->setStreamBound(pCase, m_latticeType.at(f7)->pdf());
	m_latticeType.at(f6)->pdf()->setStreamBound(pCase, m_latticeType.at(f8)->pdf());
	m_latticeType.at(f7)->pdf()->setStreamBound(pCase, m_latticeType.at(f5)->pdf());
	m_latticeType.at(f8)->pdf()->setStreamBound(pCase, m_latticeType.at(f6)->pdf());
	
	
}
D2Q9Incomp_domain::D2Q9Incomp_domain(CLbmCase* pCase, cgsize_t matIndex, cgsize_t nvertex)
{
	m_latticeType.push_back(new D2Q900_incompDomain(pCase,nvertex));
	m_latticeType.push_back(new D2Q901_incompDomain(pCase,nvertex));
	m_latticeType.push_back(new D2Q902_incompDomain(pCase,nvertex));
	m_latticeType.push_back(new D2Q903_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D2Q904_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D2Q905_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D2Q906_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D2Q907_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D2Q908_incompDomain(pCase, nvertex));
}
PdfDomain::LatticeType_t* 
D2Q9Incomp_domain::LatticePdf()
{
	return &m_latticeType;
}

cgsize_t 
D2Q9Incomp_domain::MaterialIndex()
{
	return m_materialIndex;
}

void
D2Q9Incomp_domain::updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t node, PdfDomain::LatticeType_t::size_type f)
{
	//place holder
}

D2Q9Incomp_domain::~D2Q9Incomp_domain()
{
	
}
//LEVEL THREE: invoked from 
//==============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//LbmDomain Abstract class
//Properties:
//1. Base class of the model-specific processes implementation class heirachy
/*--------------------------------------------------------------------------------------------*/
//Static data member
Srt_collider LbmDomain::f_srt;
Srt_collider_improved LbmDomain::f_srt_i;
Trt_collider LbmDomain::f_trt;
Mrt_collider LbmDomain::f_mrt;

void
LbmDomain::lbmInitialize(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index)
{
	Lbm_Initializer functor;
	functor(pdf(), domainVariables, mat_index);
}
void
LbmDomain::lbmCollide(PdfDomain* allPdfs, Domain& domainVariables, CLbmModel::RELAXTIME input, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	input = CLbmModel::MRT;//remove later		
	switch (input)
	{
	case CLbmModel::SRT: f_srt(pdf(), domainVariables, tau, index, f);
		break;
	case CLbmModel::TRT: f_trt(pdf(), domainVariables, tau, index, f);
		break;
	case CLbmModel::MRT: f_mrt(pdf(), allPdfs, domainVariables, tau, index, f);
		break;
	default:
		break;
	}
}

void
LbmDomain::lbmStream(cgsize_t node)
{
	Pdf_streamer functor;
	functor(pdf(), node);
}
void
LbmDomain::lbmWrite(std::fstream& outputFile)
{
	Lbm_writer functor;
	functor(this->pdf(), outputFile);
}

void
LbmDomain::lbmRead(std::fstream& inputFile)
{
	Lbm_reader functor;
	functor(this->pdf(), inputFile);
}

PdfBlock*
LbmDomain::pdf()
{
	return pdfDomain();
}

PdfMaker*
LbmDomain::Impl()
{
	return implDomain();
}

void
LbmDomain::solidObjectNodeBounder(Domain& domainVariables, cgsize_t index, PatchNodeType_t* solidObjectNodes, LbmDomain::BOUND b_type, FORMULA f_type)
{
	cgsize_t neighborNode(0), node(0);

	for(PatchNodeType_t::size_type j=1; j <= solidObjectNodes->size(); j++)
	{
		node = solidObjectNodes->at(j-1);
		neighborNode = pdf()->GetPdfunction()->at(node-1)->m_neighborNum;

		if(neighborNode)
		{

			if(pdf()->GetPdfunction()->at(neighborNode-1)->m_isSolid == false)
			{
				GetFullBBBounder()(index, pdf(),domainVariables,node,b_type,f_type);
			}
		}
	}
}

////fix this method
void 
LbmDomain::nodeBounder(Domain& domainVariables, cgsize_t index, cgsize_t* pPnts, cgsize_t numPnts, BCType_t boundType, LbmDomain::BOUND b_type, FORMULA f_type)
{
	//Distribution function value update flag needs to be false for an update to be performed
	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"now in node bounder \n";
	std::cin.get();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	cgsize_t node(0);
	for(cgsize_t i =0; i < numPnts; i++)
	{
		node = *(pPnts + i);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node # "<<node<<" update status "<<pdf()->GetPdfunction()->at(node-1)->m_valueUpdated << "\n";
		std::cout<<"node # "<<node<<" m_boundtype is "<<pdf()->GetPdfunction()->at(node-1)->m_boundType << "\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		if( pdf()->GetPdfunction()->at(node-1)->m_boundType     == Node::INWARD
			&& pdf()->GetPdfunction()->at(node-1)->m_valueUpdated == false)
			{
				switch (boundType)
					{
						case BCDirichlet:
							{
								if(b_type == VELOCITY || b_type == PRESSURE)
								{
									if(pdf()->GetPdfunction()->at(node-1)->m_isSolid == false)
									{
										
									/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
										std::cout<<"node = "<<node<<" is  a dirichlet node that points inward and not on a wall \n";
										std::cin.get();
									//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
										GetVelBounder()(index, pdf(),domainVariables,node,b_type,f_type);
									}
								}

								else if(b_type == EXTRAPOLATION)
											GetFullBBBounder()(index, pdf(),domainVariables,node,b_type,f_type);
							}
						break;
						case BCWall:
							{
								/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
								std::cout<<"node = "<<node<<" is  a wall node and points inward \n";
								std::cin.get();
								//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
								GetFullBBBounder()(index, pdf(),domainVariables,node,b_type,f_type);
							}
						break;
						case BCGeneral:
							{
								//place holder
							}
						break;
						case BCOutflow:
							{
								//place holder
							}
						break;
				}
		}
		else if(   pdf()->GetPdfunction()->at(node-1)->m_boundType == Node::BURIED
			      && pdf()->GetPdfunction()->at(node-1)->m_valueUpdated == false)
		{
			/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<"node = "<<node<<" is buried \n";
			std::cin.get();
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			BuriedBounder(domainVariables,node, index);
		}
	}
}

void    
LbmDomain::openBCBounder(Domain& domainVariables, cgsize_t index,  cgsize_t* pPnts, cgsize_t numPnts)
{
	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"now in open BC bounder \n";
	std::cin.get();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	cgsize_t node(0);
	//Set b_type and f_type here
	LbmDomain::FORMULA f_type  = LbmDomain::NBC;
	LbmDomain::BOUND   b_type  = LbmDomain::EXTRAPOLATION;
	for(cgsize_t i =0; i < numPnts; i++)
	{
		node = *(pPnts + i);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node # "<<node<<" update status "<<pdf()->GetPdfunction()->at(node-1)->m_valueUpdated << "\n";
		std::cout<<"node # "<<node<<" m_boundtype is "<<pdf()->GetPdfunction()->at(node-1)->m_boundType << "\n";
		std::cout<<"node # "<<node<<" m_periodicNeighborNum is "<<pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum << "\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		
		if( pdf()->GetPdfunction()->at(node-1)->m_boundType     == Node::INWARD
			&& pdf()->GetPdfunction()->at(node-1)->m_valueUpdated == false
			&& pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum == 0)
		{
			/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			std::cout<<"node = "<<node<<" is  an open node that points inward \n";
			std::cin.get();
			//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
			GetOpenBounder()(index, pdf(),domainVariables,node, b_type,f_type);
		}
		else if(   pdf()->GetPdfunction()->at(node-1)->m_boundType == Node::BURIED
			      && pdf()->GetPdfunction()->at(node-1)->m_valueUpdated == false)
		{
			BuriedBounder(domainVariables,node, index);
		}
	}
}

void    
LbmDomain::periodicBCBounder(Domain& domainVariables, cgsize_t index, cgsize_t* pPnts, cgsize_t numPnts)
{
	//apply periodic condition for nodes
	cgsize_t node(0), neighbor(0);
	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"now in periodic bounder \n";
	std::cin.get();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	for(cgsize_t i =0; i < numPnts; i++)
	{
		node = *(pPnts + i);
		neighbor = pdf()->GetPdfunction()->at(node-1)->m_periodicNeighborNum;

		if(pdf()->GetPdfunction()->at(node-1)->m_isSolid == false)
		{
			if(pdf()->GetPdfunction()->at(node-1)->m_neighborNum == 0)
			{
				pdf()->GetPdfunction()->at(neighbor-1)->m_nodeVal = pdf()->GetPdfunction()->at(node-1)->m_nodeVal;
			}
			else;
		}
		else;
	}
}
/*--------------------------------------------------------------------------------------------*/
//D2q9_incompDomain Abstract class
//Properties:
//1. Derived from LbmDomain and base to D2q901_incompDomain...
/*--------------------------------------------------------------------------------------------*/
PdfBlock*
D2q9_incompDomain::maker (PdfMaker* pimpl, CLbmCase* pCase)
{

	 return pimpl->incomp(pCase);
}

PdfBlock*
D2q9_incompDomain::maker (PdfMaker* pimpl, CLbmCase* pCase, cgsize_t nvertex)
{

	 return pimpl->incomp(pCase, nvertex);
}
/*--------------------------------------------------------------------------------------------*/
//D2Q900_incompDomain concrete class
//Properties:
//1. Derived from D2Q9Incomp 
//2. Implementation of process specific to D2q900Incomp objects
/*--------------------------------------------------------------------------------------------*/
D2Q900_incompDomain::D2Q900_incompDomain(CLbmCase* pCase)
{
	pdf_00Domain = maker(Impl(), pCase);
} 
D2Q900_incompDomain::D2Q900_incompDomain(CLbmCase* pCase, cgsize_t nvertex)
{
	pdf_00Domain = maker(Impl(), pCase, nvertex);
} 
Pdf_bounder&
D2Q900_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_setNeighbor&
D2Q900_incompDomain::GetNeighbor()
{
	return neighbor;
}

Pdf_bounder&
D2Q900_incompDomain::GetVelBounder()
{
	return v_boundary;
}

Pdf_bounder&
D2Q900_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

inline PdfBlock* 
D2Q900_incompDomain::pdfDomain() 
{
	return pdf_00Domain;
}

inline PdfMaker* 
D2Q900_incompDomain::implDomain() 
{
	return &impl_00;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q901_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q901_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q901_incompDomain::D2Q901_incompDomain(CLbmCase* pCase)
{
	pdf_01Domain = maker(Impl(), pCase);
}
D2Q901_incompDomain::D2Q901_incompDomain(CLbmCase* pCase, cgsize_t nvertex)
{
	pdf_01Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q901_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q901_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q901_incompDomain::GetVelBounder()
{
	return v_boundary;
}

Pdf_bounder&
D2Q901_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

PdfBlock*  
D2Q901_incompDomain::pdfDomain()
{
	return pdf_01Domain;
}

PdfMaker* 
D2Q901_incompDomain::implDomain()
{
	return &impl_01;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q902_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q902_incompDomain::D2Q902_incompDomain(CLbmCase* pCase) 
{	
	pdf_02Domain = maker(Impl(), pCase);
}
D2Q902_incompDomain::D2Q902_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_02Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q902_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q902_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q902_incompDomain::GetVelBounder()
{
	return v_boundary;
}

Pdf_bounder&
D2Q902_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

PdfBlock*  
D2Q902_incompDomain::pdfDomain()
{
	return pdf_02Domain;
}

PdfMaker* 
D2Q902_incompDomain::implDomain()
{
	return &impl_02;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q903_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q903_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q903_incompDomain::D2Q903_incompDomain(CLbmCase* pCase) 
{	
	pdf_03Domain = maker(Impl(), pCase);
}
D2Q903_incompDomain::D2Q903_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_03Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q903_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q903_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q903_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D2Q903_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
PdfBlock*  
D2Q903_incompDomain::pdfDomain()
{
	return pdf_03Domain;
}

PdfMaker* 
D2Q903_incompDomain::implDomain()
{
	return &impl_03;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q904_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q904_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q904_incompDomain::D2Q904_incompDomain(CLbmCase* pCase) 
{	
	pdf_04Domain = maker(Impl(), pCase);
}
D2Q904_incompDomain::D2Q904_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_04Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q904_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q904_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q904_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D2Q904_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
PdfBlock*  
D2Q904_incompDomain::pdfDomain()
{
	return pdf_04Domain;
}

PdfMaker* 
D2Q904_incompDomain::implDomain()
{
	return &impl_04;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q905_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q905_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q905_incompDomain::D2Q905_incompDomain(CLbmCase* pCase) 
{	
	pdf_05Domain = maker(Impl(), pCase);
}
D2Q905_incompDomain::D2Q905_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_05Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q905_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q905_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q905_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D2Q905_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
void 
D2Q905_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f7(7);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D2Q905_incompDomain::pdfDomain()
{
	return pdf_05Domain;
}

PdfMaker* 
D2Q905_incompDomain::implDomain()
{
	return &impl_05;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q906_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q906_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q906_incompDomain::D2Q906_incompDomain(CLbmCase* pCase) 
{	
	pdf_06Domain = maker(Impl(), pCase);
}
D2Q906_incompDomain::D2Q906_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_06Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q906_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q906_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q906_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D2Q906_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
void 
D2Q906_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f8(8);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D2Q906_incompDomain::pdfDomain()
{
	return pdf_06Domain;
}

PdfMaker* 
D2Q906_incompDomain::implDomain()
{
	return &impl_06;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q907_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q907_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q907_incompDomain::D2Q907_incompDomain(CLbmCase* pCase) 
{	
	pdf_07Domain = maker(Impl(), pCase);
}
D2Q907_incompDomain::D2Q907_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_07Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q907_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q907_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q907_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D2Q907_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
void 
D2Q907_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f5(5);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D2Q907_incompDomain::pdfDomain()
{
	return pdf_07Domain;
}

PdfMaker* 
D2Q907_incompDomain::implDomain()
{
	return &impl_07;
}

/*--------------------------------------------------------------------------------------------*/
//D2Q908_incompDomain concrete class
//Properties:
//1. Derived from D2q9_incompDomain 
//2. Implementation of process specific to D2Q908_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D2Q908_incompDomain::D2Q908_incompDomain(CLbmCase* pCase) 
{	
	pdf_08Domain = maker(Impl(), pCase);
}
D2Q908_incompDomain::D2Q908_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_08Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D2Q908_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D2Q908_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D2Q908_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D2Q908_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
void 
D2Q908_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f6(6);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D2Q908_incompDomain::pdfDomain()
{
	return pdf_08Domain;
}

PdfMaker* 
D2Q908_incompDomain::implDomain()
{
	return &impl_08;
}

//============================================================================================//
//FUNCTORS
void 
Lbm_Initializer::operator() (PdfBlock* pblock, Domain& domainVariables, cgsize_t mat_index)
	{
		//check whether it is improved collision
		pblock->initialize_Pdf(domainVariables, mat_index);
	}
void 
Lbm_writer::operator() (PdfBlock* pblock, std::fstream& outputFile)
	{
		//check whether it is improved collision
		pblock->lbm_write(outputFile);
	}

void 
Lbm_reader::operator() (PdfBlock* pblock, std::fstream& inputFile)
	{
		//check whether it is improved collision
		pblock->lbm_read(inputFile);
	}
void 
Srt_collider::operator() (PdfBlock* pblock, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
	{
		//check whether it is improved collision
		pblock->srt_collide(domainVariables, tau, index, f);
	}

void 
Srt_collider_improved::operator() (PdfBlock* pblock, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index)
	{
		//check whether it is improved collision
		pblock->srt_collide_improved(domainVariables, tau, index);
	}

void 
Trt_collider::operator() (PdfBlock* pblock, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
	{
		pblock->trt_collide(domainVariables, tau, index, f);
	}

void 
Mrt_collider::operator() (PdfBlock* pblock, PdfDomain* allPdfs, Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
	{
		pblock->mrt_collide(allPdfs, domainVariables, tau, index, f);
	}

void
Pdf_streamer::operator() (PdfBlock* pblock, cgsize_t node)
{
	pblock->stream_domain(node);
}

/*Functors derived from Pdf_SetNeighbor
void 
Pdf_setNeighbor::set_Neighbor(PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index, cgsize_t f_opp)
{
	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}*/
//Functors derived from Pdf_Bounder
/*------------------------------------------------------------------------------------------------------------------------------*/
Node::NodeValueType_t
Pdf_bounder::moving_BB(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node)
{			
	cgsize_t rho(4);
	cgsize_t xVel(0), yVel(1), zVel(2);
	Node::NodeValueType_t u_wall(0.0), v_wall(0.0), w_wall(0.0), temp1(0.0), temp(0.0);
	u_wall = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
	v_wall = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
	w_wall = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;	
	//first obtain node in fluid domain
	cgsize_t  useNode(0);
	useNode  = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;

	temp1   = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetCx() * u_wall + 
						domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetCy() * v_wall + 
						domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetCz() * w_wall ;

	temp    = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(useNode-1)->m_nodeVal
						+ 6 * domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetWk() * 
									domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal   * temp1;
	return temp;
}
Node::NodeValueType_t
Pdf_bounder::halfWay_BB(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node)
{
	Node::NodeValueType_t temp1(0.0), temp2(0.0);
	//store value
	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"node = "<<node<<" now perform halfway bb \n";
	std::cin.get();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	temp1 = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)
				 ->pdf()->GetPdfunction()->at(node-1)->m_nodeVal;
	temp2 = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)
				 ->pdf()->GetPdfunction()->at(node-1)->m_nodeVal;
	domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)
				 ->pdf()->GetPdfunction()->at(node-1)->m_nodeVal = temp1;
	return temp2;

}


Node::NodeValueType_t
Pdf_bounder::extrapolation_BB(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node)
{
	cgsize_t  useNode1(0), useNode2(0);
	Node::NodeValueType_t temp2(0.0), c(2.0);
	
	
	useNode1  = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	useNode2  = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(useNode1-1)->m_connectTo;
	//store value
	
	if(useNode2 != 0)
	{
		temp2 = c * (domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)
						->pdf()->GetPdfunction()->at(useNode1-1)->m_nodeVal) -
						(domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)
						->pdf()->GetPdfunction()->at(useNode2-1)->m_nodeVal);
	}
	else
	{
		temp2 = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)
						->pdf()->GetPdfunction()->at(useNode1-1)->m_nodeVal;
	}	
	return temp2;
}
Node::NodeValueType_t
Pdf_bounder::open_NBC(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node)
{
	Node::NodeValueType_t temp(0.0);
	cgsize_t  useNode(0);
	useNode = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	temp    =  domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(useNode-1)->m_nodeVal;
	return temp;
}
Node::NodeValueType_t
Pdf_bounder::open_CBC(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, Node::NodeValueType_t vel)
{
	Node::NodeValueType_t temp(0.0);
	cgsize_t  useNode(0);
	useNode = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	temp    = (domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal - 
		   vel * domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(useNode-1)->m_nodeVal) / 
		  (vel + 1.0);
	return temp;
}
Node::NodeValueType_t 
D2Q9_bounder::ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
													Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound)
{
	//cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);
	//cgsize_t  rho(4), xVel(0), yVel(1), zVel(2);
	Node::NodeValueType_t density(0.0), velocity(0.0), temp(0.0);
	
	//place holder
		
	temp = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal + 
				 c*density*velocity;
	return temp;
}

Node::NodeValueType_t 
D2Q9_bounder::ZhouHe_diag(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t fI, cgsize_t fII, cgsize_t node,
Node::NodeValueType_t c1, Node::NodeValueType_t c2, Node::NodeValueType_t c3,  cgsize_t xvel, cgsize_t yvel, LbmDomain::BOUND b_type, Node::LBMBOUND bound)
{
	cgsize_t rho(4);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8);
	cgsize_t xVel(0), yVel(1);
	Node::NodeValueType_t u_wall(0.0), v_wall(0.0), density(0.0), temp(0.0);
	u_wall = domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
	v_wall = domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;

	if     (b_type == LbmDomain::PRESSURE)
		density = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
	else if(b_type == LbmDomain::VELOCITY)
	{
		if( bound    == Node::WEST
			||bound    == Node::WS
			||bound    == Node::WN)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, node, c, xVel, f0, f2, f4, f3, f6, f7);
			}

		else if(bound == Node::EAST
			||    bound == Node::ES
			||    bound == Node::EN)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, node, c, xVel, f0, f2, f4, f1, f5, f8);
			}

		else if(bound == Node::SOUTH)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, node, c, yVel, f0, f1, f3, f4, f7, f8);
			}
		else if(bound == Node::NORTH)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, node, c, yVel, f0, f1, f3, f2, f6, f5);
			}
	}

	temp =		(c1 * density *u_wall) + (c2 * density * v_wall)
					+ domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal
					+ c3 * ( domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal
					       - domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal ); 
	return temp;
}

void
D2Q900_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//NULL BODY	
}
void
D2Q900_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D2Q900_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q900_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D2Q901_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f3(3);
	
}
void
D2Q901_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d2q9(mat_index); 
	cgsize_t f1(1), f3(3), xvel(0);
	Node::NodeValueType_t  c = (2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d2q9, f1, f3, node, c, xvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q901_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q901_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index); 
		cgsize_t f1(1), f3(3);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f1, f3, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q902_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f4(4);
	
}

void
D2Q902_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d2q9(mat_index);
	cgsize_t f2(2), f4(4), yvel(1);
	Node::NodeValueType_t  c = (2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d2q9, f2, f4, node, c, yvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q902_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q902_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f2(2), f4(4);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f2, f4, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q903_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f1(1);
	
}

void
D2Q903_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d2q9(mat_index);
	cgsize_t f1(1), f3(3), xvel(0);
	Node::NodeValueType_t  c = (-2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d2q9, f3, f1, node, c, xvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D2Q903_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q903_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f1(1), f3(3);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f3, f1, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q904_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f2(2);
	
}

void
D2Q904_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d2q9(mat_index);
	cgsize_t f2(2), f4(4), yvel(1);
	Node::NodeValueType_t  c = (-2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d2q9, f4, f2, node, c, yvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D2Q904_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q904_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f2(2), f4(4);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f4, f2, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q905_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f7(7);
	
}

void
D2Q905_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1);
	cgsize_t d2q9(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f7(7);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WS )
			{
				c1 = (1.0/6.0), c2 = (1.0/2.0), c3 = (-1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f5, f7, f2, f4, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH)
			{
				c1 = (1.0/2.0), c2 = (1.0/6.0), c3 = (-1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f5, f7, f1, f3, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
}
void
D2Q905_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q905_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f5(5), f7(7);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f5, f7, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D2Q906_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f8(8);
	
}

void
D2Q906_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1);
	cgsize_t d2q9(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f6(6), f8(8);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ES )
			{
				c1 = (-1.0/6.0), c2 = (1.0/2.0), c3 = (-1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f6, f8, f2, f4, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH)
			{
				c1 = (-1.0/2.0), c2 = (1.0/6.0), c3 = (1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f6, f8, f1, f3, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
}
void
D2Q906_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q906_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f6(6), f8(8);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f6, f8, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q907_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f5(5);
	
}

void
D2Q907_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1);
	cgsize_t d2q9(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f7(7);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EN )
			{
				c1 = (-1.0/6.0), c2 = (-1.0/2.0), c3 = (1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f7, f5, f2, f4, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH)
			{
				c1 = (-1.0/2.0), c2 = (-1.0/6.0), c3 = (1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f7, f5, f1, f3, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
}
void
D2Q907_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q907_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f5(5), f7(7);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f7, f5, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D2Q908_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//cgsize_t f6(6);
	
}

void
D2Q908_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1);
	cgsize_t d2q9(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f6(6), f8(8);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WN )
			{
				c1 = (1.0/6.0), c2 = (-1.0/2.0), c3 = (1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f8, f6, f2, f4, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH)
			{
				c1 = (1.0/2.0), c2 = (-1.0/6.0), c3 = (-1.0/2.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables, d2q9, f8, f6, f1, f3, node, 
																												c1, c2, c3, xVel, yVel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}

}
void
D2Q908_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D2Q908_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	if(f_type == LbmDomain::HALF_BB)
	{
		cgsize_t d2q9(mat_index);
		cgsize_t f6(6), f8(8);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d2q9, f8, f6, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}


