#include "CLbmBlock3D.h"
#include "CLbmSolve3D.h"
#include <cmath>
#include <iostream>

/*--------------------------------------------------------------------------------------------*/
//Lbm_D3QN Abstract class
//Properties:
//1. Derived from PdfBlock and base to Incomp_D3Q15...
/*--------------------------------------------------------------------------------------------*/

inline Node::NodeValueType_t 
Lbm_D3QN::GetWeight()
{
	return weightFactor();
}

void 
Lbm_D3QN::setNodeCoord(CLbmCase* pCase)
{
	cgsize_t x(0), y(1), z(2);
	// Add the boundary name(s) to boundary nodes
	for(NodeType_t::size_type i = 1; i <= GetPdfunction()->size(); i++)
	{
		GetPdfunction()->at(i-1)->m_xcoord = *(*(pCase->m_grid.m_coord + x)  + (i-1));
		GetPdfunction()->at(i-1)->m_ycoord = *(*(pCase->m_grid.m_coord + y)  + (i-1));
		GetPdfunction()->at(i-1)->m_zcoord = *(*(pCase->m_grid.m_coord + z)  + (i-1));
	}
}

void
Lbm_D3QN::setLinkage(CLbmCase* pCase)
{
	ElementType_t eltType;
	for(CLbmElement::MaterialIndex_t::iterator iter = pCase->m_element.m_materialIndex.begin();
		  iter != pCase->m_element.m_materialIndex.end(); iter++)
	{
		eltType = *(pCase->m_element.m_elementType + *iter);
		switch (eltType)
		{
		case HEXA_8: 
			setHexa_8conn(*iter);
			break;
		default:
			break;
		}
	}
}


Force3D
Lbm_D3QN::GetTotalForce(Domain& domainVariables, cgsize_t mat_index, cgsize_t node, Node::NodeValueType_t local_rho)
{
	Force3D  force;
	cgsize_t fx(1), fy(2), fz(3);
	Node::NodeValueType_t gx(0.0), gy(0.0), gz(0.0);

	domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal += ( gx * local_rho);
	domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal += ( gy * local_rho);
	domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal += ( gz * local_rho);

	force.x = domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal ;
	force.y = domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal ;
	force.z = domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal ;

	/*domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;
	domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;
	domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;*/
	/*if(node == 620000)
	{
		std::cout<<" ****************************\n";
		std::cout<<" force.x = "<<force.x<<"\n";
		std::cout<<" force.y = "<<force.y<<"\n";
		std::cout<<" force.z = "<<force.z<<"\n";
		std::cout<<" ****************************\n";
	}*/
	return force;
}

Velocity3D
Lbm_D3QN::GetComponentVel(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Velocity3D velocity;
	cgsize_t xVel(0), yVel(1), zVel(2);

	velocity.x = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	velocity.y = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	velocity.z = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

	return velocity;
}
Velocity3D
Lbm_D3QN::GetOverallVel(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Velocity3D velocity;
	cgsize_t uoverall_x(5), uoverall_y(6), uoverall_z(7);

	velocity.x = domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	velocity.y = domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	velocity.z = domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

	return velocity;
}
Velocity3D
Lbm_D3QN::ComputeVelocity(Force3D totalForce, Velocity3D velocity, Node::NodeValueType_t localRho, Node::NodeValueType_t tau, cgsize_t M)
{
	Velocity3D velInEquation;
	velInEquation.x = velocity.x + (0.5 * M * tau * totalForce.x / localRho);
	velInEquation.y = velocity.y + (0.5 * M * tau * totalForce.y / localRho);
	velInEquation.z = velocity.z + (0.5 * M * tau * totalForce.z / localRho);
	return velInEquation;
}

Velocity3D
Lbm_D3QN::ComputeVelocity(Force3D totalForce, Velocity3D velocity, Node::NodeValueType_t localRho, cgsize_t M)
{
	//not correct fix!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Just use overall velocity
	Velocity3D velInEquation;
	velInEquation.x = velocity.x + (0.5 * M * totalForce.x / localRho);
	velInEquation.y = velocity.y + (0.5 * M * totalForce.y / localRho);
	velInEquation.z = velocity.z + (0.5 * M * totalForce.z / localRho);
	return velInEquation;
}

Velocity3D
Lbm_D3QN::ComputeImproveVel(Force3D totalForce, Velocity3D velocity, Node::NodeValueType_t kinVisco, Node::NodeValueType_t psi, Node::NodeValueType_t sigma)
{
	Velocity3D velInEquation;
	velInEquation.x = velocity.x + (sigma * totalForce.x / (kinVisco * psi * psi));
	velInEquation.y = velocity.y + (sigma * totalForce.y / (kinVisco * psi * psi));
	velInEquation.z = velocity.z + (sigma * totalForce.z / (kinVisco * psi * psi));
	return velInEquation;
}
//PDFBLOCK HELPER FUNCTION
/*--------------------------------------------------------------------------------------------*/
//Incomp_D3Q15 Abstract class
//Properties:
//1. Derived from Lbm_D3QN and base to D3q1500_incomp...
/*--------------------------------------------------------------------------------------------*/
//Forcing term utilities
Node::NodeValueType_t
Incomp_D3Q15::Forcing_Part1(Force3D totalForce, Node::NodeValueType_t cs_square)
{	
	Node::NodeValueType_t temp = GetCx() * totalForce.x + GetCy() * totalForce.y + GetCz() * totalForce.z;
	return (temp / cs_square);
}

Node::NodeValueType_t
Incomp_D3Q15::Forcing_Part2(Force3D totalForce, Velocity3D vel, Node::NodeValueType_t cs_square)
{	
	Node::NodeValueType_t temp = (vel.x * totalForce.x * GetCx() * GetCx() - vel.x * totalForce.x * cs_square) +
				     (vel.y * totalForce.x * GetCx() * GetCy()) +
				     (vel.z * totalForce.x * GetCx() * GetCz()) +

				     (vel.x * totalForce.y * GetCy() * GetCx()) +
				     (vel.y * totalForce.y * GetCy() * GetCy() - vel.y * totalForce.y * cs_square) +
				     (vel.z * totalForce.y * GetCy() * GetCz()) +

				     (vel.x * totalForce.z * GetCz() * GetCx()) +
				     (vel.y * totalForce.z * GetCz() * GetCy()) +
				     (vel.z * totalForce.z * GetCz() * GetCz() - vel.z * totalForce.z * cs_square) ;
	return (temp / (cs_square * cs_square));
}

//Forcing Term
Node::NodeValueType_t
Incomp_D3Q15::ForcingTerm(Force3D totalForce, Velocity3D vel, Node::NodeValueType_t B_e, Node::NodeValueType_t cs_square)
{
	Node::NodeValueType_t temp = Forcing_Part1(totalForce, cs_square) + Forcing_Part2(totalForce, vel, cs_square); 
	return GetWeight() * temp * B_e;
}

//Equilibrium PDF utilities
Node::NodeValueType_t
Incomp_D3Q15::EquilPDF_Part1(Velocity3D equilVel, Node::NodeValueType_t cs_square)
{	
	Node::NodeValueType_t temp = GetCx() * equilVel.x + GetCy() * equilVel.y + GetCz() * equilVel.z;
	return (temp / cs_square);
}

Node::NodeValueType_t
Incomp_D3Q15::EquilPDF_Part2(Velocity3D equilVel, Node::NodeValueType_t cs_square)
{	
	Node::NodeValueType_t temp = (equilVel.x * equilVel.x * GetCx() * GetCx() - equilVel.x * equilVel.x * cs_square) +
				     (equilVel.y * equilVel.x * GetCx() * GetCy()) +
				     (equilVel.z * equilVel.x * GetCx() * GetCz()) +

				     (equilVel.x * equilVel.y * GetCy() * GetCx()) +
				     (equilVel.y * equilVel.y * GetCy() * GetCy() - equilVel.y * equilVel.y * cs_square) +
				     (equilVel.z * equilVel.y * GetCy() * GetCz()) +

				     (equilVel.x * equilVel.z * GetCz() * GetCx()) +
				     (equilVel.y * equilVel.z * GetCz() * GetCy()) +
				     (equilVel.z * equilVel.z * GetCz() * GetCz() - equilVel.z * equilVel.z * cs_square) ;
	return (temp / (2.0 * cs_square * cs_square));
}

//Equilibrium PDF
Node::NodeValueType_t
Incomp_D3Q15::EquilPDFTerm(Velocity3D equilVel, Node::NodeValueType_t cs_square, Node::NodeValueType_t local_rho)
{
	Node::NodeValueType_t temp = 1.0 + EquilPDF_Part1(equilVel, cs_square) + EquilPDF_Part2(equilVel, cs_square); 
	return GetWeight() * local_rho * temp;
}

void 
Incomp_D3Q15::srt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	cgsize_t rho(4), pseudoP(0);
	NodeType_t::size_type i;
	NodeType_t::size_type end = GetPdfunction()->size();

	Node::NodeValueType_t local_rho(0.0), equilPDF(0.0), forceTerm(0.0), cs_square(1.0/3.0), psi(0.0);
	Velocity3D  equil_vel, force_vel, compo_vel, improved_vel;
	Force3D     total_force;

	//AAAAAAATTTTTTTTTTTEEEEEEEEEEEENNNNNNNNNNNNNNTTTTTTTTIIIIIIIIOOOOOOOOOOOOOOOOOONNNNNNNNNNNNNNN?//
	//TEMPORARILY ADJUST TAU HERE AND AT THE TOP
	tau = 1.0;
	Node::NodeValueType_t B_e = (1.0 - (1.0/(tau * 2.0)));
	Node::NodeValueType_t kinVisco = (tau - 0.5);
	Node::NodeValueType_t sigma = 0.105;

	

#pragma omp for private(i)
	for(i =1; i <= end; i++)
	{
		if(GetPdfunction()->at(i-1)->m_isSolid == false)
		{
			
			local_rho   = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(index)->at(i-1)->m_nodeVal;
			
			psi         = domainVariables.GetDomainTempVariables()->at(pseudoP)->GetVariable()->at(index)->at(i-1)->m_nodeVal;
			
			compo_vel   = GetComponentVel(domainVariables, index, i);
			
			total_force = GetTotalForce(domainVariables, index, i, local_rho);
			
			equil_vel   = GetOverallVel(domainVariables, index, i);
			
			equilPDF    = EquilPDFTerm(equil_vel, cs_square, local_rho);

			force_vel   = GetOverallVel(domainVariables, index, i);
			
			improved_vel = ComputeImproveVel(total_force, force_vel, kinVisco, psi, sigma);
			
			forceTerm   = ForcingTerm(total_force, improved_vel, B_e, cs_square);

			GetPdfunction()->at(i-1)->m_nodeVal =  GetPdfunction()->at(i-1)->m_nodeVal - 
				(GetPdfunction()->at(i-1)->m_nodeVal - equilPDF) / tau +  forceTerm;

			
		}
	}

}


void 
Incomp_D3Q15::srt_collide_improved(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index)
{
	//place holder
}


void 
Incomp_D3Q15::trt_collide(Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	
}
/*--------------------------------------------------------------------------------MRT----------------------------------------------------------------------------*/
void 
Incomp_D3Q15::mrt_collide(PdfDomain* allPdfs, Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	NodeType_t::size_type i;
	NodeType_t::size_type end = GetPdfunction()->size();
#pragma omp for private(i)  
	for(i =1; i <= end; i++)
	{
		if(GetPdfunction()->at(i-1)->m_isSolid == false)
		{
			//transform pdf first
			//GetPdfunction()->at(i-1)->m_nodeVal   -=   transformPDF(domainVariables, index, i);
			//perform mrt collision
			GetPdfmoment()->at(i-1)->m_nodeVal     =   compute_collision(domainVariables, allPdfs, index, tau, i);
		}
		
	}
}
/*--------------------------------------------------------------------------------MRT----------------------------------------------------------------------------*/
/*--------------------------------------------------------------------------------------------*/
//Incomp_D3Q19 Abstract class
//Properties:
//1. Derived from Lbm_D3QN and base to D3q1900_incomp...
/*--------------------------------------------------------------------------------------------*/
//Member function that implements collision
void 
Incomp_D3Q19::srt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	//place holder
}


void 
Incomp_D3Q19::srt_collide_improved(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index)
{
	//place holder
}


void 
Incomp_D3Q19::trt_collide(Domain& domainVariables, Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	//place holder
}
void 
Incomp_D3Q19::mrt_collide(PdfDomain* allPdfs, Domain& domainVariables,  Node::NodeValueType_t tau, cgsize_t index, cgsize_t f)
{
	//place holder
}





//============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//D3q1500_incomp concrete class
//Properties:
//1. Derived from D3q15_00. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1500_incomp::D3q1500_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0), cz(0.0)
{
	weight = 2.0/9.0;
	setBlock(pCase);
}

D3q1500_incomp::D3q1500_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(0.0), cz(0.0)
{
	weight = 2.0/9.0;
	setBlock(nvertex);
}
void
D3q1500_incomp::setHexa_8conn(cgsize_t i){
	//NULL BODY
}
//copy ctor: using default copy ctor

//dtor
D3q1500_incomp::~D3q1500_incomp()
{

}
//functions that reads PdfBlock fields

PdfBlock::NodeType_t* 
D3q1500_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1500_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1500_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1500_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1500_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1500_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1500_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1500_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1500_incomp::c_z()
{
	return cz;
}
Node::NodeValueType_t
D3q1500_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = GetWk()*local_rho - 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1500_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	
	return source;
}
Node::NodeValueType_t
D3q1500_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq;
	cgsize_t rho(4);
	moment_eq = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;	

	return moment_eq;
}

Node::NodeValueType_t
D3q1500_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t   M0(1.0), M1(1.0), M2(1.0), M3(1.0), M4(1.0), M5(1.0), M6(1.0),
		 		M7(1.0), M8(1.0), M9(1.0), M10(1.0), M11(1.0), M12(1.0), M13(1.0), M14(1.0);

	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;
	/*if(node == 620000)
	{

		std::cout<<" f0_old = "<<allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f1_old = "<<allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f2_old = "<<allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f3_old = "<<allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f4_old = "<<allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f5_old = "<<allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f6_old = "<<allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f7_old = "<<allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f8_old = "<<allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f9_old = "<<allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f10_old = "<<allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f11_old = "<<allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f12_old = "<<allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f13_old = "<<allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" f14_old = "<<allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
	}*/
	return moment;
}

Node::NodeValueType_t
D3q1500_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A0(1.0);
	
	mf           = compute_moment(allPdfs, node);	
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	//sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" in f0 \n";
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A0* (mf - mf_eq) - (I - 0.5*A0) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}


Node::NodeValueType_t
D3q1500_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment = 0.0;

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1501_incomp concrete class
//Properties:
//1. Derived from D3q15_01. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1501_incomp::D3q1501_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(0.0), cz(0.0)
{
	//set block info
	weight = 1.0/9.0;
	setBlock(pCase);
}
D3q1501_incomp::D3q1501_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(0.0), cz(0.0)
{
	//set block info
	weight = 1.0/9.0;
	setBlock(nvertex);
}
void 
D3q1501_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeF, nodeB, tonodeF, tonodeB;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeF   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeF = *(*(m_pCase->m_element.m_conn + i) + j+1);
		nodeB   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonodeB = *(*(m_pCase->m_element.m_conn + i) + j+2);
			if(GetPdfunction()->at(tonodeF-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeF-1)->m_connectTo = nodeF;
				GetPdfunction()->at(tonodeF-1)->m_linkUpdated   = true;
			}
			if (GetPdfunction()->at(tonodeB-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeB-1)->m_connectTo = nodeB;
				GetPdfunction()->at(tonodeB-1)->m_linkUpdated   = true;
			}
	}
}

//copy ctor: using default copy ctor
//dtor
D3q1501_incomp::~D3q1501_incomp()
{

}

//function that read fields

PdfBlock::NodeType_t* 
D3q1501_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1501_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1501_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1501_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1501_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1501_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1501_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1501_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1501_incomp::c_z()
{
	return cz;
}

void 
D3q1501_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::WEST)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WSB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WNB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WST ||
						GetPdfunction()->at(node-1)->m_bound == Node::WNT ||
						GetPdfunction()->at(node-1)->m_bound == Node::WS  ||
						GetPdfunction()->at(node-1)->m_bound == Node::WN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::WT  ||
						GetPdfunction()->at(node-1)->m_bound == Node::WB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}

Node::NodeValueType_t
D3q1501_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1501_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	source = 3.0 * GetWk() * (force.x + (2.0*vel.x*force.x) - (vel.y*force.y) - (vel.z*force.z));
	
	return source;
}
Node::NodeValueType_t
D3q1501_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), temp(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	
	temp = 	pow(vel.x, 2) + pow(vel.y, 2) + pow(vel.z, 2);

	moment_eq =  local_rho * (-1.0 + temp);
	//moment_eq =  local_rho * (1.0 + temp);
		   	 
	return moment_eq;
}

Node::NodeValueType_t
D3q1501_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(-2.0), M1(-1.0), M2(-1.0), M3(-1.0), M4(-1.0), M5(-1.0), M6(-1.0),
		 		M7(1.0), M8(1.0), M9(1.0), M10(1.0), M11(1.0), M12(1.0), M13(1.0), M14(1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1501_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A1(1.0);//A1 related to bulk viscosity
	
	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A1* (mf - mf_eq) - (I - 0.5*A1) * sourceMoment;
 	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1501_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	sourceMoment   = 2.0 * (force.x * vel.x + force.y * vel.y + force.z * vel.z);
	
	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1502_incomp concrete class
//Properties:
//1. Derived from D3q15_02. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1502_incomp::D3q1502_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(0.0), cz(0.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(pCase);
}
D3q1502_incomp::D3q1502_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(0.0), cz(0.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(nvertex);
}
void 
D3q1502_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeF, nodeB, tonodeF, tonodeB;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeF   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeF = *(*(m_pCase->m_element.m_conn + i) + j+3);
		nodeB   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeB = *(*(m_pCase->m_element.m_conn + i) + j);
			if(GetPdfunction()->at(tonodeF-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeF-1)->m_connectTo   = nodeF;
				GetPdfunction()->at(tonodeF-1)->m_linkUpdated = true;
			}
			if(GetPdfunction()->at(tonodeB-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeB-1)->m_connectTo   = nodeB;
				GetPdfunction()->at(tonodeB-1)->m_linkUpdated = true;
			}
	}
}

//copy ctor: using default copy ctor
//dtor
D3q1502_incomp::~D3q1502_incomp()
{

}

//function that read fields

PdfBlock::NodeType_t* 
D3q1502_incomp::pdfunction()
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1502_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1502_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1502_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1502_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1502_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1502_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1502_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1502_incomp::c_z()
{
	return cz;
}

void 
D3q1502_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::EAST)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::ESB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENB ||
						GetPdfunction()->at(node-1)->m_bound == Node::EST ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENT ||
						GetPdfunction()->at(node-1)->m_bound == Node::ES  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ET  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
Node::NodeValueType_t
D3q1502_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1502_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	source = 3.0 * GetWk() * (-1.0*force.x + (2.0*vel.x*force.x) - (vel.y*force.y) - (vel.z*force.z));
	
	return source;
}
Node::NodeValueType_t
D3q1502_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), temp(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	temp = 	pow(vel.x, 2) + pow(vel.y, 2) + pow(vel.z, 2);

	moment_eq = local_rho * (1.0 - 5.0*temp);
	//moment_eq = -5.0*local_rho * (3.0 + temp);
		   	 
	return moment_eq;
}

Node::NodeValueType_t
D3q1502_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(16.0), M1(-4.0), M2(-4.0), M3(-4.0), M4(-4.0), M5(-4.0), M6(-4.0),
		 		M7(1.0), M8(1.0), M9(1.0), M10(1.0), M11(1.0), M12(1.0), M13(1.0), M14(1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1502_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A2(1.0);//free parameters that does not relate to the hydrodynamics (artibrary)

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A2* (mf - mf_eq) - (I - 0.5*A2) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}


Node::NodeValueType_t
D3q1502_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = -10.0 * (force.x * vel.x + force.y * vel.y + force.z * vel.z);
	
	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1503_incomp concrete class
//Properties:
//1. Derived from D3q15_03. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1503_incomp::D3q1503_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(1.0), cz(0.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(pCase);
}
D3q1503_incomp::D3q1503_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(1.0), cz(0.0)
{
	weight = 1.0/9.0;
	//set block info
	setBlock(nvertex);
}
void 
D3q1503_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeO, nodeI, tonodeO, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+3);
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+2);
			if(GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeO-1)->m_connectTo   = nodeO;
				GetPdfunction()->at(tonodeO-1)->m_linkUpdated = true;
			}
			if(GetPdfunction()->at(tonodeI-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeI-1)->m_connectTo   = nodeI;
				GetPdfunction()->at(tonodeI-1)->m_linkUpdated = true;
			}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1503_incomp::~D3q1503_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q1503_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1503_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1503_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1503_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1503_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1503_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1503_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1503_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1503_incomp::c_z()
{
	return cz;
}

void 
D3q1503_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::SOUTH)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WSB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ESB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WST ||
						GetPdfunction()->at(node-1)->m_bound == Node::EST ||
						GetPdfunction()->at(node-1)->m_bound == Node::WS  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ES  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ST  ||
						GetPdfunction()->at(node-1)->m_bound == Node::SB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
Node::NodeValueType_t
D3q1503_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1503_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	source = 3.0 * GetWk() * (force.y + (2.0*vel.y*force.y) - (vel.x*force.x) - (vel.z*force.z));
	/*if(node == 620000)
	{
		std::cout<<" fx   = "<<force.x<<"\n";
		std::cout<<" fy   = "<<force.y<<"\n";
		std::cout<<" fz   = "<<force.z<<"\n";
		std::cout<<" fsource = "<<source<<"\n";
	}*/
	return source;
}
Node::NodeValueType_t
D3q1503_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	
	moment_eq = local_rho * vel.x;

	return moment_eq;
}

Node::NodeValueType_t
D3q1503_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t	M0(0), M1(1.0), M2(-1.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(-1.0), M9(1.0), M10(-1.0), M11(1.0), M12(-1.0), M13(1.0), M14(-1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;	

	return moment;
}

Node::NodeValueType_t
D3q1503_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A3(1.0);//related to conservation equations

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node); 
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A3* (mf - mf_eq) - (I - 0.5*A3) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1503_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);

	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Force3D force = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = force.x;

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1504_incomp concrete class
//Properties:
//1. Derived from D3q15_04. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1504_incomp::D3q1504_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0), cy(-1.0), cz(0.0)
{
	weight = 1.0/9.0;;
	//set block info
	setBlock(pCase);
}
D3q1504_incomp::D3q1504_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0), cy(-1.0), cz(0.0)
{
	weight = 1.0/9.0;;
	//set block info
	setBlock(nvertex);
}
void 
D3q1504_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeO, nodeI, tonodeO, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j);
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+1);
			if(GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeO-1)->m_connectTo   = nodeO;
				GetPdfunction()->at(tonodeO-1)->m_linkUpdated = true;
			}
			if(GetPdfunction()->at(tonodeI-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeI-1)->m_connectTo   = nodeI;
				GetPdfunction()->at(tonodeI-1)->m_linkUpdated = true;
			}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1504_incomp::~D3q1504_incomp()
{

} 

//function that read fields

PdfBlock::NodeType_t* 
D3q1504_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1504_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1504_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1504_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1504_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1504_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1504_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1504_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1504_incomp::c_z()
{
	return cz;
}

void 
D3q1504_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::NORTH)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WNB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WNT ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENT ||
						GetPdfunction()->at(node-1)->m_bound == Node::WN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NT  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
Node::NodeValueType_t
D3q1504_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1504_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	source = 3.0 * GetWk() * (-1.0*force.y + (2.0*vel.y*force.y) - (vel.x*force.x) - (vel.z*force.z));
	
	return source;
}
Node::NodeValueType_t
D3q1504_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	
	moment_eq = -(7.0/3.0) * local_rho * vel.x;	

	return moment_eq;
}

Node::NodeValueType_t
D3q1504_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(-4.0), M2(4.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(-1.0), M9(1.0), M10(-1.0), M11(1.0), M12(-1.0), M13(1.0), M14(-1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1504_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{	
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A4(1.0);//free parameters that does not relate to the hydrodynamics (artibrary)

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A4 * (mf - mf_eq) - (I - 0.5*A4) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1504_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);

	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Force3D force = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = (-7.0/3.0) * force.x;

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1505_incomp concrete class
//Properties:
//1. Derived from D3q15_05. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1505_incomp::D3q1505_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0), cz(1.0)
{
	weight = 1.0/9.0;;
	//set block info
	setBlock(pCase);
}
D3q1505_incomp::D3q1505_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(0.0), cz(1.0)
{
	weight = 1.0/9.0;;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1505_incomp::setHexa_8conn(cgsize_t i){
	const cgsize_t hexa = 8;
	const cgsize_t quad = 4;
	cgsize_t temp(0);
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		temp = hexa -1;
		for(cgsize_t k=temp; k >= quad; k--)
		{
			node    = *(*(m_pCase->m_element.m_conn + i) + j+(k-quad));
			tonode  = *(*(m_pCase->m_element.m_conn + i) + j+k);
		
			if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonode-1)->m_connectTo     = node;
				GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
			}
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1505_incomp::~D3q1505_incomp()
{

} 

//function that read fields

PdfBlock::NodeType_t* 
D3q1505_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1505_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1505_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1505_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1505_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1505_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1505_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1505_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1505_incomp::c_z()
{
	return cz;
}

void 
D3q1505_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WNB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WSB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ESB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WB  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EB  ||
						GetPdfunction()->at(node-1)->m_bound == Node::SB  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
Node::NodeValueType_t
D3q1505_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1505_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	source = 3.0 * GetWk() * (force.z + (2.0*vel.z*force.z) - (vel.x*force.x) - (vel.y*force.y));
	/*if(node == 620000)
	{
		std::cout<<" fx   = "<<force.x<<"\n";
		std::cout<<" fy   = "<<force.y<<"\n";
		std::cout<<" fz   = "<<force.z<<"\n";
		std::cout<<" fsource = "<<source<<"\n";
	}*/
	return source;
}
Node::NodeValueType_t
D3q1505_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	
	moment_eq = local_rho * vel.y;	

	return moment_eq;
}

Node::NodeValueType_t
D3q1505_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(1.0), M4(-1.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(1.0), M9(-1.0), M10(-1.0), M11(1.0), M12(1.0), M13(-1.0), M14(-1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;	

	return moment;
}

Node::NodeValueType_t
D3q1505_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A5(1.0);//related to conservation equations

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A5 * (mf - mf_eq) - (I - 0.5*A5) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1505_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);

	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Force3D force = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = force.y;

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1506_incomp concrete class
//Properties:
//1. Derived from D3q15_06. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1506_incomp::D3q1506_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/9.0;;
	//set block info
	setBlock(pCase);
}
D3q1506_incomp::D3q1506_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/9.0;;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1506_incomp::setHexa_8conn(cgsize_t i){
	const cgsize_t hexa = 8;
	const cgsize_t quad = 4;
	cgsize_t temp(0);
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		temp = hexa -1;
		for(cgsize_t k=temp; k >= quad; k--)
		{
			node    = *(*(m_pCase->m_element.m_conn + i) + j+k);
			tonode  = *(*(m_pCase->m_element.m_conn + i) + j+(k-quad));
		
			if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonode-1)->m_connectTo     = node;
				GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
			}
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1506_incomp::~D3q1506_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q1506_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1506_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1506_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1506_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1506_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1506_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1506_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1506_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1506_incomp::c_z()
{
	return cz;
}

void 
D3q1506_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::TOP)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WNT ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENT ||
						GetPdfunction()->at(node-1)->m_bound == Node::WST ||
						GetPdfunction()->at(node-1)->m_bound == Node::EST ||
						GetPdfunction()->at(node-1)->m_bound == Node::WT  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ET  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ST  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NT
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
Node::NodeValueType_t
D3q1506_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1506_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	source = 3.0 * GetWk() * (-1.0*force.z + (2.0*vel.z*force.z) - (vel.x*force.x) - (vel.y*force.y));
	
	return source;
}
Node::NodeValueType_t
D3q1506_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);	

	moment_eq = -(7.0/3.0) * local_rho * vel.y;	

	return moment_eq;
}

Node::NodeValueType_t
D3q1506_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(-4.0), M4(4.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(1.0), M9(-1.0), M10(-1.0), M11(1.0), M12(1.0), M13(-1.0), M14(-1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1506_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A6(1.0);//free parameters that does not relate to the hydrodynamics (artibrary)

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node); 
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A6 * (mf - mf_eq) - (I - 0.5*A6) * sourceMoment ;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1506_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);

	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Force3D force = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = (-7.0/3.0) * force.y;

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1507_incomp concrete class
//Properties:
//1. Derived from D3q15_07. Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1507_incomp::D3q1507_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q1507_incomp::D3q1507_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1507_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+6);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1507_incomp::~D3q1507_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q1507_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1507_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1507_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1507_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1507_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1507_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1507_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1507_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1507_incomp::c_z()
{
	return cz;
}

void 
D3q1507_incomp::setNodeBType(cgsize_t  node)
{
	if(
				GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			||	GetPdfunction()->at(node-1)->m_bound == Node::WB
			||	GetPdfunction()->at(node-1)->m_bound == Node::SB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
			||	GetPdfunction()->at(node-1)->m_bound == Node::WSB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q1507_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1507_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = force.x + force.y + force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (vel.y*force.x + vel.z*force.x + vel.x*force.y + vel.z*force.y + vel.x*force.z + vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	/*if(node == 620000)
	{
		std::cout<<" fx   = "<<force.x<<"\n";
		std::cout<<" fy   = "<<force.y<<"\n";
		std::cout<<" fz   = "<<force.z<<"\n";
		std::cout<<" fsource = "<<source<<"\n";
	}*/
	return source;
}
Node::NodeValueType_t
D3q1507_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	

	moment_eq = local_rho * vel.z;	

	return moment_eq;
}

Node::NodeValueType_t
D3q1507_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(1.0), M6(-1.0),
		 		M7(1.0), M8(1.0), M9(1.0), M10(1.0), M11(-1.0), M12(-1.0), M13(-1.0), M14(-1.0);
	

	/*if(node == 620000)
	{
		std::cout<<" moment before = "<<moment<<"\n";
		std::cout<<" M0 = "<<M0<<" f0 = "<<allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M0 * f0 = "<<M0*allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M1 = "<<M1<<" f1 = "<<allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M1 * f1 = "<<M1*allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M2 = "<<M2<<" f2 = "<<allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M2 * f2 = "<<M2*allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M3 = "<<M3<<" f3 = "<<allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M3 * f3 = "<<M3*allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M4 = "<<M4<<" f4 = "<<allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M4 * f4 = "<<M4*allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M5 = "<<M5<<" f5 = "<<allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M5 * f5 = "<<M5*allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M6 = "<<M6<<" f6 = "<<allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M6 * f6 = "<<M6*allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M7 = "<<M7<<" f7 = "<<allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M7 * f7 = "<<M7*allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M8 = "<<M8<<" f8 = "<<allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M8 * f8 = "<<M8*allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M9 = "<<M9<<" f9 = "<<allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M9 * f9 = "<<M9*allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M10 = "<<M10<<" f10 = "<<allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M10 * f10 = "<<M10*allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M11 = "<<M11<<" f11 = "<<allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M11 * f11 = "<<M11*allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M12 = "<<M12<<" f12 = "<<allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M12 * f12 = "<<M12*allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M13 = "<<M13<<" f13 = "<<allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M13 * f13 = "<<M13*allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cout<<" M14 = "<<M14<<" f14 = "<<allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<" M14 * f14 = "<<M14*allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
		std::cin.get();
	}*/

	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1507_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A7(1.0);//related to conservation equations

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A7 * (mf - mf_eq) - (I - 0.5*A7) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1507_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);

	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Force3D force = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = force.z;

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1508_incomp concrete class
//Properties:
//1. Derived from D3q15_08 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q15014_incomp::D3q15014_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(-1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q15014_incomp::D3q15014_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(-1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q15014_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+6);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q15014_incomp::~D3q15014_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q15014_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q15014_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q15014_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q15014_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q15014_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q15014_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q15014_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q15014_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q15014_incomp::c_z()
{
	return cz;
}

void 
D3q15014_incomp::setNodeBType(cgsize_t  node)
{
	if(    GetPdfunction()->at(node-1)->m_bound == Node::EAST
			|| GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::TOP
			|| GetPdfunction()->at(node-1)->m_bound == Node::ET
			|| GetPdfunction()->at(node-1)->m_bound == Node::NT
			|| GetPdfunction()->at(node-1)->m_bound == Node::EN
			|| GetPdfunction()->at(node-1)->m_bound == Node::ENT
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q15014_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q15014_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node) //fix later
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = force.x + force.y + force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (vel.y*force.x + vel.z*force.x + vel.x*force.y + vel.z*force.y + vel.x*force.z + vel.y*force.z);

	source = 3.0 * GetWk() * ( -1.0*term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q15014_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq;
	moment_eq = 0.0;	

	return moment_eq;
}

Node::NodeValueType_t
D3q15014_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment;

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(-1.0), M9(-1.0), M10(1.0), M11(-1.0), M12(1.0), M13(1.0), M14(-1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;	

	return moment;
}

Node::NodeValueType_t
D3q15014_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A14(1.0);//free parameters that does not relate to the hydrodynamics (artibrary)
	
	mf           = compute_moment(allPdfs, node);
	//mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	//sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A14 * (mf - mf_eq) - (I - 0.5*A14) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q15014_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment = 0.0;
	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1509_incomp concrete class
//Properties:
//1. Derived from D3q15_09 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q15011_incomp::D3q15011_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q15011_incomp::D3q15011_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q15011_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+4);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+2);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q15011_incomp::~D3q15011_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q15011_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q15011_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q15011_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q15011_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q15011_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q15011_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q15011_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q15011_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q15011_incomp::c_z()
{
	return cz;
}

void 
D3q15011_incomp::setNodeBType(cgsize_t  node)
{
	if(    GetPdfunction()->at(node-1)->m_bound == Node::WEST
			|| GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::TOP
			|| GetPdfunction()->at(node-1)->m_bound == Node::WT
			|| GetPdfunction()->at(node-1)->m_bound == Node::ST
			|| GetPdfunction()->at(node-1)->m_bound == Node::WS
			|| GetPdfunction()->at(node-1)->m_bound == Node::WST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q15011_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q15011_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = force.x + force.y - force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (vel.y*force.x - vel.z*force.x + vel.x*force.y - vel.z*force.y - vel.x*force.z - vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q15011_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);

	moment_eq = local_rho  * vel.x * vel.y;	

	return moment_eq;
}

Node::NodeValueType_t
D3q15011_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(-1.0), M9(-1.0), M10(1.0), M11(1.0), M12(-1.0), M13(-1.0), M14(1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q15011_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A11(1.0/tau);//related to kinematic viscosity

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node); 
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A11 * (mf - mf_eq) - (I - 0.5*A11) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q15011_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = force.y * vel.x + force.x * vel.y;

	return sourceMoment;	
}
/*--------------------------------------------------------------------------------------------*/
//D3q15010_incomp concrete class
//Properties:
//1. Derived from D3q15_010 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q15010_incomp::D3q15010_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(-1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q15010_incomp::D3q15010_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(-1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q15010_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+4);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q15010_incomp::~D3q15010_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q15010_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q15010_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q15010_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q15010_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q15010_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q15010_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q15010_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q15010_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q15010_incomp::c_z()
{
	return cz;
}

void 
D3q15010_incomp::setNodeBType(cgsize_t  node)
{
	if(		 GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			|| GetPdfunction()->at(node-1)->m_bound == Node::EAST
			|| GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::EB
			|| GetPdfunction()->at(node-1)->m_bound == Node::NB
			|| GetPdfunction()->at(node-1)->m_bound == Node::EN
			|| GetPdfunction()->at(node-1)->m_bound == Node::ENB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q15010_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q15010_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = -1.0*force.x - force.y + force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (vel.y*force.x - vel.z*force.x + vel.x*force.y - vel.z*force.y - vel.x*force.z - vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q15010_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), temp(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);

	temp =	pow(vel.y, 2) - pow(vel.z, 2); 

	moment_eq = local_rho * temp;
	
	return moment_eq;
}

Node::NodeValueType_t
D3q15010_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(1.0), M4(1.0), M5(-1.0), M6(-1.0),
		 		M7(0.0), M8(0.0), M9(0.0), M10(0.0), M11(0.0), M12(0.0), M13(0.0), M14(0.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q15010_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A10(1.0/tau);//related to kinematic viscosity

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A10* (mf - mf_eq) - (I - 0.5*A10) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q15010_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = 2.0 * (force.y * vel.y - force.z * vel.z);

	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q15011_incomp concrete class
//Properties:
//1. Derived from D3q15_011 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1509_incomp::D3q1509_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(-1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q1509_incomp::D3q1509_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(-1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1509_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+5);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1509_incomp::~D3q1509_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q1509_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1509_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1509_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1509_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1509_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1509_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1509_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1509_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1509_incomp::c_z()
{
	return cz;
}

void 
D3q1509_incomp::setNodeBType(cgsize_t  node)
{
	if(		 GetPdfunction()->at(node-1)->m_bound == Node::WEST
			|| GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			|| GetPdfunction()->at(node-1)->m_bound == Node::WB
			|| GetPdfunction()->at(node-1)->m_bound == Node::NB
			|| GetPdfunction()->at(node-1)->m_bound == Node::WN
			|| GetPdfunction()->at(node-1)->m_bound == Node::WNB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q1509_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1509_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = force.x - force.y + force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (-1.0*vel.y*force.x + vel.z*force.x - vel.x*force.y - vel.z*force.y + vel.x*force.z - vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q1509_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), temp(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	
	temp =	2.0 *	pow(vel.x, 2) - pow(vel.y, 2) - pow(vel.z, 2);

	moment_eq = local_rho * temp;
		    	   	
	return moment_eq;
}

Node::NodeValueType_t
D3q1509_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(2.0), M2(2.0), M3(-1.0), M4(-1.0), M5(-1.0), M6(-1.0),
		 		M7(0.0), M8(0.0), M9(0.0), M10(0.0), M11(0.0), M12(0.0), M13(0.0), M14(0.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1509_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A9(1.0/tau);//related to kinematic viscosity

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node); 
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A9 * (mf - mf_eq) - (I - 0.5*A9) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;	
}

Node::NodeValueType_t
D3q1509_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   =       (4 * force.x * vel.x - 2 * force.y * vel.y - 2 * force.z * vel.z );

	return sourceMoment;	
}
/*--------------------------------------------------------------------------------------------*/
//D3q15012_incomp concrete class
//Properties:
//1. Derived from D3q15_012 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q15012_incomp::D3q15012_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q15012_incomp::D3q15012_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q15012_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+5);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+3);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q15012_incomp::~D3q15012_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q15012_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q15012_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q15012_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q15012_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q15012_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q15012_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q15012_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q15012_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q15012_incomp::c_z()
{
	return cz;
}

void 
D3q15012_incomp::setNodeBType(cgsize_t  node)
{
	if(    GetPdfunction()->at(node-1)->m_bound == Node::EAST
			|| GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::TOP
			|| GetPdfunction()->at(node-1)->m_bound == Node::ET
			|| GetPdfunction()->at(node-1)->m_bound == Node::ST
			|| GetPdfunction()->at(node-1)->m_bound == Node::ES
			|| GetPdfunction()->at(node-1)->m_bound == Node::EST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q15012_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q15012_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = -1.0*force.x + force.y - force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (-1.0*vel.y*force.x + vel.z*force.x - vel.x*force.y - vel.z*force.y + vel.x*force.z - vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q15012_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);

	moment_eq = local_rho  * vel.y * vel.z;	

	return moment_eq;
}


Node::NodeValueType_t
D3q15012_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(1.0), M9(-1.0), M10(-1.0), M11(-1.0), M12(-1.0), M13(1.0), M14(1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q15012_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A12(1.0/tau);//related to kinematic viscosity

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A12 * (mf - mf_eq) - (I - 0.5*A12) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q15012_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = force.z * vel.y + force.y * vel.z;
	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q15013_incomp concrete class
//Properties:
//1. Derived from D3q15_013 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q15013_incomp::D3q15013_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(-1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q15013_incomp::D3q15013_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(-1.0), cz(-1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q15013_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+7);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+1);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q15013_incomp::~D3q15013_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q15013_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q15013_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q15013_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q15013_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q15013_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q15013_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q15013_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q15013_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q15013_incomp::c_z()
{
	return cz;
}

void 
D3q15013_incomp::setNodeBType(cgsize_t  node)
{
	if(    GetPdfunction()->at(node-1)->m_bound == Node::WEST
			|| GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::TOP
			|| GetPdfunction()->at(node-1)->m_bound == Node::WT
			|| GetPdfunction()->at(node-1)->m_bound == Node::NT
			|| GetPdfunction()->at(node-1)->m_bound == Node::WN
			|| GetPdfunction()->at(node-1)->m_bound == Node::WNT
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q15013_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q15013_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = force.x - force.y - force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (-1.0*vel.y*force.x - vel.z*force.x - vel.x*force.y + vel.z*force.y - vel.x*force.z + vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q15013_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);

	moment_eq = local_rho  * vel.x * vel.z;	

	return moment_eq;
}


Node::NodeValueType_t
D3q15013_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(0.0), M6(0.0),
		 		M7(1.0), M8(-1.0), M9(1.0), M10(-1.0), M11(-1.0), M12(1.0), M13(-1.0), M14(1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;	

	return moment;
}

Node::NodeValueType_t
D3q15013_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{	
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A13(1.0/tau);//related to kinematic viscosity

	mf           = compute_moment(allPdfs, node);
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node);
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node); 
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A13 * (mf - mf_eq) - (I - 0.5*A13) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q15013_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = force.z * vel.x + force.x * vel.z;
	return sourceMoment;
}
/*--------------------------------------------------------------------------------------------*/
//D3q15014_incomp concrete class
//Properties:
//1. Derived from D3q15_014 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1508_incomp::D3q1508_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(pCase);
}
D3q1508_incomp::D3q1508_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(1.0), cz(1.0)
{
	weight = 1.0/72.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1508_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+7);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1508_incomp::~D3q1508_incomp()
{

} 
//function that read fields

PdfBlock::NodeType_t* 
D3q1508_incomp::pdfunction() 
{
	return &pdf;
}
PdfBlock::NodeType_t* 
D3q1508_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}
PdfBlock::NodeType_t* 
D3q1508_incomp::pdfmoment() 
{
	return &moment;
}
 
PdfBlock::NodeType_t* 
D3q1508_incomp::pdfmoment_eq() 
{
	return &moment_eq;
}
CLbmCase*   
D3q1508_incomp::caseInfo()
{
	return m_pCase;
}

Node::NodeValueType_t 
D3q1508_incomp::weightFactor()
{
	return weight;
}

Node::NodeValueType_t 
D3q1508_incomp::c_x()
{
	return cx;
}

Node::NodeValueType_t 
D3q1508_incomp::c_y()
{
	return cy;
}

Node::NodeValueType_t 
D3q1508_incomp::c_z()
{
	return cz;
}


void 
D3q1508_incomp::setNodeBType(cgsize_t  node)
{
	if(    GetPdfunction()->at(node-1)->m_bound == Node::EAST
			|| GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			|| GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			|| GetPdfunction()->at(node-1)->m_bound == Node::EB
			|| GetPdfunction()->at(node-1)->m_bound == Node::SB
			|| GetPdfunction()->at(node-1)->m_bound == Node::ES
			|| GetPdfunction()->at(node-1)->m_bound == Node::ESB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
Node::NodeValueType_t
D3q1508_incomp::transformPDF(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t local_rho(0.0), transformedPdf(0.0);
	cgsize_t rho(4);
	local_rho      = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	
	transformedPdf = 0.5*computeSource(domainVariables, mat_index, local_rho, node);  
	
	return transformedPdf;
}
Node::NodeValueType_t
D3q1508_incomp::computeSource(Domain& domainVariables, cgsize_t mat_index, Node::NodeValueType_t localRho, cgsize_t node)
{
	Node::NodeValueType_t source(0.0), term1(0.0), term2(0.0), term3(0.0);
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);
	Force3D    force  = GetTotalForce(domainVariables, mat_index, node, localRho);

	term1 = -1.0*force.x + force.y + force.z ;
	term2 = 2.0 * (vel.x*force.x + vel.y*force.y + vel.z*force.z);
	term3 = 3.0 * (-1.0*vel.y*force.x - vel.z*force.x - vel.x*force.y + vel.z*force.y - vel.x*force.z + vel.y*force.z);

	source = 3.0 * GetWk() * ( term1 + term2 + term3);
	
	return source;
}
Node::NodeValueType_t
D3q1508_incomp::compute_moment_eq(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t moment_eq(0.0), local_rho(0.0);
	cgsize_t rho(4);
	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Velocity3D vel    = GetOverallVel(domainVariables, mat_index, node);
	//Velocity3D vel    = GetComponentVel(domainVariables, mat_index, node);

	moment_eq = (-7.0/3.0) * local_rho * vel.z;

	return moment_eq;
}


Node::NodeValueType_t
D3q1508_incomp::compute_moment(PdfDomain* allPdfs, cgsize_t node)
{
	Node::NodeValueType_t moment(0.0);

	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	Node::NodeValueType_t 	M0(0.0), M1(0.0), M2(0.0), M3(0.0), M4(0.0), M5(-4.0), M6(4.0),
		 		M7(1.0), M8(1.0), M9(1.0), M10(1.0), M11(-1.0), M12(-1.0), M13(-1.0), M14(-1.0);
	
	moment = (M0 * allPdfs->GetLatticePdf()->at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) + 
		 (M1 * allPdfs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M2 * allPdfs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M3 * allPdfs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M4 * allPdfs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M5 * allPdfs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M6 * allPdfs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M7 * allPdfs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M8 * allPdfs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M9 * allPdfs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +

		 (M10 * allPdfs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M11 * allPdfs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M12 * allPdfs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M13 * allPdfs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) +
		 (M14 * allPdfs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal) ;

	return moment;
}

Node::NodeValueType_t
D3q1508_incomp::compute_collision(Domain& domainVariables, PdfDomain* allPdfs, cgsize_t mat_index,  Node::NodeValueType_t tau, cgsize_t node)
{
	Node::NodeValueType_t rhs_coll(0.0), mf(0.0), mf_eq(0.0), sourceMoment(0.0);
	Node::NodeValueType_t I(1.0), A8(1.0);//free parameters that does not relate to the hydrodynamics (artibrary)

	mf           = compute_moment(allPdfs, node);	
	mf_eq        = compute_moment_eq(domainVariables, mat_index, node); 
	sourceMoment = computeSourceMoment(domainVariables, mat_index, node);
	/*if(node == 620000)
	{
		std::cout<<" mf   = "<<mf<<"\n";
		std::cout<<" mfeq = "<<mf_eq<<"\n";
		std::cout<<" mS   = "<<sourceMoment<<"\n";
		std::cin.get();
	}*/
	rhs_coll = A8 * (mf - mf_eq) - (I - 0.5*A8) * sourceMoment;
	/*if(node == 620000)
	{
		std::cout<<" rhs_coll   = "<<rhs_coll<<"\n";
		std::cin.get();
	}*/
	return rhs_coll;
}

Node::NodeValueType_t
D3q1508_incomp::computeSourceMoment(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	Node::NodeValueType_t sourceMoment(0.0), local_rho(0.0);
	cgsize_t rho(4);

	local_rho = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
	Force3D force = GetTotalForce(domainVariables, mat_index, node, local_rho);
	
	sourceMoment   = (-7.0/3.0) * force.z;

	return sourceMoment;
}
//============================================================================================//

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_00 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_00::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1500_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_00::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1500_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_01 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_01::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1501_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_01::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1501_incomp(pCase,nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_02 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_02::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1502_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_02::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1502_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_03 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_03::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1503_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_03::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1503_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_04 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_04::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1504_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_04::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1504_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_05 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_05::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1505_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_05::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1505_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_06 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_06::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1506_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_06::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1506_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_07 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_07::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1507_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_07::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1507_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_08 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_08::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1508_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_08::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1508_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_09 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_09::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1509_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_09::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1509_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_10 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_10::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q15010_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_10::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q15010_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_11 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_11::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q15011_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_11::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q15011_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_12 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_12::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q15012_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_12::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q15012_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_13 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_13::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q15013_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_13::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q15013_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q15_14 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q15_14::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q15014_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q15_14::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q15014_incomp(pCase,nvertex);
	return pblock;
}


//============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//D3q1900_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1900_incomp::D3q1900_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0), cz(0.0)
{
	weight = 1.0/3.0;
	setBlock(pCase);
}

D3q1900_incomp::D3q1900_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(0.0), cz(0.0)
{
	weight = 1.0/3.0;
	setBlock(nvertex);
}
void
D3q1900_incomp::setHexa_8conn(cgsize_t i){
	//NULL BODY
}
//copy ctor: using default copy ctor

//dtor
D3q1900_incomp::~D3q1900_incomp()
{

}
//functions that reads PdfBlock fields

inline PdfBlock::NodeType_t* 
D3q1900_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1900_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1900_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1900_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1900_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1900_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1900_incomp::c_z()
{
	return cz;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1901_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1901_incomp::D3q1901_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(0.0), cz(0.0)
{
	//set block info
	weight = 1.0/18.0;
	setBlock(pCase);
}
D3q1901_incomp::D3q1901_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(0.0), cz(0.0)
{
	//set block info
	weight = 1.0/18.0;
	setBlock(nvertex);
}
void 
D3q1901_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeF, nodeB, tonodeF, tonodeB;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeF   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeF = *(*(m_pCase->m_element.m_conn + i) + j+1);
		nodeB   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonodeB = *(*(m_pCase->m_element.m_conn + i) + j+2);
			if(GetPdfunction()->at(tonodeF-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeF-1)->m_connectTo = nodeF;
				GetPdfunction()->at(tonodeF-1)->m_linkUpdated   = true;
			}
			if (GetPdfunction()->at(tonodeB-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeB-1)->m_connectTo = nodeB;
				GetPdfunction()->at(tonodeB-1)->m_linkUpdated   = true;
			}
	}
}

//copy ctor: using default copy ctor
//dtor
D3q1901_incomp::~D3q1901_incomp()
{

}

//function that read fields

inline PdfBlock::NodeType_t* 
D3q1901_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1901_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1901_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1901_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1901_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1901_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1901_incomp::c_z()
{
	return cz;
}

void 
D3q1901_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::WEST)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WSB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WNB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WST ||
						GetPdfunction()->at(node-1)->m_bound == Node::WNT ||
						GetPdfunction()->at(node-1)->m_bound == Node::WS  ||
						GetPdfunction()->at(node-1)->m_bound == Node::WN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::WT  ||
						GetPdfunction()->at(node-1)->m_bound == Node::WB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1902_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1902_incomp::D3q1902_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(0.0), cz(0.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(pCase);
}
D3q1902_incomp::D3q1902_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(0.0), cz(0.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(nvertex);
}
void 
D3q1902_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeF, nodeB, tonodeF, tonodeB;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeF   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeF = *(*(m_pCase->m_element.m_conn + i) + j+3);
		nodeB   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeB = *(*(m_pCase->m_element.m_conn + i) + j);
			if(GetPdfunction()->at(tonodeF-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeF-1)->m_connectTo   = nodeF;
				GetPdfunction()->at(tonodeF-1)->m_linkUpdated = true;
			}
			if(GetPdfunction()->at(tonodeB-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeB-1)->m_connectTo   = nodeB;
				GetPdfunction()->at(tonodeB-1)->m_linkUpdated = true;
			}
	}
}

//copy ctor: using default copy ctor
//dtor
D3q1902_incomp::~D3q1902_incomp()
{

}

//function that read fields

inline PdfBlock::NodeType_t* 
D3q1902_incomp::pdfunction()
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1902_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1902_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1902_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1902_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1902_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1902_incomp::c_z()
{
	return cz;
}

void 
D3q1902_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::EAST)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::ESB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENB ||
						GetPdfunction()->at(node-1)->m_bound == Node::EST ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENT ||
						GetPdfunction()->at(node-1)->m_bound == Node::ES  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ET  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1903_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1903_incomp::D3q1903_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(1.0), cz(0.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(pCase);
}
D3q1903_incomp::D3q1903_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(1.0), cz(0.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(nvertex);
}
void 
D3q1903_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeO, nodeI, tonodeO, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+3);
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+2);
			if(GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeO-1)->m_connectTo   = nodeO;
				GetPdfunction()->at(tonodeO-1)->m_linkUpdated = true;
			}
			if(GetPdfunction()->at(tonodeI-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeI-1)->m_connectTo   = nodeI;
				GetPdfunction()->at(tonodeI-1)->m_linkUpdated = true;
			}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1903_incomp::~D3q1903_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1903_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1903_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1903_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1903_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1903_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1903_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1903_incomp::c_z()
{
	return cz;
}

void 
D3q1903_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::SOUTH)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WSB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ESB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WST ||
						GetPdfunction()->at(node-1)->m_bound == Node::EST ||
						GetPdfunction()->at(node-1)->m_bound == Node::WS  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ES  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ST  ||
						GetPdfunction()->at(node-1)->m_bound == Node::SB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1904_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1904_incomp::D3q1904_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0), cy(-1.0), cz(0.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(pCase);
}
D3q1904_incomp::D3q1904_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0), cy(-1.0), cz(0.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(nvertex);
}
void 
D3q1904_incomp::setHexa_8conn(cgsize_t i){
	cgsize_t nodeO, nodeI, tonodeO, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += 4)
	{
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j);
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+1);
			if(GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonodeO-1)->m_connectTo   = nodeO;
				GetPdfunction()->at(tonodeO-1)->m_linkUpdated = true;
			}
			if(GetPdfunction()->at(tonodeI-1)->m_linkUpdated == false)
			{
				GetPdfunction()->at(tonodeI-1)->m_connectTo   = nodeI;
				GetPdfunction()->at(tonodeI-1)->m_linkUpdated = true;
			}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1904_incomp::~D3q1904_incomp()
{

} 

//function that read fields

inline PdfBlock::NodeType_t* 
D3q1904_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1904_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1904_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1904_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1904_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1904_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1904_incomp::c_z()
{
	return cz;
}

void 
D3q1904_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::NORTH)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WNB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WNT ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENT ||
						GetPdfunction()->at(node-1)->m_bound == Node::WN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EN  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NT  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1905_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1905_incomp::D3q1905_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0), cz(1.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(pCase);
}
D3q1905_incomp::D3q1905_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(0.0), cz(1.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1905_incomp::setHexa_8conn(cgsize_t i){
	const cgsize_t hexa = 8;
	const cgsize_t quad = 4;
	cgsize_t temp(0);
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		temp = hexa -1;
		for(cgsize_t k=temp; k >= quad; k--)
		{
			node    = *(*(m_pCase->m_element.m_conn + i) + j+(k-quad));
			tonode  = *(*(m_pCase->m_element.m_conn + i) + j+k);
		
			if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonode-1)->m_connectTo     = node;
				GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
			}
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1905_incomp::~D3q1905_incomp()
{

} 

//function that read fields

inline PdfBlock::NodeType_t* 
D3q1905_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1905_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1905_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1905_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1905_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1905_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1905_incomp::c_z()
{
	return cz;
}

void 
D3q1905_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WNB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WSB ||
						GetPdfunction()->at(node-1)->m_bound == Node::ESB ||
						GetPdfunction()->at(node-1)->m_bound == Node::WB  ||
						GetPdfunction()->at(node-1)->m_bound == Node::EB  ||
						GetPdfunction()->at(node-1)->m_bound == Node::SB  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NB
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1906_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1906_incomp::D3q1906_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(pCase);
}
D3q1906_incomp::D3q1906_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/18.0;;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1906_incomp::setHexa_8conn(cgsize_t i){
	const cgsize_t hexa = 8;
	const cgsize_t quad = 4;
	cgsize_t temp(0);
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		temp = hexa -1;
		for(cgsize_t k=temp; k >= quad; k--)
		{
			node    = *(*(m_pCase->m_element.m_conn + i) + j+k);
			tonode  = *(*(m_pCase->m_element.m_conn + i) + j+(k-quad));
		
			if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
			{
				GetPdfunction()->at(tonode-1)->m_connectTo     = node;
				GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
			}
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1906_incomp::~D3q1906_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1906_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1906_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1906_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1906_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1906_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1906_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1906_incomp::c_z()
{
	return cz;
}

void 
D3q1906_incomp::setNodeBType(cgsize_t  node)
{
	if(GetPdfunction()->at(node-1)->m_bound == Node::TOP)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else if(GetPdfunction()->at(node-1)->m_periodicNeighborNum != 0
					&&
					(
						GetPdfunction()->at(node-1)->m_bound == Node::WNT ||
						GetPdfunction()->at(node-1)->m_bound == Node::ENT ||
						GetPdfunction()->at(node-1)->m_bound == Node::WST ||
						GetPdfunction()->at(node-1)->m_bound == Node::EST ||
						GetPdfunction()->at(node-1)->m_bound == Node::WT  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ET  ||
						GetPdfunction()->at(node-1)->m_bound == Node::ST  ||
						GetPdfunction()->at(node-1)->m_bound == Node::NT
					)
		)
		GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;

	else
		GetPdfunction()->at(node-1)->m_boundType = Node::NONE;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1907_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1907_incomp::D3q1907_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1907_incomp::D3q1907_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1907_incomp::setHexa_8conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+2);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1907_incomp::~D3q1907_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1907_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1907_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1907_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1907_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1907_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1907_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1907_incomp::c_z()
{
	return cz;
}

void 
D3q1907_incomp::setNodeBType(cgsize_t  node)
{
	if(
					GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
			||	GetPdfunction()->at(node-1)->m_bound == Node::WB
			||	GetPdfunction()->at(node-1)->m_bound == Node::SB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WT
			||	GetPdfunction()->at(node-1)->m_bound == Node::ST
			||	GetPdfunction()->at(node-1)->m_bound == Node::WSB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1908_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1908_incomp::D3q1908_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(-1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1908_incomp::D3q1908_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(-1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1908_incomp::setHexa_8conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+1);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1908_incomp::~D3q1908_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1908_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1908_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1908_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1908_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1908_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1908_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1908_incomp::c_z()
{
	return cz;
}

void 
D3q1908_incomp::setNodeBType(cgsize_t  node)
{
	if(  
					GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::WN
			||	GetPdfunction()->at(node-1)->m_bound == Node::WB
			||	GetPdfunction()->at(node-1)->m_bound == Node::NB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WT
			||	GetPdfunction()->at(node-1)->m_bound == Node::NT
			||	GetPdfunction()->at(node-1)->m_bound == Node::WNB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WNT
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1909_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1909_incomp::D3q1909_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1909_incomp::D3q1909_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1909_incomp::setHexa_8conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j+3);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1909_incomp::~D3q1909_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1909_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1909_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1909_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1909_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1909_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1909_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1909_incomp::c_z()
{
	return cz;
}

void 
D3q1909_incomp::setNodeBType(cgsize_t  node)
{
	if( 
					GetPdfunction()->at(node-1)->m_bound == Node::EAST
			||	GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::ES
			||	GetPdfunction()->at(node-1)->m_bound == Node::EB
			||	GetPdfunction()->at(node-1)->m_bound == Node::SB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ET
			||	GetPdfunction()->at(node-1)->m_bound == Node::ST
			||	GetPdfunction()->at(node-1)->m_bound == Node::ESB
			||	GetPdfunction()->at(node-1)->m_bound == Node::EST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1910_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1910_incomp::D3q1910_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(-1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1910_incomp::D3q1910_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(-1.0), cz(0.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1910_incomp::setHexa_8conn(cgsize_t i){
	const int quad = 4;
	cgsize_t node, tonode;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += quad)
	{
		node   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonode = *(*(m_pCase->m_element.m_conn + i) + j);
	
		if( GetPdfunction()->at(tonode-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonode-1)->m_connectTo     = node;
			GetPdfunction()->at(tonode-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1910_incomp::~D3q1910_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1910_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1910_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1910_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1910_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1910_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1910_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1910_incomp::c_z()
{
	return cz;
}

void 
D3q1910_incomp::setNodeBType(cgsize_t  node)
{
	if(	
					GetPdfunction()->at(node-1)->m_bound == Node::EAST
			||	GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::EN
			||	GetPdfunction()->at(node-1)->m_bound == Node::NB
			||	GetPdfunction()->at(node-1)->m_bound == Node::EB
			||	GetPdfunction()->at(node-1)->m_bound == Node::NT
			||	GetPdfunction()->at(node-1)->m_bound == Node::ET
			||	GetPdfunction()->at(node-1)->m_bound == Node::ENB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ENT
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1911_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1911_incomp::D3q1911_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(0.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1911_incomp::D3q1911_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(0.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1911_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeO, tonodeO, nodeI, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);
	
	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+5);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+6);
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1911_incomp::~D3q1911_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1911_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1911_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1911_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1911_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1911_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1911_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1911_incomp::c_z()
{
	return cz;
}

void 
D3q1911_incomp::setNodeBType(cgsize_t  node)
{
	if(	
					GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			||	GetPdfunction()->at(node-1)->m_bound == Node::WB
			||	GetPdfunction()->at(node-1)->m_bound == Node::NB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WN
			||	GetPdfunction()->at(node-1)->m_bound == Node::SB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
			||	GetPdfunction()->at(node-1)->m_bound == Node::WNB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WSB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1912_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1912_incomp::D3q1912_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(0.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1912_incomp::D3q1912_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(0.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1912_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeO, tonodeO, nodeI, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+4);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+7);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1912_incomp::~D3q1912_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1912_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1912_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1912_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1912_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1912_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1912_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1912_incomp::c_z()
{
	return cz;
}

void 
D3q1912_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::EAST
			||	GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			||	GetPdfunction()->at(node-1)->m_bound == Node::EB
			||	GetPdfunction()->at(node-1)->m_bound == Node::NB
			||	GetPdfunction()->at(node-1)->m_bound == Node::EN
			||	GetPdfunction()->at(node-1)->m_bound == Node::SB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ES
			||	GetPdfunction()->at(node-1)->m_bound == Node::ENB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ESB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1913_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1913_incomp::D3q1913_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(1.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1913_incomp::D3q1913_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(1.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1913_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeO, tonodeO, nodeI, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+4);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+1);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+7);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+2);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1913_incomp::~D3q1913_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1913_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1913_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1913_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1913_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1913_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1913_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1913_incomp::c_z()
{
	return cz;
}

void 
D3q1913_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::WEST
			||	GetPdfunction()->at(node-1)->m_bound == Node::TOP
			||	GetPdfunction()->at(node-1)->m_bound == Node::WT
			||	GetPdfunction()->at(node-1)->m_bound == Node::NT
			||	GetPdfunction()->at(node-1)->m_bound == Node::WN
			||	GetPdfunction()->at(node-1)->m_bound == Node::ST
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
			||	GetPdfunction()->at(node-1)->m_bound == Node::WNT
			||	GetPdfunction()->at(node-1)->m_bound == Node::WST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
/*--------------------------------------------------------------------------------------------*/
//D3q1914_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1914_incomp::D3q1914_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(-1.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1914_incomp::D3q1914_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(-1.0), cy(0.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1914_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeO, tonodeO, nodeI, tonodeI;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+5);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+6);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+3);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1914_incomp::~D3q1914_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1914_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1914_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1914_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1914_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1914_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1914_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1914_incomp::c_z()
{
	return cz;
}

void 
D3q1914_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::EAST
			||	GetPdfunction()->at(node-1)->m_bound == Node::TOP
			||	GetPdfunction()->at(node-1)->m_bound == Node::ET
			||	GetPdfunction()->at(node-1)->m_bound == Node::NT
			||	GetPdfunction()->at(node-1)->m_bound == Node::EN
			||	GetPdfunction()->at(node-1)->m_bound == Node::ST
			||	GetPdfunction()->at(node-1)->m_bound == Node::ES
			||	GetPdfunction()->at(node-1)->m_bound == Node::ENT
			||	GetPdfunction()->at(node-1)->m_bound == Node::EST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1915_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1915_incomp::D3q1915_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(1.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1915_incomp::D3q1915_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(1.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1915_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeI, tonodeI, nodeO, tonodeO;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+7);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+1);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+6);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1915_incomp::~D3q1915_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1915_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1915_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1915_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1915_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1915_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1915_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1915_incomp::c_z()
{
	return cz;
}

void 
D3q1915_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			||	GetPdfunction()->at(node-1)->m_bound == Node::SB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
			||	GetPdfunction()->at(node-1)->m_bound == Node::EB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ES
			||	GetPdfunction()->at(node-1)->m_bound == Node::WSB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ESB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1916_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1916_incomp::D3q1916_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(1.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1916_incomp::D3q1916_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(1.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1916_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeI, tonodeI, nodeO, tonodeO;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+4);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+3);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+5);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+2);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1916_incomp::~D3q1916_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1916_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1916_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1916_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1916_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1916_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1916_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1916_incomp::c_z()
{
	return cz;
}

void 
D3q1916_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::SOUTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::TOP
			||	GetPdfunction()->at(node-1)->m_bound == Node::ST
			||	GetPdfunction()->at(node-1)->m_bound == Node::WS
			||	GetPdfunction()->at(node-1)->m_bound == Node::WT
			||	GetPdfunction()->at(node-1)->m_bound == Node::ES
			||	GetPdfunction()->at(node-1)->m_bound == Node::ET
			||	GetPdfunction()->at(node-1)->m_bound == Node::WST
			||	GetPdfunction()->at(node-1)->m_bound == Node::EST
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1917_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1917_incomp::D3q1917_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(-1.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1917_incomp::D3q1917_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(-1.0), cz(1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1917_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeI, tonodeI, nodeO, tonodeO;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+2);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+5);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+3);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j+4);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1917_incomp::~D3q1917_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1917_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1917_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1917_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1917_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1917_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1917_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1917_incomp::c_z()
{
	return cz;
}

void 
D3q1917_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::BOTTOM
			||	GetPdfunction()->at(node-1)->m_bound == Node::NB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WN
			||	GetPdfunction()->at(node-1)->m_bound == Node::WB
			||	GetPdfunction()->at(node-1)->m_bound == Node::EN
			||	GetPdfunction()->at(node-1)->m_bound == Node::EB
			||	GetPdfunction()->at(node-1)->m_bound == Node::ENB
			||	GetPdfunction()->at(node-1)->m_bound == Node::WNB
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}

/*--------------------------------------------------------------------------------------------*/
//D3q1918_incomp concrete class
//Properties:
//1. Derived from Incomp_D3Q19 Implementation level class
/*--------------------------------------------------------------------------------------------*/
//non-default ctor
D3q1918_incomp::D3q1918_incomp (CLbmCase* pCase): m_pCase(pCase)
	, cx(0.0), cy(-1.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(pCase);
}
D3q1918_incomp::D3q1918_incomp (CLbmCase* pCase, cgsize_t nvertex): m_pCase(pCase),
cx(0.0), cy(-1.0), cz(-1.0)
{
	weight = 1.0/36.0;
	//set block info
	setBlock(nvertex);
}
void //remove later does not apply to this model
D3q1918_incomp::setHexa_8conn(cgsize_t i){
	const int hexa = 8;
	cgsize_t nodeI, tonodeI, nodeO, tonodeO;
	cgsize_t  dim  = *(m_pCase->m_element.m_elementDataSize + i);

	for(cgsize_t j=0; j < dim; j += hexa)
	{
		nodeI   = *(*(m_pCase->m_element.m_conn + i) + j+6);
		tonodeI = *(*(m_pCase->m_element.m_conn + i) + j+1);
		nodeO   = *(*(m_pCase->m_element.m_conn + i) + j+7);
		tonodeO = *(*(m_pCase->m_element.m_conn + i) + j);
	
		if( GetPdfunction()->at(tonodeI-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeI-1)->m_connectTo     = nodeI;
			GetPdfunction()->at(tonodeI-1)->m_linkUpdated   = true;
		}
		if( GetPdfunction()->at(tonodeO-1)->m_linkUpdated   == false)
		{
			GetPdfunction()->at(tonodeO-1)->m_connectTo     = nodeO;
			GetPdfunction()->at(tonodeO-1)->m_linkUpdated   = true;
		}
	}
}
//copy ctor: using default copy ctor
//dtor
D3q1918_incomp::~D3q1918_incomp()
{

} 
//function that read fields

inline PdfBlock::NodeType_t* 
D3q1918_incomp::pdfunction() 
{
	return &pdf;
}
inline PdfBlock::NodeType_t* 
D3q1918_incomp::pdfunction_eq() 
{
	return &pdf_eq;
}

inline CLbmCase*   
D3q1918_incomp::caseInfo()
{
	return m_pCase;
}

inline Node::NodeValueType_t 
D3q1918_incomp::weightFactor()
{
	return weight;
}

inline Node::NodeValueType_t 
D3q1918_incomp::c_x()
{
	return cx;
}

inline Node::NodeValueType_t 
D3q1918_incomp::c_y()
{
	return cy;
}

inline Node::NodeValueType_t 
D3q1918_incomp::c_z()
{
	return cz;
}

void 
D3q1918_incomp::setNodeBType(cgsize_t  node)
{
	if( 		GetPdfunction()->at(node-1)->m_bound == Node::NORTH
			||	GetPdfunction()->at(node-1)->m_bound == Node::TOP
			||	GetPdfunction()->at(node-1)->m_bound == Node::NT
			||	GetPdfunction()->at(node-1)->m_bound == Node::WN
			||	GetPdfunction()->at(node-1)->m_bound == Node::WT
			||	GetPdfunction()->at(node-1)->m_bound == Node::EN
			||	GetPdfunction()->at(node-1)->m_bound == Node::ET
			||	GetPdfunction()->at(node-1)->m_bound == Node::WNT
			||	GetPdfunction()->at(node-1)->m_bound == Node::ENT
			)
			GetPdfunction()->at(node-1)->m_boundType = Node::INWARD;
		else
			GetPdfunction()->at(node-1)->m_boundType = Node::BURIED;
}
//============================================================================================//
//============================================================================================//

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_00 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_00::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1900_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_00::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1900_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_01 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_01::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1901_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_01::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1901_incomp(pCase,nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_02 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_02::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1902_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_02::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1902_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_03 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_03::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1903_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_03::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1903_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_04 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_04::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1904_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_04::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1904_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_05 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_05::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1905_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_05::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1905_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_06 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_06::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1906_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_06::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1906_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_07 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_07::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1907_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_07::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1907_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_08 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_08::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1908_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_08::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1908_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_09 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_09::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1909_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_09::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1909_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_10 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_10::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1910_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_10::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1910_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_11 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_11::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1911_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_11::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1911_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_12 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_12::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1912_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_12::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1912_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_13 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_13::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1913_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_13::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1913_incomp(pCase, nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_14 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_14::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1914_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_14::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1914_incomp(pCase,nvertex);
	return pblock;
}

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_15 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_15::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1915_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_15::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1915_incomp(pCase,nvertex);
	return pblock;
}

/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_16 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_16::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1916_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_16::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1916_incomp(pCase,nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_17 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_17::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1917_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_17::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1917_incomp(pCase,nvertex);
	return pblock;
}
/*--------------------------------------------------------------------------------------------*/
//MakeD3Q19_18 derived concrete class
//Properties:
//1. implementation/leaf level class; 
/*--------------------------------------------------------------------------------------------*/
PdfBlock*  
MakeD3Q19_18::incomp(CLbmCase* pCase)
{
	PdfBlock* pblock = new D3q1918_incomp(pCase);
	return pblock;
}
PdfBlock*  
MakeD3Q19_18::incomp(CLbmCase* pCase, cgsize_t nvertex)
{
	PdfBlock* pblock = new D3q1918_incomp(pCase,nvertex);
	return pblock;
}
