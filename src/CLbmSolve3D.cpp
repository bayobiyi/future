
#include "CLbmSolve3D.h"
#include <cmath>
#include <iostream>

//Domain_D3QN
void 
Domain_D3QN::convergence(Domain& domainVariables, CLbmCase* pCase)
{
	cgsize_t xvel(0), yvel(1), zvel(2), density(4), uoverall_x(5), uoverall_y(6), uoverall_z(7);
	clock_t time;
	time = clock();
	
	m_convergenceData.x_residual = domainVariables.GetDomainVariables()->at(uoverall_x)->residual();
	m_convergenceData.y_residual = domainVariables.GetDomainVariables()->at(uoverall_y)->residual();
	m_convergenceData.z_residual = domainVariables.GetDomainVariables()->at(uoverall_z)->residual();
	m_convergenceData.density_residual = domainVariables.GetDomainVariables()->at(density)->residual();
	
	time = (clock() - time);/*CLOCKS_PER_SEC;*/
	///////////////////////////////////////////////////////////////////////////////////////////////////////////
	//std::cout<<"time to compute convergence = "<< time <<"\n";
	//std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}
void 
Domain_D3QN::computeCompositeVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	Node::NodeValueType_t ux_sum(0.0), uy_sum(0.0), uz_sum(0.0), temp(0.0), tau_primary(0.0), tau_secondary(0.0);
	cgsize_t xvel(0), yvel(1), zvel(2), rho(4), uprime_x(8), uprime_y(9), uprime_z(10);
	cgsize_t primary(0), secondary(1);

	//obtain tau0 and tau1 first
	std::map <cgsize_t,CLbmElement::MaterialRecord_t>::const_iterator iter = pCase->m_element.m_sectionsMaterial.find(primary);
	tau_primary     = (3.0 * iter->second.second.k_viscosity + 0.5);
	
	iter = pCase->m_element.m_sectionsMaterial.find(secondary);
	tau_secondary     = (3.0 * iter->second.second.k_viscosity + 0.5);

	if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_isSolid == false)
	{
		//denominator of eqn 95 first
					

		ux_sum = (domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal / tau_primary )
			+      (domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal / tau_secondary );

		uy_sum = (domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal / tau_primary )
			+      (domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal / tau_secondary );

		uz_sum = (domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal / tau_primary )
			+      (domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal / tau_secondary );

		if((domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal) != 0)
		{
			temp = (domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal  / tau_primary) + 
							(domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal  / tau_secondary) ;

			domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = ux_sum / temp;
			domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = ux_sum / temp;
			domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = uy_sum / temp;
			domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = uy_sum / temp;
			domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = uz_sum / temp;
			domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = uz_sum / temp;
		}
		else
		{
			//set values to zero
			domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
			domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 0.0;
			domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
			domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 0.0;
			domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
			domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 0.0;
		}

	}
	else;
	/*update m_sol (try to update values of uprime_x... in m_sol array in preparation for export)
	pCase->m_sol.m_sol.at(uprime_x)->at(primary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_x)->at(secondary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_y)->at(primary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_y)->at(secondary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_z)->at(primary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_z)->at(secondary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;*/


	/*if(node == 8000)
	{
		std::cout<<" uprime_x: "<<uprime_x+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" uprime_y: "<<uprime_y+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" uprime_z: "<<uprime_z+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";

		std::cout<<" uprime_x: "<<uprime_x+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" uprime_y: "<<uprime_y+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" uprime_z: "<<uprime_z+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal<<"\n";
	}*/
}

void 
Domain_D3QN::computeOverallVars(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	cgsize_t primary(0), secondary(1), xvel(0), yvel(1), zvel(2), rho(4), uoverall_x(5), uoverall_y(6), uoverall_z(7);
	cgsize_t fx(1), fy(2), fz(3), fsurf_x(4), fsurf_y(5), fsurf_z(6);
	Node::NodeValueType_t c(0.5);
	Node::NodeValueType_t overall_density(0.0);
	
	//reset node vlaue to zero first
	domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;


	domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal += 
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;/* +
			c * domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;*/

	domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal += 
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(primary)->at(node-1)->m_nodeVal; /*+
			c * domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;*/

	domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal += 
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(primary)->at(node-1)->m_nodeVal; /*+
			c * domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;*/

		
	domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal += 
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal; /*+
			c * domainVariables.GetDomainTempVariables()->at(fsurf_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;*/

	domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal += 
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal; /*+
			c * domainVariables.GetDomainTempVariables()->at(fsurf_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;*/

	domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal += 
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal; /*+
			c * domainVariables.GetDomainTempVariables()->at(fsurf_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;*/


	overall_density = (domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal);

	if(overall_density != 0)
	{
		domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal /= overall_density;
		domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal /= overall_density;
		domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal /= overall_density;
	}
	else;
	
	/*update m_sol (try to update values of overall_x... in m_sol array in preparation for export)
	pCase->m_sol.m_sol.at(uoverall_x)->at(primary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	
	pCase->m_sol.m_sol.at(uoverall_y)->at(primary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	
	pCase->m_sol.m_sol.at(uoverall_z)->at(primary)->at(node-1) = 
		domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;*/

}
void 
Domain_D3QN::computeOverallVarsSgl(CLbmCase* pCase, Domain& domainVariables, cgsize_t node)
{
	cgsize_t primary(0), xvel(0), yvel(1), zvel(2), rho(4), uoverall_x(5), uoverall_y(6), uoverall_z(7);
	cgsize_t fx(1), fy(2), fz(3);
	Node::NodeValueType_t c(0.5);
	Node::NodeValueType_t overall_density(0.0);
	
	//reset node vlaue to zero first
	domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;

	/*if(node == 620000)
	{
		std::cout<<" u_sum = "<<domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" v_sum = "<<domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" w_sum = "<<domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";

		std::cout<<" fx = "<<domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" fy = "<<domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" fz = "<<domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
	}*/


	domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fx)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

	domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fy)->GetVariable()->at(primary)->at(node-1)->m_nodeVal; 

	domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal = 
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(primary)->at(node-1)->m_nodeVal + 
			c * domainVariables.GetDomainTempVariables()->at(fz)->GetVariable()->at(primary)->at(node-1)->m_nodeVal; 


	overall_density = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;

	if(overall_density != 0)
	{
		domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal /= overall_density;
		domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal /= overall_density;
		domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal /= overall_density;
	}
	else;
	/*if(node == 620000)
	{
		std::cout<<" u_overall = "<<domainVariables.GetDomainVariables()->at(uoverall_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" v_overall = "<<domainVariables.GetDomainVariables()->at(uoverall_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" w_overall = "<<domainVariables.GetDomainVariables()->at(uoverall_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal<<"\n";
	}*/
}
void 
Domain_D3QN::computeIndividualVars(CLbmCase* pCase, Domain& domainVariables, MaterialType_t::size_type mat_index, cgsize_t node)
{
	cgsize_t xvel(0), yvel(1), zvel(2), velmag(3), rho(4);

	if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_isSolid == false)
	{
		/*store density for export///////////////////////////////WHY? in case I need to use export to initialize continuing iteration
		pCase->m_sol.m_sol.at(rho)->at(mat_index)->at(node-1) = 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;*/

		//check division by zero
		if(domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal != 0)
		{
			//update x-velocity	
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    /= 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
				

			//pCase->m_sol.m_sol.at(xvel)->at(mat_index)->at(node-1) = 
				//domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
						

			//Update y-velocity	
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    /= 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;


			//pCase->m_sol.m_sol.at(yvel)->at(mat_index)->at(node-1) = 
				//domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
						
			//Update z-velocity	
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    /= 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

			//pCase->m_sol.m_sol.at(zvel)->at(mat_index)->at(node-1) = 
				//domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
						
			//update velocity magnitude
			domainVariables.GetDomainVariables()->at(velmag)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 
				sqrt(  pow(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal , 2) + 
					pow(domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal , 2) +
					pow(domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal , 2));
		}
		else
		{
			domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;
			domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;
			domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;
		}
	}
	else;
	/*if(node == 620000)
	{
		std::cout<<" u_comp = "<<domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" v_comp = "<<domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" w_comp = "<<domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
	}*/
	/*if(node == 8000)
	{
		std::cout<<" xvelocity: "<<xvel+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" yvelocity: "<<yvel+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
		std::cout<<" zvelocity: "<<zvel+1<<" for node = "<<node<<" (before) uploading is "<< domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	}*/
}

void 
Domain_D3QN::setNewToOld(Domain& domainVariables, cgsize_t node)
{
	cgsize_t primary(0), secondary(1), uprime_x(8), uprime_y(9), uprime_z(10);

	domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal   = 
				domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal ;
	domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal =
				domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
	domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal   = 
				domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
	domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_oldNodeVal   = 
				domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
}
void 
Domain_D3QN::setToZero(Domain& domainVariables, CLbmCase* pCase, cgsize_t node)
{
	cgsize_t primary(0), secondary(1), uprime_x(8), uprime_y(9), uprime_z(10);
	//first store data for export
	pCase->m_sol.m_sol.at(uprime_x)->at(primary)->at(node-1) = 
		domainVariables.GetDomainTempVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_x)->at(secondary)->at(node-1) = 
		domainVariables.GetDomainTempVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_y)->at(primary)->at(node-1) = 
		domainVariables.GetDomainTempVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_y)->at(secondary)->at(node-1) = 
		domainVariables.GetDomainTempVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_z)->at(primary)->at(node-1) = 
		domainVariables.GetDomainTempVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal;
	pCase->m_sol.m_sol.at(uprime_z)->at(secondary)->at(node-1) = 
		domainVariables.GetDomainTempVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal;

	//Now reset values to zero
	domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uprime_x)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 0.0;
	domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uprime_y)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 0.0;
	domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(primary)->at(node-1)->m_nodeVal   = 0.0;
	domainVariables.GetDomainVariables()->at(uprime_z)->GetVariable()->at(secondary)->at(node-1)->m_nodeVal = 0.0;
}
Convergence*   
Domain_D3QN::GetConvergence() 
{
	return &m_convergenceData;
}

inline void 
Domain_D3QN::SetNumThreads(cgsize_t  num)
{
	m_numThreads = num;//do some checking later
}

inline cgsize_t 
Domain_D3QN::GetNumThreads() const
{
	return m_numThreads;
}
Domain_D3QN::~Domain_D3QN()
{
	
}

/*****************************************************************************************************************************/
//D3Q15
Domain_D3Q15::Domain_D3Q15(CLbmCase* pCase)
{
	//Initialize materials in domain
	if( pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		cgsize_t twoPhase(2);
		for(cgsize_t i=0 ; i < twoPhase; i++)
		{
			//if material is incompressible
			m_domainMaterial.push_back(new D3Q15Incomp_domain(pCase, i));
		}
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainMaterial (D3Q15) contains = " << m_domainMaterial.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

		//initialize variables
		this->initializeVars(pCase, this->GetDomainVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainVariables contains = " << m_domainVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//Initialize Temporary variables into a different container
		this->initializeTemps(pCase, this->GetDomainTempVariables());
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" initialized materials and variables for multiphase..\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else
	{
		//if material is incompressible
		m_domainMaterial.push_back(new D3Q15Incomp_domain(pCase, 0));
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainMaterial (D3Q15) contains = " << m_domainMaterial.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//initialize variables
		this->initializeVars(pCase, this->GetDomainVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainVariables contains = " << m_domainVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

		//SCMP--------------------------------------------------------------------------------------------------------
		//Initialize Temporary variables into a different container
		this->initializeTemps(pCase, this->GetDomainTempVariables());
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" initialized materials and variables for SCMP\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}

	this->mapBounds(pCase);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" Finished all bound mapping\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	this->mapPeriodicNeighbor(pCase);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" Finished all periodic mapping\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	/*for(MaterialType_t::size_type i=0 ; i < GetDomainMaterial()->size(); i++)
	{
		for(PdfDomain::LatticeType_t::size_type j =1; j < GetDomainMaterial()->at(i)->GetLatticePdf()->size(); j++)
		{
			for(PdfBlock::NodeType_t::size_type k = 1; k <= GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->size(); k++)
			{
				///////////////////////////////////////////////////////////////////////////////////////////////////////////
				std::cout<<"node = "<< k <<"periodic neighbor = "<<GetDomainMaterial()->at(i)->GetLatticePdf()->at(j)->pdf()->GetPdfunction()->at(k-1)->m_periodicNeighborNum<<"\n";
				std::cin.get();
				/////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
		}
	}*/
	this->mapBoundType(pCase);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" Finished all Boundtype mapping\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/				
	this->mapStreamingStart(pCase);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" Finished all streamstart mapping\n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
}

Domain_D3Q15::Domain_D3Q15(CLbmCase* pCase, cgsize_t nvertex)
{
	//Initialize materials in domain
	//if material is incompressible
	m_domainMaterial.push_back(new D3Q15Incomp_domain(pCase, nvertex));
	//Initialize domain variables
	m_domainVariables.push_back(new CXVelocity(pCase, nvertex));
	m_domainVariables.push_back(new CYVelocity(pCase, nvertex));
	m_domainVariables.push_back(new CZVelocity(pCase, nvertex));
	m_domainVariables.push_back(new CVelocityMag(pCase, nvertex));
	m_domainVariables.push_back(new CDensity(pCase, nvertex));
}
//initializing methods-maps boundary nodes

void
Domain_D3Q15::mapStreamingStart(CLbmCase* pCase)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" map begining of streaming \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7),
					 f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), 
					 f14(14);
	for(MaterialType_t::size_type i=0 ; i < m_domainMaterial.size(); i++)
	{
		m_domainMaterial.at(i)->GetLatticePdf()->at(f1)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f2)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f2)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f1)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f3)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f4)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f4)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f3)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f5)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f6)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f6)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f5)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f7)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f14)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f14)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f7)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f11)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f10)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f10)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f11)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f9)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f12)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f12)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f9)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f13)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f8)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		m_domainMaterial.at(i)->GetLatticePdf()->at(f8)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f13)->pdf());//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	}
}

void 
Domain_D3Q15::mapCardinalDirection(cgsize_t node)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" map cardinal directions for D3Q15 \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//if node is a surface node
	if(mapCardinalSurface(node));
	//else if node is an edge node
	else if(mapCardinalEdge(node));
	//else if node is a corner node
	else if(mapCardinalCorner(node));
	else;
}

bool
Domain_D3Q15::mapCardinalSurface(cgsize_t node)
{
	cgsize_t first (0);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6);
	if(mapSurface(f1,f3,f4,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WEST;
		copyToOthers(node, Node::WEST);
		return true;
	}
	else if(mapSurface(f2,f3,f4,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EAST;
		copyToOthers(node, Node::EAST);
		return true;
	}
	else if(mapSurface(f3,f1,f2,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::SOUTH;
		copyToOthers(node, Node::SOUTH);
		return true;
	}
	else if(mapSurface(f4,f1,f2,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::NORTH;
		copyToOthers(node, Node::NORTH);
		return true;
	}
	else if(mapSurface(f5,f1,f2,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::BOTTOM;
		copyToOthers(node, Node::BOTTOM);
		return true;
	}
	else if(mapSurface(f6,f1,f2,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::TOP;
		copyToOthers(node, Node::TOP);
		return true;
	}
	else
		return false;
}

bool 
Domain_D3Q15::mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t node)
{
	cgsize_t primary(0);
	
	if(  m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
		)
		return true;
	else
		return false;
}


bool 
Domain_D3Q15::mapCardinalEdge(cgsize_t node)
{
	cgsize_t first(0);
	cgsize_t   f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	if(mapEdge(f7,f14,f9,f12,f11,f13,f3,f4,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WT;
		copyToOthers(node, Node::WT);
		return true;
	}

	else if(mapEdge(f11,f10,f13,f8,f14,f12,f3,f4,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ET;
		copyToOthers(node, Node::ET);
		return true;
	}
	else if(mapEdge(f11,f10,f13,f8,f7,f9,f3,f4,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WB;
		copyToOthers(node, Node::WB);
		return true;
	}
	else if(mapEdge(f7,f14,f9,f12,f10,f8,f3,f4,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EB;
		copyToOthers(node, Node::EB);
		return true;
	}
	/////////////////////////////////////////////////////////////////////////
	else if(mapEdge(f11,f10,f9,f12,f14,f13,f1,f2,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::NT;
		copyToOthers(node, Node::NT);
		return true;
	}
	else if(mapEdge(f7,f14,f13,f8,f11,f12,f1,f2,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ST;
		copyToOthers(node, Node::ST);
		return true;
	}
	else if(mapEdge(f7,f14,f13,f8,f10,f9,f1,f2,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::NB;
		copyToOthers(node, Node::NB);
		return true;
	}
	else if(mapEdge(f11,f10,f9,f12,f7,f8,f1,f2,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::SB;
		copyToOthers(node, Node::SB);
		return true;
	}
	//////////////////////////////////////////////////////////////////////////
	else if(mapEdge(f7,f14,f11,f10,f9,f13,f5,f6,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WN;
		copyToOthers(node, Node::WN);
		return true;
	}
	else if(mapEdge(f9,f12,f13,f8,f7,f11,f5,f6,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WS;
		copyToOthers(node, Node::WS);
		return true;
	}
	else if(mapEdge(f9,f12,f13,f8,f14,f10,f5,f6,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EN;
		copyToOthers(node, Node::EN);
		return true;
	}
	else if(mapEdge(f7,f14,f11,f10,f12,f8,f5,f6,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ES;
		copyToOthers(node, Node::ES);
		return true;
	}
	else
		return false;
}

bool 
Domain_D3Q15::mapEdge(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									 cgsize_t fV, cgsize_t fVI, cgsize_t f_orthoI,cgsize_t f_orthoII, cgsize_t node)
{
	cgsize_t primary(0);
	if(  m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0 
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fVI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(f_orthoI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(f_orthoII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0)
		return true;
	else
		return false;
}

bool 
Domain_D3Q15::mapCardinalCorner(cgsize_t node)
{
	cgsize_t first (0);
	cgsize_t  f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	if(mapCorner(f7,f14,f9,f12,f13,f8,f10, node)) //transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WST;
		copyToOthers(node, Node::WST);
		return true;
	}
	else if(mapCorner(f7,f14,f11,f10,f13,f8,f9,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EST;
		copyToOthers(node, Node::EST);
		return true;
	}
	else if(mapCorner(f7,f14,f11,f10,f9,f12,f8,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WNT;
		copyToOthers(node, Node::WNT);
		return true;
	}
	else if(mapCorner(f11,f10,f9,f12,f13,f8,f7,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ENT;
		copyToOthers(node, Node::ENT);
		return true;
	}
	else if(mapCorner(f7,f14,f11,f10,f13,f8,f12,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WNB;
		copyToOthers(node, Node::WNB);
		return true;
	}
	else if(mapCorner(f7,f14,f11,f10,f9,f12,f13,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ESB;
		copyToOthers(node, Node::ESB);
		return true;
	}
	else if(mapCorner(f11,f10,f9,f12,f13,f8,f14,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WSB;
		copyToOthers(node, Node::WSB);
		return true;
	}
	else if(mapCorner(f7,f14,f9,f12,f13,f8,f11,node))//transposed to match new lattice definition for d3q15 (d'humieres 2002)
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ENB;
		copyToOthers(node, Node::ENB);
		return true;
	}
	else
		return false;
}

bool 
Domain_D3Q15::mapCorner(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV,
									   cgsize_t fV, cgsize_t fVI, cgsize_t fVII, cgsize_t node)
{
	cgsize_t primary(0);
	if(  m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0 
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fVI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fVII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0)
		return true;
	else
		return false;
}


///44444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444444//
//Getters
Domain::MaterialType_t*
Domain_D3Q15::DomainMaterial()
{
	return &m_domainMaterial;
}

Domain::VariableType_t*
Domain_D3Q15::DomainVariables()
{
	return &m_domainVariables;
}

Domain::VariableType_t*
Domain_D3Q15::DomainTempVariables()
{
	return &m_domainTempVariables;
}

Domain::PatchNodeType_t*
Domain_D3Q15::DomainSolidObjectNodes()
{
	return &m_domainSolidObjectNodes;
}

Domain::MicroPillarsType_t*
Domain_D3Q15::DomainMicroPillars()
{
	return &m_domainMicroPillars;
}

Domain_D3Q15::~Domain_D3Q15()
{
	//place holder
}

Domain_D3Q19::Domain_D3Q19(CLbmCase* pCase)
{
	//Initialize materials in domain
	if( pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		cgsize_t twoPhase(2);
		for(cgsize_t i=0 ; i < twoPhase; i++)
		{
			//if material is incompressible
			m_domainMaterial.push_back(new D3Q19Incomp_domain(pCase, i));
		}
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainMaterial contains = " << m_domainMaterial.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

		//initialize variables
		this->initializeVars(pCase, this->GetDomainVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainVariables contains = " << m_domainVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//Initialize Temporary variables into a different container
		this->initializeTemps(pCase, this->GetDomainTempVariables());
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" initialized materials and variables for multiphase..\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else
	{
		//if material is incompressible
		m_domainMaterial.push_back(new D3Q19Incomp_domain(pCase, 0));
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainMaterial contains = " << m_domainMaterial.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//initialize variables
		this->initializeVars(pCase, this->GetDomainVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainVariables contains = " << m_domainVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}

	this->mapBounds(pCase);
	this->mapPeriodicNeighbor(pCase);
	this->mapBoundType(pCase);				
	this->mapStreamingStart(pCase);
}

Domain_D3Q19::Domain_D3Q19(CLbmCase* pCase, cgsize_t nvertex)
{
	//Initialize materials in domain
	if( pCase->m_model.m_phaseModel == CLbmModel::MULTIPHASE)
	{
		cgsize_t twoPhase(2);
		for(cgsize_t i=0 ; i < twoPhase; i++)
		{
			//if material is incompressible
			m_domainMaterial.push_back(new D3Q19Incomp_domain(pCase,i,nvertex));
		}
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainMaterial contains = " << m_domainMaterial.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

		//initialize variables
		this->initializeVars(pCase, this->GetDomainVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainVariables contains = " << m_domainVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//Initialize Temporary variables into a different container
		this->initializeTemps(pCase, this->GetDomainTempVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainTempVariables contains = " << m_domainTempVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" initialized materials and variables for multiphase..\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else
	{
		//if material is incompressible
		m_domainMaterial.push_back(new D3Q19Incomp_domain(pCase, 0, nvertex));
		/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainMaterial contains = " << m_domainMaterial.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		//initialize variables
		this->initializeVars(pCase, this->GetDomainVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainVariables contains = " << m_domainVariables.size() << " elements\n";
		std::cin.get();
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////*/
		this->initializeTemps(pCase, this->GetDomainTempVariables());
		/*////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"domainTempVariables contains = " << m_domainTempVariables.size() << " elements\n";
		std::cin.get();
		////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
}

//initializing methods-maps boundary nodes

void
Domain_D3Q19::mapStreamingStart(CLbmCase* pCase)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" map begining of streaming \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7),
					 f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), 
					 f14(14), f15(15), f16(16), f17(17), f18(18);
	for(MaterialType_t::size_type i=0 ; i < m_domainMaterial.size(); i++)
	{
		m_domainMaterial.at(i)->GetLatticePdf()->at(f1)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f2)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f2)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f1)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f3)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f4)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f4)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f3)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f5)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f6)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f6)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f5)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f7)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f10)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f8)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f9)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f9)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f8)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f10)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f7)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f11)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f14)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f12)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f13)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f13)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f12)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f14)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f11)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f15)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f18)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f16)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f17)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f17)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f16)->pdf());
		m_domainMaterial.at(i)->GetLatticePdf()->at(f18)->pdf()->setStreamBound(pCase, m_domainMaterial.at(i)->GetLatticePdf()->at(f15)->pdf());
	}
}
void 
Domain_D3Q19::mapCardinalDirection(cgsize_t node)
{
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" map cardinal directions \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//if node is a surface node
	if(mapCardinalSurface(node));
	//else if node is an edge node
	else if(mapCardinalEdge(node));
	//else if node is a corner node
	else if(mapCardinalCorner(node));
	else;
}

bool
Domain_D3Q19::mapCardinalSurface(cgsize_t node)
{
	cgsize_t first (0);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6);
	if(mapSurface(f1,f3,f4,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WEST;
		copyToOthers(node, Node::WEST);
		return true;
	}
	else if(mapSurface(f2,f3,f4,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EAST;
		copyToOthers(node, Node::EAST);
		return true;
	}
	else if(mapSurface(f3,f1,f2,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::SOUTH;
		copyToOthers(node, Node::SOUTH);
		return true;
	}
	else if(mapSurface(f4,f1,f2,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::NORTH;
		copyToOthers(node, Node::NORTH);
		return true;
	}
	else if(mapSurface(f5,f1,f2,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::BOTTOM;
		copyToOthers(node, Node::BOTTOM);
		return true;
	}
	else if(mapSurface(f6,f1,f2,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::TOP;
		copyToOthers(node, Node::TOP);
		return true;
	}
	else
		return false;
}

bool
Domain_D3Q19::mapSurface(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t node)
{
	cgsize_t primary(0);
	
	if(  m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   == 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
		)
		return true;
	else
		return false;
}


bool 
Domain_D3Q19::mapCardinalEdge(cgsize_t node)
{
	cgsize_t first(0);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7),
					 f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), 
					 f14(14), f15(15), f16(16), f17(17), f18(18);

	if(mapEdge(f1,f6,f13,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WT;
		copyToOthers(node, Node::WT);
		return true;
	}

	else if(mapEdge(f2,f6,f14,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ET;
		copyToOthers(node, Node::ET);
		return true;
	}
	else if(mapEdge(f1,f5,f11,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WB;
		copyToOthers(node, Node::WB);
		return true;
	}
	else if(mapEdge(f2,f5,f12,f3,f4,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EB;
		copyToOthers(node, Node::EB);
		return true;
	}
	/////////////////////////////////////////////////////////////////////////
	else if(mapEdge(f4,f6,f18,f1,f2,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::NT;
		copyToOthers(node, Node::NT);
		return true;
	}
	else if(mapEdge(f3,f6,f16,f1,f2,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ST;
		copyToOthers(node, Node::ST);
		return true;
	}
	else if(mapEdge(f4,f5,f17,f1,f2,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::NB;
		copyToOthers(node, Node::NB);
		return true;
	}
	else if(mapEdge(f3,f5,f15,f1,f2,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::SB;
		copyToOthers(node, Node::SB);
		return true;
	}
	//////////////////////////////////////////////////////////////////////////
	else if(mapEdge(f1,f4,f8,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WN;
		copyToOthers(node, Node::WN);
		return true;
	}
	else if(mapEdge(f1,f3,f7,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WS;
		copyToOthers(node, Node::WS);
		return true;
	}
	else if(mapEdge(f2,f4,f10,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EN;
		copyToOthers(node, Node::EN);
		return true;
	}
	else if(mapEdge(f2,f3,f9,f5,f6,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ES;
		copyToOthers(node, Node::ES);
		return true;
	}
	else
		return false;
}

bool 
Domain_D3Q19::mapEdge (cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t fIV, cgsize_t fV, cgsize_t node)
{
	cgsize_t primary(0);
	
	if( (			m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   == 0
				&&	m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
				&&	m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
			) //concave
			||
			(
						m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
				&&	m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
				&&	m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
			) //convex
		)
	{
		if(m_domainMaterial.at(primary)->GetLatticePdf()->at(fIV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
		&& m_domainMaterial.at(primary)->GetLatticePdf()->at(fV)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
		)
			return true;
		else
			return false;
	}
	else
		return false;
}
bool 
Domain_D3Q19::mapCardinalCorner(cgsize_t node)
{
	cgsize_t first (0);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7),
					 f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), 
					 f14(14), f15(15), f16(16), f17(17), f18(18);
	if(mapCorner(f1,f3,f6,f7,f13,f16,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WST;
		copyToOthers(node, Node::WST);
		return true;
	}
	else if(mapCorner(f2,f3,f6,f9,f14,f16,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::EST;
		copyToOthers(node, Node::EST);
		return true;
	}
	else if(mapCorner(f1,f4,f6,f8,f13,f18,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WNT;
		copyToOthers(node, Node::WNT);
		return true;
	}
	else if(mapCorner(f2,f4,f6,f10,f14,f18,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ENT;
		copyToOthers(node, Node::ENT);
		return true;
	}
	else if(mapCorner(f1,f4,f5,f8,f11,f18,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WNB;
		copyToOthers(node, Node::WNB);
		return true;
	}
	else if(mapCorner(f2,f3,f5,f9,f12,f15,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ESB;
		copyToOthers(node, Node::ESB);
		return true;
	}
	else if(mapCorner(f1,f3,f5,f7,f11,f15,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::WSB;
		copyToOthers(node, Node::WSB);
		return true;
	}
	else if(mapCorner(f2,f4,f5,f10,f12,f17,node))
	{
		m_domainVariables.at(first)->GetVariable()->at(first)->at(node-1)->m_bound = Node::ENB;
		copyToOthers(node, Node::ENB);
		return true;
	}
	else
		return false;
}

bool 
Domain_D3Q19::mapCorner(cgsize_t fI, cgsize_t fII, cgsize_t fIII, cgsize_t ffI, cgsize_t ffII, cgsize_t ffIII, cgsize_t node)
{
	cgsize_t primary(0);
	
		if(		m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   == 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  == 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo == 0
			)
			return true;

			  //concave 
		/*else if(
					m_domainMaterial.at(primary)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo    != 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(fIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(ffI)->pdf()->GetPdfunction()->at(node-1)->m_connectTo   != 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(ffII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo  != 0
			&&  m_domainMaterial.at(primary)->GetLatticePdf()->at(ffIII)->pdf()->GetPdfunction()->at(node-1)->m_connectTo != 0
			)    //convex            //PROBLEM IMPLEMENTING CONVEX CORNER FIX LATER!
			return true;*/
	else
		return false;
}

//Getters
inline Domain::MaterialType_t*
Domain_D3Q19::DomainMaterial()
{
	return &m_domainMaterial;
}

inline Domain::VariableType_t*
Domain_D3Q19::DomainVariables()
{
	return &m_domainVariables;
}

inline Domain::VariableType_t*
Domain_D3Q19::DomainTempVariables()
{
	return &m_domainTempVariables;
}

Domain::PatchNodeType_t*
Domain_D3Q19::DomainSolidObjectNodes()
{
	return &m_domainSolidObjectNodes;
}

Domain::MicroPillarsType_t*
Domain_D3Q19::DomainMicroPillars()
{
	return &m_domainMicroPillars;
}

Domain_D3Q19::~Domain_D3Q19()
{
	//place holder
}

//LEVEL TWO:invoked from Domain
//==============================================================================================//
//Derived class
//==============//
//D3QNIncomp_domain
void 
D3QNIncomp_domain::mapBoundAndNeighbor(CLbmCase* pCase, LatticeType_t* lattice)
{
	//map neighbor nodes
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" map neighbors \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	for(PdfDomain::LatticeType_t::size_type i =0; i < lattice->size(); i++)
		{
			lattice->at(i)->GetNeighbor()(GetLatticePdf()->at(i)->pdf(),this,GetMaterialIndex());
		}
}

void
D3QNIncomp_domain::updateAllSummations(CLbmCase* pCase, Domain& domainVariables, cgsize_t mat_index, cgsize_t node, PdfDomain::LatticeType_t::size_type f)
{
	cgsize_t xvel(0), yvel(1), zvel(2), rho(4);
	
	/*if(node == 620000)
	{
	std::cout<<" previous sum for density"<<domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" previous sum for ux(comp)"<<domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" previous sum for uy(comp)"<<domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" previous sum for uz(comp)"<<domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cin.get();
	}*/
	
	if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_isSolid == false)
	{
		//if(node == 8000)
			//std::cout<<" pdf: "<<f<<" for node = "<<node<<" (before) uploading is "<< GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			
		//update density
		domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal +=
			GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal;
			
		//sum fix * cx
		domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    += 
			(GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal * GetLatticePdf()->at(f)->pdf()->GetCx());
		//sum fiy * cy
		domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    += 
			(GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal * GetLatticePdf()->at(f)->pdf()->GetCy());
		//sum fiz * cz
		domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    += 
			(GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal * GetLatticePdf()->at(f)->pdf()->GetCz());

	}
	else;
	
	GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = false;
	GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeStreamed = false;
	
	/*if(node == 620000)
	{
	std::cout<<" pdf_new_"<<f<<" = "<<GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" pdf_old_"<<f<<" = "<<GetLatticePdf()->at(f)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal<<"\n";
	std::cout<<" running sum for density"<<domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" running sum for ux(comp)"<<domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" running sum for uy(comp)"<<domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cout<<" running sum for uz(comp)"<<domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal<<"\n";
	std::cin.get();
	}*/
}

void 
D3QNIncomp_domain::setNewToOld(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	cgsize_t xvel(0), yvel(1), zvel(2), rho(4);

	domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

	domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

	domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;

	domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_oldNodeVal = 
				domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal;
}

void 
D3QNIncomp_domain::resetToZero(Domain& domainVariables, cgsize_t mat_index, cgsize_t node)
{
	cgsize_t xvel(0), yvel(1), zvel(2), density(4);
	domainVariables.GetDomainVariables()->at(density)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal = 0.0;
	//reset before summation old = new, new = 0.0
	domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    = 0.0;
	//reset before summation old = new, new = 0.0
	domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    = 0.0;
	//reset before summation old = new, new = 0.0
	domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(mat_index)->at(node-1)->m_nodeVal    = 0.0;
}

D3QNIncomp_domain::~D3QNIncomp_domain()
{
}

//D3Q15Incomp_domain

D3Q15Incomp_domain::D3Q15Incomp_domain(CLbmCase* pCase, cgsize_t matIndex)
{
	m_materialIndex = matIndex;
	m_latticeType.push_back(new D3Q1500_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1501_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1502_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1503_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1504_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1505_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1506_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1507_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1508_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1509_incompDomain(pCase));
	m_latticeType.push_back(new D3Q15010_incompDomain(pCase));
	m_latticeType.push_back(new D3Q15011_incompDomain(pCase));
	m_latticeType.push_back(new D3Q15012_incompDomain(pCase));
	m_latticeType.push_back(new D3Q15013_incompDomain(pCase));
	m_latticeType.push_back(new D3Q15014_incompDomain(pCase));

	//map Nodes
	this->mapBoundAndNeighbor(pCase, LatticePdf());
}

D3Q15Incomp_domain::D3Q15Incomp_domain(CLbmCase* pCase, cgsize_t matIndex, cgsize_t  nvertex)
{
	m_latticeType.push_back(new D3Q1500_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1501_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1502_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1503_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1504_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1505_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1506_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1507_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1508_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1509_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q15010_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q15011_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q15012_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q15013_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q15014_incompDomain(pCase, nvertex));
}

inline PdfDomain::LatticeType_t*
D3Q15Incomp_domain::LatticePdf()
{
	return &m_latticeType;
}

inline cgsize_t 
D3Q15Incomp_domain::MaterialIndex()
{
	return m_materialIndex;
}

void 
D3Q15Incomp_domain::mapStreamStart(CLbmCase* pCase)
{
	//map start of streaming
	for(LatticeType_t::size_type i =1; i < m_latticeType.size(); i++)
		{
			m_latticeType.at(i)->pdf()->setStreamBound(pCase, m_latticeType.at(i+1)->pdf());
			m_latticeType.at(i+1)->pdf()->setStreamBound(pCase, m_latticeType.at(i)->pdf());
			i += 1;
		}
}
void 
D3Q15Incomp_domain::momentToVelocitySpace(CLbmCase* pCase, Domain& domainVariables)
{
	cgsize_t node(0), xvel(0), single(0);
	cgsize_t end = pCase->m_grid.m_NVertex;

#pragma omp for private(node)
	for(node = 1; node <= end; node++)
	{
		cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 	f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);

	
		if(domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(single)->at(node-1)->m_isSolid == false)
		{
			m_latticeType.at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn0(node-1);
			m_latticeType.at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn1(node-1);	
			m_latticeType.at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn2(node-1);
			m_latticeType.at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn3(node-1);
			m_latticeType.at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn4(node-1);
			m_latticeType.at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn5(node-1);
			m_latticeType.at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn6(node-1);
			m_latticeType.at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn7(node-1);
			m_latticeType.at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn8(node-1);
			m_latticeType.at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn9(node-1);

			m_latticeType.at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn10(node-1);
			m_latticeType.at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn11(node-1);
			m_latticeType.at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn12(node-1);
			m_latticeType.at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn13(node-1);
			m_latticeType.at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal -=  convertColumn14(node-1);

			/*if(node == 620000)
			{
			std::cout<< " f0_new = "<<m_latticeType.at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f1_new = "<<m_latticeType.at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f2_new = "<<m_latticeType.at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f3_new = "<<m_latticeType.at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f4_new = "<<m_latticeType.at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f5_new = "<<m_latticeType.at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f6_new = "<<m_latticeType.at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f7_new = "<<m_latticeType.at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f8_new = "<<m_latticeType.at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f9_new = "<<m_latticeType.at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f10_new = "<<m_latticeType.at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f11_new = "<<m_latticeType.at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f12_new = "<<m_latticeType.at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f13_new = "<<m_latticeType.at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cout<< " f14_new = "<<m_latticeType.at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
			std::cin.get();
			}*/
			//store old value for f0 since it does not stream 
			//m_latticeType.at(f0)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal = m_latticeType.at(f0)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal;
		}
	}
}

Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn0(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2);

	sum =	(0.0667) * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  //GetPdfmoment()->at(i-1)->m_nodeVal
		(-0.1111)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0444) * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal;

	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn1(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f9(9);

	sum =	(0.0667)   * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0111)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)   * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000)  * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1667)   * m_latticeType.at(f9)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn2(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f9(9);

	sum =	(0.0667)   * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.01111) * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000)  * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)   * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1667)   * m_latticeType.at(f9)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn3(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f5(5), f6(6), f9(9), f10(10);


	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0556) * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0111) * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)  * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0833) * m_latticeType.at(f9)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.2500)  * m_latticeType.at(f10)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn4(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f5(5), f6(6), f9(9), f10(10);

	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0556) * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0111) * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000) * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0833) * m_latticeType.at(f9)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.2500)  * m_latticeType.at(f10)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn5(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f7(7), f8(8), f9(9), f10(10);
	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0556) * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0111) * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)  * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0833) * m_latticeType.at(f9)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.2500) * m_latticeType.at(f10)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn6(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f7(7), f8(8), f9(9), f10(10);

	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0556) * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(-0.0111) * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000) * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0833) * m_latticeType.at(f9)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.2500) * m_latticeType.at(f10)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn7(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);
	
	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)  * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(0.1250) * m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn8(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);	

	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000) * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(-0.1250)* m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn9(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);
	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)  * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(-0.1250)* m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn10(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);
	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000) * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(0.1250) * m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn11(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);
	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)  * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(0.1250) * m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn12(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);
	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000) * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1000)  * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(-0.1250)* m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn13(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);

	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(0.1000)  * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.0250)  * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(-0.1250)* m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
Node::NodeValueType_t
D3Q15Incomp_domain::convertColumn14(cgsize_t node)
{
	Node::NodeValueType_t sum(0.0);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6),
		 f7(7), f8(8), f11(11), f12(12), f13(13), f14(14);

	
	sum =	(0.0667)  * m_latticeType.at(f0)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0556)  * m_latticeType.at(f1)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +  
		(0.0028)  * m_latticeType.at(f2)->pdf()->GetPdfmoment()->at(node)->m_nodeVal + 
		(-0.1000) * m_latticeType.at(f3)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f4)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f5)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f6)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1000) * m_latticeType.at(f7)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.0250) * m_latticeType.at(f8)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +

		(0.1250) * m_latticeType.at(f11)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f12)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(0.1250) * m_latticeType.at(f13)->pdf()->GetPdfmoment()->at(node)->m_nodeVal +
		(-0.1250)* m_latticeType.at(f14)->pdf()->GetPdfmoment()->at(node)->m_nodeVal ;
	/*if(node == 619999)
	{
		std::cout<<" Minv*DeltaPDF   = "<<sum<<"\n";
		std::cin.get();
	}*/
	return sum;
}
D3Q15Incomp_domain::~D3Q15Incomp_domain()
{
}

//D3Q19Incomp_domain
D3Q19Incomp_domain::D3Q19Incomp_domain(CLbmCase* pCase, cgsize_t matIndex)
{
	m_materialIndex = matIndex;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" material =  "<< m_materialIndex <<"\n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_latticeType.push_back(new D3Q1900_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1901_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1902_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1903_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1904_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1905_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1906_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1907_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1908_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1909_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1910_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1911_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1912_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1913_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1914_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1915_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1916_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1917_incompDomain(pCase));
	m_latticeType.push_back(new D3Q1918_incompDomain(pCase));

	//map neighbor Nodes
	this->mapBoundAndNeighbor(pCase, LatticePdf());
}

//D3Q19
D3Q19Incomp_domain::D3Q19Incomp_domain(CLbmCase* pCase, cgsize_t matIndex, cgsize_t  nvertex)
{
	//place holder
	m_materialIndex = matIndex;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" material =  "<< m_materialIndex <<"\n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	m_latticeType.push_back(new D3Q1900_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1901_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1902_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1903_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1904_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1905_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1906_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1907_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1908_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1909_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1910_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1911_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1912_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1913_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1914_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1915_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1916_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1917_incompDomain(pCase, nvertex));
	m_latticeType.push_back(new D3Q1918_incompDomain(pCase, nvertex));
}

//velocities at boundary for pressure outlet boundary condition
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline PdfDomain::LatticeType_t*
D3Q19Incomp_domain::LatticePdf()
{
	return &m_latticeType;
}

inline cgsize_t 
D3Q19Incomp_domain::MaterialIndex()
{
	return m_materialIndex;
}

void 
D3Q19Incomp_domain::mapStreamStart(CLbmCase* pCase)
{
	//place holder
}

D3Q19Incomp_domain::~D3Q19Incomp_domain()
{
}
//LEVEL THREE: invoked from 
//==============================================================================================//
/*--------------------------------------------------------------------------------------------*/
//D3q15_incompDomain Abstract class
//Properties:
//1. Derived from LbmDomain and base to D3q1501_incompDomain...
/*--------------------------------------------------------------------------------------------*/
PdfBlock*
D3qN_incompDomain::maker (PdfMaker* pimpl, CLbmCase* pCase)
{

	 return pimpl->incomp(pCase);
}
PdfBlock*
D3qN_incompDomain::maker (PdfMaker* pimpl, CLbmCase* pCase, cgsize_t nvertex)
{

	 return pimpl->incomp(pCase, nvertex);
}
/*--------------------------------------------------------------------------------------------*/
//D3Q1500_incompDomain concrete class
//Properties:
//1. Derived from D3Q15Incomp 
//2. Implementation of process specific to D3q1500Incomp objects
/*--------------------------------------------------------------------------------------------*/
D3Q1500_incompDomain::D3Q1500_incompDomain(CLbmCase* pCase)
{
	pdf_00Domain = maker(Impl(), pCase);
}

D3Q1500_incompDomain::D3Q1500_incompDomain(CLbmCase* pCase, cgsize_t nvertex)
{
	pdf_00Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1500_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}
Pdf_setNeighbor&
D3Q1500_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D3Q1500_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}
Pdf_bounder&
D3Q1500_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1500_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

PdfBlock* 
D3Q1500_incompDomain::pdfDomain() 
{
	return pdf_00Domain;
}

PdfMaker* 
D3Q1500_incompDomain::implDomain() 
{
	return &impl_00;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1501_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1501_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1501_incompDomain::D3Q1501_incompDomain(CLbmCase* pCase)
{
	pdf_01Domain = maker(Impl(), pCase);
}
D3Q1501_incompDomain::D3Q1501_incompDomain(CLbmCase* pCase, cgsize_t nvertex)
{
	pdf_01Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1501_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_setNeighbor&
D3Q1501_incompDomain::GetNeighbor()
{
	return neighbor;
}

Pdf_bounder&
D3Q1501_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1501_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1501_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
void 
D3Q1501_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f2(2);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

PdfBlock*  
D3Q1501_incompDomain::pdfDomain()
{
	return pdf_01Domain;
}

PdfMaker* 
D3Q1501_incompDomain::implDomain()
{
	return &impl_01;
}


/*--------------------------------------------------------------------------------------------*/
//D3Q1502_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1502_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1502_incompDomain::D3Q1502_incompDomain(CLbmCase* pCase) 
{	
	pdf_02Domain = maker(Impl(), pCase);
}
D3Q1502_incompDomain::D3Q1502_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_02Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1502_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1502_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1502_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1502_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1502_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1502_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f1(1);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1502_incompDomain::pdfDomain()
{
	return pdf_02Domain;
}

PdfMaker* 
D3Q1502_incompDomain::implDomain()
{
	return &impl_02;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1503_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1503_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1503_incompDomain::D3Q1503_incompDomain(CLbmCase* pCase) 
{	
	pdf_03Domain = maker(Impl(), pCase);
}
D3Q1503_incompDomain::D3Q1503_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_03Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1503_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1503_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1503_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1503_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1503_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1503_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f4(4);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1503_incompDomain::pdfDomain()
{
	return pdf_03Domain;
}

PdfMaker* 
D3Q1503_incompDomain::implDomain()
{
	return &impl_03;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1504_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1504_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1504_incompDomain::D3Q1504_incompDomain(CLbmCase* pCase) 
{	
	pdf_04Domain = maker(Impl(), pCase);
}
D3Q1504_incompDomain::D3Q1504_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_04Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1504_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1504_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1504_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1504_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1504_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1504_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f3(3);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1504_incompDomain::pdfDomain()
{
	return pdf_04Domain;
}

PdfMaker* 
D3Q1504_incompDomain::implDomain()
{
	return &impl_04;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1505_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1505_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1505_incompDomain::D3Q1505_incompDomain(CLbmCase* pCase) 
{	
	pdf_05Domain = maker(Impl(), pCase);
}
D3Q1505_incompDomain::D3Q1505_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_05Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1505_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1505_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1505_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1505_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1505_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1505_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f6(6);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1505_incompDomain::pdfDomain()
{
	return pdf_05Domain;
}

PdfMaker* 
D3Q1505_incompDomain::implDomain()
{
	return &impl_05;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1506_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1506_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1506_incompDomain::D3Q1506_incompDomain(CLbmCase* pCase) 
{	
	pdf_06Domain = maker(Impl(), pCase);
}
D3Q1506_incompDomain::D3Q1506_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_06Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1506_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1506_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1506_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1506_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1506_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1506_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f5(5);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1506_incompDomain::pdfDomain()
{
	return pdf_06Domain;
}

PdfMaker* 
D3Q1506_incompDomain::implDomain()
{
	return &impl_06;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1507_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1507_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1507_incompDomain::D3Q1507_incompDomain(CLbmCase* pCase) 
{	
	pdf_07Domain = maker(Impl(), pCase);
}
D3Q1507_incompDomain::D3Q1507_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_07Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1507_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1507_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1507_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1507_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1507_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1507_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f8(8);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}
PdfBlock*  
D3Q1507_incompDomain::pdfDomain()
{
	return pdf_07Domain;
}

PdfMaker* 
D3Q1507_incompDomain::implDomain()
{
	return &impl_07;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1508_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1508_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1508_incompDomain::D3Q1508_incompDomain(CLbmCase* pCase) 
{	
	pdf_08Domain = maker(Impl(), pCase);
}
D3Q1508_incompDomain::D3Q1508_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_08Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1508_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1508_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1508_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1508_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1508_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1508_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f7(7);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1508_incompDomain::pdfDomain()
{
	return pdf_08Domain;
}

PdfMaker* 
D3Q1508_incompDomain::implDomain()
{
	return &impl_08;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1509_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q1509_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1509_incompDomain::D3Q1509_incompDomain(CLbmCase* pCase) 
{	
	pdf_09Domain = maker(Impl(), pCase);
}
D3Q1509_incompDomain::D3Q1509_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_09Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1509_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q1509_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1509_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q1509_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q1509_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1509_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f10(10);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1509_incompDomain::pdfDomain()
{
	return pdf_09Domain;
}

PdfMaker* 
D3Q1509_incompDomain::implDomain()
{
	return &impl_09;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q15010_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15010_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q15010_incompDomain::D3Q15010_incompDomain(CLbmCase* pCase) 
{	
	pdf_010Domain = maker(Impl(), pCase);
}
D3Q15010_incompDomain::D3Q15010_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_010Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q15010_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q15010_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q15010_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q15010_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
Pdf_setNeighbor&
D3Q15010_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q15010_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f9(9);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

PdfBlock*  
D3Q15010_incompDomain::pdfDomain()
{
	return pdf_010Domain;
}

PdfMaker* 
D3Q15010_incompDomain::implDomain()
{
	return &impl_010;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q15011_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15011_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q15011_incompDomain::D3Q15011_incompDomain(CLbmCase* pCase) 
{	
	pdf_011Domain = maker(Impl(), pCase);
}
D3Q15011_incompDomain::D3Q15011_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_011Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q15011_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q15011_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q15011_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q15011_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q15011_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q15011_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f12(12);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q15011_incompDomain::pdfDomain()
{
	return pdf_011Domain;
}

PdfMaker* 
D3Q15011_incompDomain::implDomain()
{
	return &impl_011;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q15012_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15012_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q15012_incompDomain::D3Q15012_incompDomain(CLbmCase* pCase) 
{	
	pdf_012Domain = maker(Impl(), pCase);
}
D3Q15012_incompDomain::D3Q15012_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_012Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q15012_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q15012_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q15012_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q15012_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q15012_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q15012_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f11(11);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q15012_incompDomain::pdfDomain()
{
	return pdf_012Domain;
}

PdfMaker* 
D3Q15012_incompDomain::implDomain()
{
	return &impl_012;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q15013_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15013_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q15013_incompDomain::D3Q15013_incompDomain(CLbmCase* pCase) 
{	
	pdf_013Domain = maker(Impl(), pCase);
}
D3Q15013_incompDomain::D3Q15013_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_013Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q15013_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q15013_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q15013_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q15013_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q15013_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q15013_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f14(14);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q15013_incompDomain::pdfDomain()
{
	return pdf_013Domain;
}

PdfMaker* 
D3Q15013_incompDomain::implDomain()
{
	return &impl_013;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q15014_incompDomain concrete class
//Properties:
//1. Derived from D3q15_incompDomain 
//2. Implementation of process specific to D3Q15014_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q15014_incompDomain::D3Q15014_incompDomain(CLbmCase* pCase) 
{	
	pdf_014Domain = maker(Impl(), pCase);
}
D3Q15014_incompDomain::D3Q15014_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_014Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q15014_incompDomain::GetFullBBBounder()
{
	return wall_fbb;
}

Pdf_bounder&
D3Q15014_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q15014_incompDomain::GetVelBounder()
{
	return v_boundary;
}
Pdf_bounder&
D3Q15014_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q15014_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q15014_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f13(13);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q15014_incompDomain::pdfDomain()
{
	return pdf_014Domain;
}

PdfMaker* 
D3Q15014_incompDomain::implDomain()
{
	return &impl_014;
}
//============================================================================================//
//FUNCTORS
//base class
Node::NodeValueType_t 
D3Q15_bounder::ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
													Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound)
{
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);
	cgsize_t  rho(4), xVel(0), yVel(1), zVel(2);
	Node::NodeValueType_t density(0.0), velocity(0.0), temp(0.0);
	
	if     (b_type == LbmDomain::PRESSURE)
	{
		density = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		if(bound  == Node::WEST)
		{
			Node::NodeValueType_t c1 = -1.0;
			Node::NodeValueType_t c2 =  1.0;
			velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f3, f4, f2, f14, f10, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		}
		else if(bound == Node::EAST)
			{
				Node::NodeValueType_t c1 =  1.0;
				Node::NodeValueType_t c2 = -1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f3, f4, f1, f7, f11, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::SOUTH)
			{
				Node::NodeValueType_t c1 = -1.0;
				Node::NodeValueType_t c2 =  1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f4, f14, f10, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::NORTH)
			{
				Node::NodeValueType_t c1 =  1.0;
				Node::NodeValueType_t c2 = -1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f3, f7, f11, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::BOTTOM)
			{
				Node::NodeValueType_t c1 =  -1.0;
				Node::NodeValueType_t c2 =   1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f6, f14, f11, f3, f4, f12, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::TOP)
			{
				Node::NodeValueType_t c1 =  1.0;
				Node::NodeValueType_t c2 = -1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f5, f7, f10, f3, f4, f9, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
	}

	else if(b_type == LbmDomain::VELOCITY)
	{
		velocity = domainVariables.GetDomainVariables()->at(vel_index)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;

		//now determine the density
		if(bound  == Node::WEST)
		{
			Node::NodeValueType_t c = -1.0;
			density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f2, f14, f10, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
		}
		else if(bound == Node::EAST)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f1, f7, f11, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::SOUTH)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f4, f14, f10, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::NORTH)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f3, f7, f11, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::BOTTOM)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f6, f14, f11, f3, f4, f12, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::TOP)
			{	
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f5, f7, f10, f3, f4, f9, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
	}
	
	temp = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal + 
				 c*density*velocity;
	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"node = "<<node<<" pdf = "<< f << " Value = "<<temp <<"\n";
	std::cin.get();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	return temp;
}

Node::NodeValueType_t 
D3Q15_bounder::ZhouHe_diag(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t fI, cgsize_t fII, cgsize_t ffI, cgsize_t ffII, 
		cgsize_t node, Node::NodeValueType_t c1, Node::NodeValueType_t c2, Node::NodeValueType_t c3, Node::NodeValueType_t c4, Node::NodeValueType_t c5,
		cgsize_t xvel, cgsize_t yvel, cgsize_t zvel, LbmDomain::BOUND b_type, Node::LBMBOUND bound)
{
	cgsize_t rho(4);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);
	cgsize_t xVel(0), yVel(1), zVel(2);
	Node::NodeValueType_t u_wall(0.0), v_wall(0.0), w_wall(0.0), density(0.0), temp(0.0);
	

	if     (b_type == LbmDomain::PRESSURE)
	{
		Node::NodeValueType_t a1 (0.0), a2(0.0);
		density = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		if(  bound == Node::WEST)
			{
				a1=-1.0, a2=1.0;
				u_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f3, f4, f2, f14, f10, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(	bound == Node::EAST )
			{
				a1=1.0, a2=-1.0;
				u_wall  = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f3, f4, f1, f7, f11, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::SOUTH)
			{
				a1=-1.0, a2=1.0;
				v_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f1, f2, f4, f14, f10, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::NORTH )
			{
				a1=1.0, a2=-1.0;
				v_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f1, f2, f3, f7, f11, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::BOTTOM)
			{
				a1=-1.0, a2=1.0;
				w_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1,a2, f0, f1, f2, f6, f14, f11, f3, f4, f12, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::TOP)
			{
				a1=1.0, a2=-1.0;
				w_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1,a2, f0, f1, f2, f5, f7, f10, f3, f4, f9, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
	}
	else if(b_type == LbmDomain::VELOCITY)
	{//velocities are specified therefore we need to determine density (rho) at boundary
		u_wall = domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		v_wall = domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		w_wall = domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		if(  bound == Node::WEST || bound == Node::WS  || bound == Node::WN  || bound == Node::WB || bound == Node::WT
			|| bound == Node::WSB  || bound == Node::WST || bound == Node::WNB || bound == Node::WNT)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f2, f14, f10, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}

		else if(	bound == Node::EAST || bound == Node::ES  || bound == Node::EN  || bound == Node::EB || bound == Node::ET
			||			bound == Node::ESB  || bound == Node::EST || bound == Node::ENB || bound == Node::ENT)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f1, f7, f11, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}

		else if(bound == Node::SOUTH || bound == Node::SB  || bound == Node::ST)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f4, f14, f10, f5, f6, f9, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::NORTH || bound == Node::NB  || bound == Node::NT)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f3, f7, f11, f5, f6, f12, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::BOTTOM)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f6, f14, f11, f3, f4, f12, f13);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
		else if(bound == Node::TOP)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f5, f7, f10, f3, f4, f9, f8);//transposed to match new lattice definition for d3q15 (d'humieres 2002)
			}
	}

	temp =		(c1 * density *u_wall) + (c2 * density * v_wall) + (c3 * density * w_wall)
					+ domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal

					+ c4 * ( domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal
					       - domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal )
								 
					+ c5 * ( domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(ffI)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal
					       - domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(ffII)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal ); 	

	/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"node = "<<node<<" pdf = "<< f << " Value = "<<temp <<"\n";
	std::cin.get();
	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	return temp;
}

//derived classes
void
D3Q1500_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//NULL BODY
}
void
D3Q1500_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}
void
D3Q1500_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D3Q1500_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D3Q1500_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D3Q1501_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f2(2);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1501_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), xvel(0);
	Node::NodeValueType_t  c = (2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f1, f2, node, c, xvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1501_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), xVel(0), useNode(0);
	Node::NodeValueType_t vel(0.0);
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f1, f2, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f1, f2, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1501_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f1, f2, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f1, f2, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f1, f2, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1501_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1502_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f1(1);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}	
}
void
D3Q1502_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), xvel(0);
	Node::NodeValueType_t  c = (-2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f2, f1, node, c, xvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1502_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), xVel(0), useNode(0);
	Node::NodeValueType_t vel(0.0);
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f2, f1, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f2, f1, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D3Q1502_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index); 
	cgsize_t f1(1), f2(2);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f2, f1, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f2, f1, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f2, f1, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1502_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1503_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f4(4);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}
void
D3Q1503_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f3(3), f4(4), yvel(1);
	Node::NodeValueType_t  c = (2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f3, f4, node, c, yvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1503_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f3(3), f4(4), yVel(1), useNode(0);
	Node::NodeValueType_t vel(0.0);
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f3, f4, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f3, f4, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D3Q1503_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f3(3), f4(4);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f3, f4, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f3, f4, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f3, f4, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1503_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1504_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f3(3);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}			
}
void
D3Q1504_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f3(3), f4(4), yvel(1);
	Node::NodeValueType_t  c = (-2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f4, f3, node, c, yvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1504_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f3(3), f4(4), yVel(1), useNode(0);
	Node::NodeValueType_t vel(0.0);
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f4, f3, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f4, f3, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

	
}

void
D3Q1504_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f3(3), f4(4);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f4, f3, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f4, f3, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f4, f3, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1504_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1505_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f6(6);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}			
}
void
D3Q1505_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f5(5), f6(6), zvel(2);
	Node::NodeValueType_t  c = (2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f5, f6, node, c, zvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1505_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f5(5), f6(6), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f5, f6, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f5, f6, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D3Q1505_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f5(5), f6(6);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f5, f6, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f5, f6, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f5, f6, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1505_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1506_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f5(5);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}
void
D3Q1506_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f5(5), f6(6), zvel(2);
	Node::NodeValueType_t  c = (-2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f6, f5, node, c, zvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1506_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f5(5), f6(6), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f6, f5, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f6, f5, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D3Q1506_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f5(5), f6(6);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f6, f5, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f6, f5, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f6, f5, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1506_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1507_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f14(14);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}
void
D3Q1507_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f14(14);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WS
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WB
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WSB)
			{
				c1 = (1.0/12.0), c2 = (1.0/4.0), c3 = (1.0/4.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f7, f14, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SB)
			{
				c1 = (1.0/4.0), c2 = (1.0/12.0), c3 = (1.0/4.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f7, f14, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
			{
				c1 = (1.0/4.0), c2 = (1.0/4.0), c3 = (1.0/12.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f7, f14, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1507_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f7(7), f2(2), f4(4), f6(6), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WS   || 
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WB   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WSB    )
	{
		f_opp = f2;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH || 
		 pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SB      )
	{
		f_opp = f4;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
	{
		f_opp = f6;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;			
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f7, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f7, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D3Q1507_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f7(7), f14(14);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f7, f14, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f7, f14, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f7, f14, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1507_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15014_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f7(7);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}			
}
void
D3Q15014_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f14(14);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EN
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ET
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ENT)
			{
				c1 = (-1.0/12.0), c2 = (-1.0/4.0), c3 = (-1.0/4.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f14, f7, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NT)
			{
				c1 = (-1.0/4.0), c2 = (-1.0/12.0), c3 = (-1.0/4.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f14, f7, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
			{
				c1 = (-1.0/4.0), c2 = (-1.0/4.0), c3 = (-1.0/12.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f14, f7, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15014_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f14(14), f1(1), f3(3), f5(5), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EN   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ET   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ENT    )
	{
		f_opp = f1;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH || 
		 pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NT      )
	{
		f_opp = f3;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP )
	{
		f_opp = f5;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f14, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f14, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}

void
D3Q15014_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f7(7), f14(14);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f14, f7, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f14, f7, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f14, f7, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15014_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15011_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f10(10);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	
	}			
}

void
D3Q15011_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f11(11), f10(10);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WS
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WT
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WST)
			{
				c1 = (1.0/12.0), c2 = (1.0/4.0), c3 = (-1.0/4.0), c4 = (-1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f11, f10, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ST)
			{
				c1 = (1.0/4.0), c2 = (1.0/12.0), c3 = (-1.0/4.0), c4 = (-1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f11, f10, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
			{
				c1 = (1.0/4.0), c2 = (1.0/4.0), c3 = (-1.0/12.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f11, f10, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15011_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f11(11), f2(2), f4(4), f5(5), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WS   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WT   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WST )
	{
		f_opp = f2;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH || 
	         pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ST )
	{
		f_opp = f4;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
	{
		f_opp = f5;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f11, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f11, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15011_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f11(11), f10(10);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f11, f10, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f11, f10, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f11, f10, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15011_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15010_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f11(11);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		
	}			
}
void
D3Q15010_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f11(11), f10(10);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EN
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EB
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ENB)
			{
				c1 = (-1.0/12.0), c2 = (-1.0/4.0), c3 = (1.0/4.0), c4 = (1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f10, f11, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NB)
			{
				c1 = (-1.0/4.0), c2 = (-1.0/12.0), c3 = (1.0/4.0), c4 = (1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f10, f11, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
			{
				c1 = (-1.0/4.0), c2 = (-1.0/4.0), c3 = (1.0/12.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f10, f11, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15010_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f10(10), f1(1), f3(3), f6(6), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EN   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EB   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ENB )
	{
		f_opp = f1;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH || 
	         pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NB )
	{
		f_opp = f3;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
	{
		f_opp = f6;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f10, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f10, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15010_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f11(11), f10(10);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f10, f11, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f10, f11, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f10, f11, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15010_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1509_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f12(12);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}

void
D3Q1509_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f9(9), f12(12);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WN
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WB
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WNB)
			{
				c1 = (1.0/12.0), c2 = (-1.0/4.0), c3 = (1.0/4.0), c4 = (1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f9, f12, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NB)
			{
				c1 = (1.0/4.0), c2 = (-1.0/12.0), c3 = (1.0/4.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f9, f12, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
			{
				c1 = (1.0/4.0), c2 = (-1.0/4.0), c3 = (1.0/12.0), c4 = (-1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f9, f12, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1509_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f9(9), f2(2), f3(3), f6(6), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WN   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WB   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WNB )
	{
		f_opp = f2;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH || 
	         pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NB )
	{
		f_opp = f3;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
	{
		f_opp = f6;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f9, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f9, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1509_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f9(9), f12(12);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f9, f12, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		cgsize_t d3q15(0), f9(9), f12(12);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f9, f12, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		cgsize_t d3q15(0), f9(9), f12(12);
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f9, f12, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1509_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15012_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f9(9);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}
void
D3Q15012_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f9(9), f12(12);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ES
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ET
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EST)
			{
				c1 = (-1.0/12.0), c2 = (1.0/4.0), c3 = (-1.0/4.0), c4 = (-1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f12, f9, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ST)
			{
				c1 = (-1.0/4.0), c2 = (1.0/12.0), c3 = (-1.0/4.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f12, f9, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
			{
				c1 = (-1.0/4.0), c2 = (1.0/4.0), c3 = (-1.0/12.0), c4 = (1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f12, f9, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15012_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f12(12), f1(1), f4(4), f5(5), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ES   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ET   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EST )
	{
		f_opp = f1;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH || 
	         pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ST )
	{
		f_opp = f4;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
	{
		f_opp = f5;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f12, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f12, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15012_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f9(9), f12(12);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f12, f9, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f12, f9, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f12, f9, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15012_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15013_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f8(8);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}

void
D3Q15013_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f13(13), f8(8);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WN
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WT
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WNT)
			{
				c1 = (1.0/12.0), c2 = (-1.0/4.0), c3 = (-1.0/4.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f13, f8, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NT)
			{
				c1 = (1.0/4.0), c2 = (-1.0/12.0), c3 = (-1.0/4.0), c4 = (-1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f13, f8, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
			{
				c1 = (1.0/4.0), c2 = (-1.0/4.0), c3 = (-1.0/12.0), c4 = (-1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f13, f8, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q15013_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f13(13), f2(2), f3(3), f5(5), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WN   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WT   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WNT )
	{
		f_opp = f2;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH || 
	         pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NT )
	{
		f_opp = f3;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
	{
		f_opp = f5;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f13, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f13, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15013_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f13(13), f8(8);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f13, f8, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f13, f8, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f13, f8, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q15013_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1508_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f13(13);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}			
}

void
D3Q1508_velBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f13(13), f8(8);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ES
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EB
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ESB)
			{
				c1 = (-1.0/12.0), c2 = (1.0/4.0), c3 = (1.0/4.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f8, f13, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SB)
			{
				c1 = (-1.0/4.0), c2 = (1.0/12.0), c3 = (1.0/4.0), c4 = (1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f8, f13, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
			{
				c1 = (-1.0/4.0), c2 = (1.0/4.0), c3 = (1.0/12.0), c4 = (1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f8, f13, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1508_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f_opp(0), f8(8), f1(1), f4(4), f6(6), xVel(0), yVel(1), zVel(2), useNode(0);
	Node::NodeValueType_t vel(0.0);

	if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ES   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EB   ||
	   pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ESB )
	{
		f_opp = f1;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(xVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH || 
	         pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SB )
	{
		f_opp = f4;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(yVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}
	else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
	{
		f_opp = f6;
		useNode = domainVariables.GetDomainMaterial()->at(d3q15)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
		vel = domainVariables.GetDomainVariables()->at(zVel)->GetVariable()->at(d3q15)->at(useNode-1)->m_nodeVal;
	}

	else;
		
	if(f_type == LbmDomain::NBC)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_NBC(domainVariables, d3q15, f8, f_opp, node);
		
	}
	else if(f_type == LbmDomain::CBC_LV)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = open_CBC(domainVariables, d3q15, f8, f_opp, node, vel);
	}
	
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}


void
D3Q1508_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q15(mat_index);
	cgsize_t f13(13), f8(8);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q15, f8, f13, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q15, f8, f13, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q15, f8, f13, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1508_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*--------------------------------------------------------------------------------------------*/
//D3Q1900_incompDomain concrete class
//Properties:
//1. Derived from D3QNIncomp 
//2. Implementation of process specific to D3q1500Incomp objects
/*--------------------------------------------------------------------------------------------*/
D3Q1900_incompDomain::D3Q1900_incompDomain(CLbmCase* pCase)
{
	pdf_00Domain = maker(Impl(), pCase);
}

D3Q1900_incompDomain::D3Q1900_incompDomain(CLbmCase* pCase, cgsize_t nvertex)
{
	pdf_00Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1900_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}
Pdf_setNeighbor&
D3Q1900_incompDomain::GetNeighbor()
{
	return neighbor;
}
Pdf_bounder&
D3Q1900_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}
Pdf_bounder&
D3Q1900_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1900_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}
PdfBlock* 
D3Q1900_incompDomain::pdfDomain() 
{
	return pdf_00Domain;
}

PdfMaker* 
D3Q1900_incompDomain::implDomain() 
{
	return &impl_00;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1901_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1901_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1901_incompDomain::D3Q1901_incompDomain(CLbmCase* pCase)
{
	pdf_01Domain = maker(Impl(), pCase);
}
D3Q1901_incompDomain::D3Q1901_incompDomain(CLbmCase* pCase, cgsize_t nvertex)
{
	pdf_01Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1901_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_setNeighbor&
D3Q1901_incompDomain::GetNeighbor()
{
	return neighbor;
}

Pdf_bounder&
D3Q1901_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1901_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1901_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

void 
D3Q1901_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f2(2);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

PdfBlock*  
D3Q1901_incompDomain::pdfDomain()
{
	return pdf_01Domain;
}

PdfMaker* 
D3Q1901_incompDomain::implDomain()
{
	return &impl_01;
}


/*--------------------------------------------------------------------------------------------*/
//D3Q1902_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1902_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1902_incompDomain::D3Q1902_incompDomain(CLbmCase* pCase) 
{	
	pdf_02Domain = maker(Impl(), pCase);
}
D3Q1902_incompDomain::D3Q1902_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_02Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1902_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1902_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1902_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1902_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1902_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1902_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f1(1);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1902_incompDomain::pdfDomain()
{
	return pdf_02Domain;
}

PdfMaker* 
D3Q1902_incompDomain::implDomain()
{
	return &impl_02;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1903_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1903_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1903_incompDomain::D3Q1903_incompDomain(CLbmCase* pCase) 
{	
	pdf_03Domain = maker(Impl(), pCase);
}
D3Q1903_incompDomain::D3Q1903_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_03Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1903_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1903_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1903_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1903_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1903_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1903_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f4(4);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1903_incompDomain::pdfDomain()
{
	return pdf_03Domain;
}

PdfMaker* 
D3Q1903_incompDomain::implDomain()
{
	return &impl_03;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1904_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1904_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1904_incompDomain::D3Q1904_incompDomain(CLbmCase* pCase) 
{	
	pdf_04Domain = maker(Impl(), pCase);
}
D3Q1904_incompDomain::D3Q1904_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_04Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1904_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1904_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1904_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1904_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1904_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1904_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f3(3);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1904_incompDomain::pdfDomain()
{
	return pdf_04Domain;
}

PdfMaker* 
D3Q1904_incompDomain::implDomain()
{
	return &impl_04;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1905_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1905_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1905_incompDomain::D3Q1905_incompDomain(CLbmCase* pCase) 
{	
	pdf_05Domain = maker(Impl(), pCase);
}
D3Q1905_incompDomain::D3Q1905_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_05Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1905_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1905_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1905_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1905_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1905_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1905_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f6(6);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1905_incompDomain::pdfDomain()
{
	return pdf_05Domain;
}

PdfMaker* 
D3Q1905_incompDomain::implDomain()
{
	return &impl_05;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1906_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1906_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1906_incompDomain::D3Q1906_incompDomain(CLbmCase* pCase) 
{	
	pdf_06Domain = maker(Impl(), pCase);
}
D3Q1906_incompDomain::D3Q1906_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_06Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1906_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1906_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1906_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1906_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1906_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1906_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f5(5);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1906_incompDomain::pdfDomain()
{
	return pdf_06Domain;
}

PdfMaker* 
D3Q1906_incompDomain::implDomain()
{
	return &impl_06;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1907_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1907_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1907_incompDomain::D3Q1907_incompDomain(CLbmCase* pCase) 
{	
	pdf_07Domain = maker(Impl(), pCase);
}
D3Q1907_incompDomain::D3Q1907_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_07Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1907_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1907_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1907_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1907_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1907_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1907_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f10(10);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;

}
PdfBlock*  
D3Q1907_incompDomain::pdfDomain()
{
	return pdf_07Domain;
}

PdfMaker* 
D3Q1907_incompDomain::implDomain()
{
	return &impl_07;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1908_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1908_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1908_incompDomain::D3Q1908_incompDomain(CLbmCase* pCase) 
{	
	pdf_08Domain = maker(Impl(), pCase);
}
D3Q1908_incompDomain::D3Q1908_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_08Domain = maker(Impl(), pCase, nvertex);
}

Pdf_bounder&
D3Q1908_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1908_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1908_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1908_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1908_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1908_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f9(9);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1908_incompDomain::pdfDomain()
{
	return pdf_08Domain;
}

PdfMaker* 
D3Q1908_incompDomain::implDomain()
{
	return &impl_08;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1909_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1909_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1909_incompDomain::D3Q1909_incompDomain(CLbmCase* pCase) 
{	
	pdf_09Domain = maker(Impl(), pCase);
}
D3Q1909_incompDomain::D3Q1909_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_09Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1909_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1909_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1909_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1909_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1909_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1909_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f8(8);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1909_incompDomain::pdfDomain()
{
	return pdf_09Domain;
}

PdfMaker* 
D3Q1909_incompDomain::implDomain()
{
	return &impl_09;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1910_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1910_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1910_incompDomain::D3Q1910_incompDomain(CLbmCase* pCase) 
{	
	pdf_10Domain = maker(Impl(), pCase);
}
D3Q1910_incompDomain::D3Q1910_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_10Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1910_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1910_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1910_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1910_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1910_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1910_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f7(7);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

PdfBlock*  
D3Q1910_incompDomain::pdfDomain()
{
	return pdf_10Domain;
}

PdfMaker* 
D3Q1910_incompDomain::implDomain()
{
	return &impl_10;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1911_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1911_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1911_incompDomain::D3Q1911_incompDomain(CLbmCase* pCase) 
{	
	pdf_11Domain = maker(Impl(), pCase);
}
D3Q1911_incompDomain::D3Q1911_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_11Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1911_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1911_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1911_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1911_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1911_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1911_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f14(14);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1911_incompDomain::pdfDomain()
{
	return pdf_11Domain;
}

PdfMaker* 
D3Q1911_incompDomain::implDomain()
{
	return &impl_11;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1912_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1912_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1912_incompDomain::D3Q1912_incompDomain(CLbmCase* pCase) 
{	
	pdf_12Domain = maker(Impl(), pCase);
}
D3Q1912_incompDomain::D3Q1912_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_12Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1912_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1912_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1912_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1912_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1912_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1912_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f13(13);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1912_incompDomain::pdfDomain()
{
	return pdf_12Domain;
}

PdfMaker* 
D3Q1912_incompDomain::implDomain()
{
	return &impl_12;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1913_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1913_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1913_incompDomain::D3Q1913_incompDomain(CLbmCase* pCase) 
{	
	pdf_13Domain = maker(Impl(), pCase);
}
D3Q1913_incompDomain::D3Q1913_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_13Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1913_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1913_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1913_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1913_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1913_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1913_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f12(12);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1913_incompDomain::pdfDomain()
{
	return pdf_13Domain;
}

PdfMaker* 
D3Q1913_incompDomain::implDomain()
{
	return &impl_13;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1914_incompDomain concrete class
//Properties:
//1. Derived from D3qN15_incompDomain 
//2. Implementation of process specific to D3Q1914_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1914_incompDomain::D3Q1914_incompDomain(CLbmCase* pCase) 
{	
	pdf_14Domain = maker(Impl(), pCase);
}
D3Q1914_incompDomain::D3Q1914_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_14Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1914_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1914_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1914_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1914_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1914_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1914_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f11(11);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1914_incompDomain::pdfDomain()
{
	return pdf_14Domain;
}

PdfMaker* 
D3Q1914_incompDomain::implDomain()
{
	return &impl_14;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1915_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1915_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1915_incompDomain::D3Q1915_incompDomain(CLbmCase* pCase) 
{	
	pdf_15Domain = maker(Impl(), pCase);
}
D3Q1915_incompDomain::D3Q1915_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_15Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1915_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1915_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1915_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1915_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1915_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1915_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t  f18(18);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f18)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal);
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1915_incompDomain::pdfDomain()
{
	return pdf_15Domain;
}

PdfMaker* 
D3Q1915_incompDomain::implDomain()
{
	return &impl_15;
}

/*--------------------------------------------------------------------------------------------*/
//D3Q1916_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1916_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1916_incompDomain::D3Q1916_incompDomain(CLbmCase* pCase) 
{	
	pdf_16Domain = maker(Impl(), pCase);
}
D3Q1916_incompDomain::D3Q1916_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_16Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1916_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1916_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1916_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1916_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1916_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1916_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f17(17);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f17)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1916_incompDomain::pdfDomain()
{
	return pdf_16Domain;
}

PdfMaker* 
D3Q1916_incompDomain::implDomain()
{
	return &impl_16;
}


/*--------------------------------------------------------------------------------------------*/
//D3Q1917_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1917_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1917_incompDomain::D3Q1917_incompDomain(CLbmCase* pCase) 
{	
	pdf_17Domain = maker(Impl(), pCase);
}
D3Q1917_incompDomain::D3Q1917_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_17Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1917_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1917_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1917_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1917_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1917_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1917_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f16(16);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f16)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1917_incompDomain::pdfDomain()
{
	return pdf_17Domain;
}

PdfMaker* 
D3Q1917_incompDomain::implDomain()
{
	return &impl_17;
}


/*--------------------------------------------------------------------------------------------*/
//D3Q1918_incompDomain concrete class
//Properties:
//1. Derived from D3qN_incompDomain 
//2. Implementation of process specific to D3Q1918_incompDomain objects
/*--------------------------------------------------------------------------------------------*/
D3Q1918_incompDomain::D3Q1918_incompDomain(CLbmCase* pCase) 
{	
	pdf_18Domain = maker(Impl(), pCase);
}
D3Q1918_incompDomain::D3Q1918_incompDomain(CLbmCase* pCase, cgsize_t nvertex) 
{	
	pdf_18Domain = maker(Impl(), pCase, nvertex);
}
Pdf_bounder&
D3Q1918_incompDomain::GetFullBBBounder()
{
	return wall_bb;
}

Pdf_bounder&
D3Q1918_incompDomain::GetPeriodicBounder()
{
	return io_periodic;
}

Pdf_bounder&
D3Q1918_incompDomain::GetVelBounder()
{
	return dirich_boundary;
}
Pdf_bounder&
D3Q1918_incompDomain::GetOpenBounder()
{
	return open_bdry_f;
}

Pdf_setNeighbor&
D3Q1918_incompDomain::GetNeighbor()
{
	return neighbor;
}

void 
D3Q1918_incompDomain::BuriedBounder(Domain& domainVariables,cgsize_t node, cgsize_t dNqN)
{
	cgsize_t f15(15);
	const Node::NodeValueType_t c1(0.5);
	pdf()->GetPdfunction()->at(node-1)->m_nodeVal = c1 * (pdf()->GetPdfunction()->at(node-1)->m_nodeVal +
				domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f15)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal );
	pdf()->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
PdfBlock*  
D3Q1918_incompDomain::pdfDomain()
{
	return pdf_18Domain;
}

PdfMaker* 
D3Q1918_incompDomain::implDomain()
{
	return &impl_18;
}



//============================================================================================//
//FUNCTORS D3Q19
//base class
Node::NodeValueType_t 
D3Q19_bounder::ZhouHe_ortho(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t node, 
													Node::NodeValueType_t c, cgsize_t vel_index, LbmDomain::BOUND b_type, Node::LBMBOUND bound)
{
	/*cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);
	cgsize_t  rho(4), xVel(0), yVel(1), zVel(2);
	Node::NodeValueType_t density(0.0), velocity(0.0), temp(0.0);
	
	if     (b_type == LbmDomain::PRESSURE)
	{
		density = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		if(bound  == Node::WEST)
		{
			Node::NodeValueType_t c1 = -1.0;
			Node::NodeValueType_t c2 =  1.0;
			velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f3, f4, f2, f8, f10, f5, f6, f12, f14);
		}
		else if(bound == Node::EAST)
			{
				Node::NodeValueType_t c1 =  1.0;
				Node::NodeValueType_t c2 = -1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f3, f4, f1, f7, f9, f5, f6, f11, f13);
			}
		else if(bound == Node::SOUTH)
			{
				Node::NodeValueType_t c1 = -1.0;
				Node::NodeValueType_t c2 =  1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f4, f8, f10, f5, f6, f11, f13);
			}
		else if(bound == Node::NORTH)
			{
				Node::NodeValueType_t c1 =  1.0;
				Node::NodeValueType_t c2 = -1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f3, f7, f9, f5, f6, f12, f14);
			}
		else if(bound == Node::BOTTOM)
			{
				Node::NodeValueType_t c1 =  -1.0;
				Node::NodeValueType_t c2 =   1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f6, f8, f9, f3, f4, f12, f13);
			}
		else if(bound == Node::TOP)
			{
				Node::NodeValueType_t c1 =  1.0;
				Node::NodeValueType_t c2 = -1.0;
				velocity = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, c1, c2, f0, f1, f2, f5, f7, f10, f3, f4, f11, f14);
			}
	}

	else if(b_type == LbmDomain::VELOCITY)
	{
		velocity = domainVariables.GetDomainVariables()->at(vel_index)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;

		//now determine the density
		if(bound  == Node::WEST)
		{
			Node::NodeValueType_t c = -1.0;
			density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f2, f8, f10, f5, f6, f12, f14);
		}
		else if(bound == Node::EAST)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f1, f7, f9, f5, f6, f11, f13);
			}
		else if(bound == Node::SOUTH)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f4, f8, f10, f5, f6, f11, f13);
			}
		else if(bound == Node::NORTH)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f3, f7, f9, f5, f6, f12, f14);
			}
		else if(bound == Node::BOTTOM)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f6, f8, f9, f3, f4, f12, f13);
			}
		else if(bound == Node::TOP)
			{	
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f5, f7, f10, f3, f4, f11, f14);
			}
	}
	
	temp = domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_oldNodeVal + 
				 c*density*velocity;
	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"node = "<<node<<" pdf = "<< f << " Value = "<<temp <<"\n";
	std::cin.get();
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return temp;*/
	return 0;
}

Node::NodeValueType_t 
D3Q19_bounder::ZhouHe_diag(Domain& domainVariables, cgsize_t dNqN, cgsize_t f, cgsize_t f_opp, cgsize_t fI, cgsize_t fII, cgsize_t ffI, cgsize_t ffII, 
		cgsize_t node, Node::NodeValueType_t c1, Node::NodeValueType_t c2, Node::NodeValueType_t c3, Node::NodeValueType_t c4, Node::NodeValueType_t c5,
		cgsize_t xvel, cgsize_t yvel, cgsize_t zvel, LbmDomain::BOUND b_type, Node::LBMBOUND bound)
{
	/*cgsize_t rho(4);
	cgsize_t f0(0), f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8), f9(9), f10(10), f11(11), f12(12), f13(13), f14(14);
	cgsize_t xVel(0), yVel(1), zVel(2);
	Node::NodeValueType_t u_wall(0.0), v_wall(0.0), w_wall(0.0), density(0.0), temp(0.0);
	

	if     (b_type == LbmDomain::PRESSURE)
	{
		Node::NodeValueType_t a1 (0.0), a2(0.0);
		density = domainVariables.GetDomainVariables()->at(rho)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		if(  bound == Node::WEST)
			{
				a1=-1.0, a2=1.0;
				u_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f3, f4, f2, f8, f10, f5, f6, f12, f14);
			}
		else if(	bound == Node::EAST )
			{
				a1=1.0, a2=-1.0;
				u_wall  = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f3, f4, f1, f7, f9, f5, f6, f11, f13);
			}
		else if(bound == Node::SOUTH)
			{
				a1=-1.0, a2=1.0;
				v_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f1, f2, f4, f8, f10, f5, f6, f11, f13);
			}
		else if(bound == Node::NORTH )
			{
				a1=1.0, a2=-1.0;
				v_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1, a2, f0, f1, f2, f3, f7, f9, f5, f6, f12, f14);
			}
		else if(bound == Node::BOTTOM)
			{
				a1=-1.0, a2=1.0;
				w_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1,a2, f0, f1, f2, f6, f8, f9, f3, f4, f12, f13);
			}
		else if(bound == Node::TOP)
			{
				a1=1.0, a2=-1.0;
				w_wall = domainVariables.GetDomainMaterial()->at(dNqN)->velocity(domainVariables, dNqN, node, a1,a2, f0, f1, f2, f5, f7, f10, f3, f4, f11, f14);
			}
	}
	else if(b_type == LbmDomain::VELOCITY)
	{//velocities are specified therefore we need to determine density (rho) at boundary
		u_wall = domainVariables.GetDomainVariables()->at(xvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		v_wall = domainVariables.GetDomainVariables()->at(yvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		w_wall = domainVariables.GetDomainVariables()->at(zvel)->GetVariable()->at(dNqN)->at(node-1)->m_nodeVal;
		if(  bound == Node::WEST || bound == Node::WS  || bound == Node::WN  || bound == Node::WB || bound == Node::WT
			|| bound == Node::WSB  || bound == Node::WST || bound == Node::WNB || bound == Node::WNT)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f2, f8, f10, f5, f6, f12, f14);
			}

		else if(	bound == Node::EAST || bound == Node::ES  || bound == Node::EN  || bound == Node::EB || bound == Node::ET
			||			bound == Node::ESB  || bound == Node::EST || bound == Node::ENB || bound == Node::ENT)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, xVel, f0, f3, f4, f1, f7, f9, f5, f6, f11, f13);
			}

		else if(bound == Node::SOUTH || bound == Node::SB  || bound == Node::ST)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f4, f8, f10, f5, f6, f11, f13);
			}
		else if(bound == Node::NORTH || bound == Node::NB  || bound == Node::NT)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, yVel, f0, f1, f2, f3, f7, f9, f5, f6, f12, f14);
			}
		else if(bound == Node::BOTTOM)
			{
				Node::NodeValueType_t c = -1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f6, f8, f9, f3, f4, f12, f13);
			}
		else if(bound == Node::TOP)
			{
				Node::NodeValueType_t c = 1.0;
				density = domainVariables.GetDomainMaterial()->at(dNqN)->density(domainVariables, dNqN, node, c, zVel, f0, f1, f2, f5, f7, f10, f3, f4, f11, f14);
			}
	}

	temp =		(c1 * density *u_wall) + (c2 * density * v_wall) + (c3 * density * w_wall)
					+ domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(f_opp)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal

					+ c4 * ( domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(fI)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal
					       - domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(fII)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal )
								 
					+ c5 * ( domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(ffI)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal
					       - domainVariables.GetDomainMaterial()->at(dNqN)->GetLatticePdf()->at(ffII)->pdf()->GetPdfunction()->at(node-1)->m_nodeVal ); 	

	////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<"node = "<<node<<" pdf = "<< f << " Value = "<<temp <<"\n";
	std::cin.get();
	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
	return temp;*/
	return 0;
}

//derived classes
void
D3Q1900_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	//NULL BODY
}
void
D3Q1900_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}
void
D3Q1900_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1900_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D3Q1900_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//NULL BODY
}

void
D3Q1901_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f2(2);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f2)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1901_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	/*cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), xvel(0);
	Node::NodeValueType_t  c = (2.0/3.0);
	if(f_type == LbmDomain::ZHOU_HE)
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_ortho(domainVariables, d3q15, f1, f2, node, c, xvel, b_type, pblock->GetPdfunction()->at(node-1)->m_bound);

	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;*/
}
void
D3Q1901_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1901_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f1(1), f2(2);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f1, f2, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f1, f2, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f1, f2, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1901_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1902_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f1(1);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f1)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1902_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}

void
D3Q1902_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1902_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f1(1), f2(2);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f2, f1, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f2, f1, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f2, f1, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1902_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1903_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f4(4);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f4)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1903_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1903_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1903_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f3(3), f4(4);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f3, f4, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f3, f4, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f3, f4, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1903_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1904_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f3(3);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f3)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1904_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1904_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1904_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f4(4), f3(3);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f4, f3, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f4, f3, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f4, f3, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1904_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1905_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f6(6);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f6)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1905_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1905_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1905_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f5(5), f6(6);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f5, f6, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f5, f6, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f5, f6, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1905_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1906_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f5(5);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f5)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1906_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1906_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1906_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f6(6), f5(5);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f6, f5, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfBB = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f6, f5, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f6, f5, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1906_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1907_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f10(10);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f10)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1907_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	/*cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WEST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WS
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WB
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::WSB)
			{
				c1 = (1.0/12.0), c2 = (1.0/4.0), c3 = (1.0/4.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f7, f8, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SOUTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::SB)
			{
				c1 = (1.0/4.0), c2 = (1.0/12.0), c3 = (1.0/4.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f7, f8, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::BOTTOM)
			{
				c1 = (1.0/4.0), c2 = (1.0/4.0), c3 = (1.0/12.0), c4 = (-1.0/4.0), c5 = (-1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f7, f8, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;*/
}
void
D3Q1907_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1907_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f7(7), f10(10);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f7, f10, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f7, f10, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f7, f10, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1907_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1908_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f9(9);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f9)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}		
}
void
D3Q1908_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	/*cgsize_t xVel(0), yVel(1),  zVel(2);
	cgsize_t d3q15(mat_index);
	cgsize_t f1(1), f2(2), f3(3), f4(4), f5(5), f6(6), f7(7), f8(8);
	Node::NodeValueType_t c1(0.0), c2(0.0), c3(0.0), c4(0.0), c5(0.0);
	//for ZHOU_HE
	if(f_type == LbmDomain::ZHOU_HE)
	{
			if(	pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EAST
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::EN
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ET
				||pblock->GetPdfunction()->at(node-1)->m_bound    == Node::ENT)
			{
				c1 = (-1.0/12.0), c2 = (-1.0/4.0), c3 = (-1.0/4.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f8, f7, f3,f4,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if( pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NORTH
						|| pblock->GetPdfunction()->at(node-1)->m_bound    == Node::NT)
			{
				c1 = (-1.0/4.0), c2 = (-1.0/12.0), c3 = (-1.0/4.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f8, f7, f1,f2,f5,f6,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}

			else if(pblock->GetPdfunction()->at(node-1)->m_bound    == Node::TOP)
			{
				c1 = (-1.0/4.0), c2 = (-1.0/4.0), c3 = (-1.0/12.0), c4 = (1.0/4.0), c5 = (1.0/4.0);
				pblock->GetPdfunction()->at(node-1)->m_nodeVal = ZhouHe_diag(domainVariables,d3q15,f8, f7, f1,f2,f3,f4,node, 
					c1,c2,c3,c4,c5,xVel,yVel,zVel,b_type,pblock->GetPdfunction()->at(node-1)->m_bound);
			}
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;*/
}
void
D3Q1908_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1908_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f8(8), f9(9);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f8, f9, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f8, f9, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f8, f9, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1908_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1909_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f8(8);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f8)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1909_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1909_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1909_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f9(9), f8(8);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f9, f8, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f9, f8, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f9, f8, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1909_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1910_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f7(7);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f7)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1910_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1910_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1910_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f10(10), f7(7);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f10, f7, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f10, f7, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f10, f7, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1910_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1911_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f14(14);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f14)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1911_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1911_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1911_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f11(11), f14(14);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f11, f14, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f11, f14, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f11, f14, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1911_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
void
D3Q1912_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f13(13);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f13)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}
void
D3Q1912_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1912_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1912_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f12(12), f13(13);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f12, f13, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f12, f13, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f12, f13, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1912_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1913_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f12(12);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f12)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1913_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1913_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}

void
D3Q1913_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f13(13), f12(12);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f13, f12, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f13, f12, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f13, f12, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1913_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1914_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f11(11);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f11)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1914_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1914_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}


void
D3Q1914_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f14(14), f11(11);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f14, f11, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f14, f11, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f14, f11, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1914_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1915_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f18(18);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f18)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1915_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1915_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}


void
D3Q1915_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f15(15), f18(18);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f15, f18, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f15, f18, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f15, f18, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1915_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1916_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f17(17);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f17)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1916_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1916_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}


void
D3Q1916_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f16(16), f17(17);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f16, f17, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f16, f17, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f16, f17, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1916_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1917_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f16(16);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f16)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1917_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1917_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}


void
D3Q1917_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f17(17), f16(16);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f17, f16, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f17, f16, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f17, f16, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1917_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1918_setNeighbor::operator() (PdfBlock* pblock, PdfDomain* allFs, cgsize_t mat_index)
{
	cgsize_t f15(15);

	for(PdfBlock::NodeType_t::size_type node = 1; node <= pblock->GetPdfunction()->size(); node++)
	{
		pblock->GetPdfunction()->at(node-1)->m_neighborNum = allFs->GetLatticePdf()->at(f15)->pdf()->GetPdfunction()->at(node-1)->m_connectTo;
	}
}

void
D3Q1918_Dirichlet::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	//place holder
}
void
D3Q1918_openBound::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND   b_type, LbmDomain::FORMULA f_type)
{
	
}


void
D3Q1918_BBack::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t d3q19(mat_index);
	cgsize_t f18(18), f15(15);
	if(f_type == LbmDomain::HALF_BB)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = halfWay_BB(domainVariables, d3q19, f18, f15, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node halfbb = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::MOVING)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = moving_BB(domainVariables, d3q19, f18, f15, node);
		/*///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<"node moving = "<<node<<" value = "<<pblock->GetPdfunction()->at(node-1)->m_nodeVal<<"\n";
		std::cin.get();
		//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	}
	else if(f_type == LbmDomain::EXTRAPOLATE)
	{
		pblock->GetPdfunction()->at(node-1)->m_nodeVal = extrapolation_BB(domainVariables, d3q19, f18, f15, node);
	}
	else;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}

void
D3Q1918_Periodic::operator() (cgsize_t mat_index, PdfBlock* pblock, Domain& domainVariables, cgsize_t node, LbmDomain::BOUND b_type, LbmDomain::FORMULA f_type)
{
	cgsize_t neighbor(0), useNode(0);
	neighbor = pblock->GetPdfunction()->at(node-1)->m_neighborNum;

	while(neighbor)
	{
		useNode  = neighbor;
		neighbor = pblock->GetPdfunction()->at(useNode-1)->m_neighborNum;
	};
	pblock->GetPdfunction()->at(node-1)->m_nodeVal = pblock->GetPdfunction()->at(useNode-1)->m_nodeVal;
	pblock->GetPdfunction()->at(node-1)->m_valueUpdated = true;
}
