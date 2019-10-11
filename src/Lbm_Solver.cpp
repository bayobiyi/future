#include <iostream>
#include <iomanip>
#include <time.h>
#include <cmath>
#include <algorithm>
#include "CLbmCase.h"
#include "CLbmSolve.h"
#include "CLbmSolve3D.h"
#include <sstream>
#include <string.h>
#include <omp.h>


void Process(Domain* pdomain, CLbmCase* pCase, cgsize_t iteration);
void outputLine (Convergence* pconvergence, cgsize_t iteration,double time/*clock_t time*/);
void InitializeCase(CLbmCase* pCase);
void writeData(CLbmCase* pCase, char* transfername, char* dataname, cgsize_t iteration, cgsize_t save_iteration);
void setDataName(CLbmCase* pCase, char* dataname, cgsize_t iteration, cgsize_t save_iteration);

int 
main(int argc, char** argv)
{
	double runStart(0), currentRun(0);
	double save_time = 82800;
	runStart = omp_get_wtime();
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" starting";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	CLbmCase* pCase = new CLbmCase;
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" about to import case information \n";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	pCase->ImportCaseInfo(argv[1]);
	InitializeCase(pCase);

	//try to patch wettability (turn on/off)
	//cgsize_t rho(4);
	//pCase->m_pDomain->GetDomainVariables()->at(rho)->setPatchedNodesWettability(pCase);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" done with wetting patch";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	//try to patch solid object (turn on/off) if object is not moving	
	//pCase->m_pDomain->setPatchedSolidObjects(pCase);

	//try to patch micropillars	
	//pCase->m_pDomain->setPatchedMicroPillars(pCase);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
	std::cout<<" done with micropillar patch";
	std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/

	//set number of threads
	omp_set_num_threads(22);
	//std::cout<<"Number of processes = "<<omp_get_num_procs()<<"\n";
	pCase->m_pDomain->SetNumThreads(omp_get_num_threads());
	//std::cout<<"Number of threads = "<<pCase->m_pDomain->GetNumThreads()<<"\n";
	//std::cin.get();

	cgsize_t iteration = pCase->m_sol.m_iteration;
	cgsize_t save_iteration = atoi(argv[5]);
	std::cout<<std::setiosflags(std::ios::left)<<std::setw(10)<<"Iteration"
		<< std::setw(13) << "x-velocity"<<std::setprecision(7)
		<< std::setw(13) << "y-velocity"<<std::setprecision(7)
		<< std::setw(13) << "z-velocity"<<std::setprecision(7)
		<< std::setw(13) << "time/iter"<<std::setprecision(7)
		<<'\n';

	pCase->m_pDomain->GetConvergence()->x_residual = pCase->m_sol.m_convergenceData.x_residual;
	pCase->m_pDomain->GetConvergence()->y_residual = pCase->m_sol.m_convergenceData.y_residual;
	pCase->m_pDomain->GetConvergence()->z_residual = pCase->m_sol.m_convergenceData.z_residual;

 	while(      iteration < atoi(argv[2]) && save_time > (currentRun - runStart) 
		&& 
	           (  pCase->m_pDomain->GetConvergence()->x_residual > 0.000001
	           || pCase->m_pDomain->GetConvergence()->y_residual > 0.000001
	           || pCase->m_pDomain->GetConvergence()->z_residual > 0.000001
		   )
	    )
	{
		iteration += 1;
		Process(pCase->m_pDomain, pCase, iteration);
		if(iteration%save_iteration == 0)
			writeData(pCase,argv[3], argv[4], iteration, save_iteration);
		else;
		currentRun=omp_get_wtime();	
	}
	writeData(pCase,argv[3], argv[4], iteration, save_iteration);
	return(0);
}

void 
Process(Domain* pdomain, CLbmCase* pCase, cgsize_t iteration)
{
	//double starttime;
	//double endtime;
	double /*clock_t*/ time;
	//time = clock();
	time = omp_get_wtime();
	pCase->m_sol.m_iteration  = iteration;
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now start collision \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//starttime = omp_get_wtime();
	pdomain->collide(pCase, *pdomain);
	/*/////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now start post collision \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//endtime = omp_get_wtime();
	//printf("Collision took %f sec. time.\n", endtime-starttime);
	//starttime = omp_get_wtime();
	pdomain->post_collide(pCase, *pdomain);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now stream nodes \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//endtime = omp_get_wtime();
	//printf("Post collision took %f sec. time.\n", endtime-starttime);
	//starttime = omp_get_wtime();
	pdomain->stream (pCase, *pdomain);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now start post stream nodes \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//endtime = omp_get_wtime();
	//printf("Streaming took %f sec. time.\n", endtime-starttime);
	//starttime = omp_get_wtime();
	pdomain->post_stream (pCase, *pdomain);
	/*//////////////////////////////////////////////////////////////////////////////////////////////////////////
		std::cout<<" now update variables \n";
		std::cin.get();
	////////////////////////////////////////////////////////////////////////////////////////////////////////////*/
	//endtime = omp_get_wtime();
	//printf("Post streaming took %f sec. time.\n", endtime-starttime);
	//starttime = omp_get_wtime();
	//MCMP
	//pdomain->updateVariables(pCase, *pdomain);
	//SCMP
	pdomain->updateVariablesSglComp(pCase, *pdomain);
 	//endtime = omp_get_wtime();
	//printf("Updating took %f sec. time.\n", endtime-starttime);
	//starttime = omp_get_wtime();
	if(iteration%1 == 0)
	{
		pdomain->convergence(*pdomain, pCase);
		//time   = clock() - time;
		time = omp_get_wtime() - time;
		//endtime = omp_get_wtime();
		//printf("Convergence took %f sec. time.\n", endtime-starttime);
		outputLine (pdomain->GetConvergence(), iteration, time);
	}
	else;
	pCase->m_sol.m_convergenceData.x_residual = pdomain->GetConvergence()->x_residual;
	pCase->m_sol.m_convergenceData.y_residual = pdomain->GetConvergence()->y_residual;
	pCase->m_sol.m_convergenceData.z_residual = pdomain->GetConvergence()->z_residual;
}

void 
outputLine (Convergence* pconvergence, cgsize_t iteration, double time /*clock_t time*/)
{
	std::cout<<std::setiosflags(std::ios::left)<<std::setw(10)<<iteration
		<< std::setw(13) << pconvergence->x_residual<<std::setprecision(7)
		<< std::setw(13) << pconvergence->y_residual<<std::setprecision(7)
		<< std::setw(13) << pconvergence->z_residual<<std::setprecision(7)
		<< std::setw(13) << /*((float)time)/CLOCKS_PER_SEC*/time<<std::setprecision(7)
		<<'\n';
}

void
InitializeCase(CLbmCase* pCase)
{
	std::cout<<" Starting  initialization...\n";
	CLbmModel::SOLVER input = pCase->m_model.m_solver;		
	switch (input)
	{
		case CLbmModel::D2Q9:
			{
				pCase->m_pDomain =  new Domain_D2Q9(pCase);
			}
		break;
		case CLbmModel::D3Q15: 
			{
				pCase->m_pDomain =  new Domain_D3Q15(pCase);
			}
		break;
		case CLbmModel::D3Q19: 
			{
				pCase->m_pDomain =  new Domain_D3Q19(pCase);
			}
		break;
		default:
		break;
	}
	
	if(pCase->m_model.m_initStatus == 'i')
	{
		//update variables before initializing if you are continuing iterations
		//pCase->m_pDomain->Lbm_read(pCase->m_importfs);
		pCase->m_pDomain->Update_variables(pCase);
	}
	else;
	pCase->m_pDomain->initialize(pCase,*pCase->m_pDomain); 
	std::cout<<" Completed initialization...\n";
}

void 
writeData(CLbmCase* pCase, char* transfername, char* dataname, cgsize_t iteration, cgsize_t save_iteration)
{
	cgsize_t rho(4);
	//write solution from the different containers into an array
	pCase->m_pDomain->writeSolution(pCase);
	pCase->m_model.m_initStatus = 'i';
	setDataName(pCase, dataname, iteration, save_iteration);
	pCase->ExportCaseInfo(transfername);
	Node::NodeValueType_t sum_rho(0.0);
	//I will be outputing the sum of rho for each component temporarily here
	std::cout<<" *------------------------------SUMMARY--------------------------------------------------------*\n";
	for(CVariable::MaterialVarType_t::size_type j=0; j<pCase->m_pDomain->GetDomainVariables()->at(rho)->GetVariable()->size(); j++)
	{
		sum_rho = 0.0;
		for(CVariable::NodeType_t::size_type k=0; k < pCase->m_pDomain->GetDomainVariables()->at(rho)->GetVariable()->at(j)->size(); k++)
		{
			sum_rho += pCase->m_sol.m_sol.at(rho)->at(j)->at(k);	
		}
		std::cout<<" Average density for component "<<j<<" = "<<sum_rho/(pCase->m_pDomain->GetDomainVariables()->at(rho)->GetVariable()->at(j)->size())<<"\n";
	}
	
	std::cout<<" *------------------------------SUMMARY--------------------------------------------------------*\n";
}

void
setDataName(CLbmCase* pCase, char* dataname, cgsize_t iteration, cgsize_t save_iteration)
{
	static cgsize_t prev_iteration = 0;
	std::string iteration_str = static_cast<std::ostringstream*>( &(std::ostringstream() << iteration) )->str();
	std::string prev_iteration_str = static_cast<std::ostringstream*>( &(std::ostringstream() << prev_iteration) )->str();
	cgsize_t iteration_length = strlen(iteration_str.c_str());
	cgsize_t prev_iteration_length = strlen(prev_iteration_str.c_str());
	char* Dpch1 ;
	Dpch1 = strstr (dataname,".dat");
	if(prev_iteration == 0)
		prev_iteration_length -= 1;
	else;

	strncpy (Dpch1-prev_iteration_length, iteration_str.c_str(), iteration_length);	
	strncpy (Dpch1-prev_iteration_length+iteration_length,".dat",4);
	strncpy (Dpch1-prev_iteration_length+iteration_length+4,"\0",1);
	pCase->Lbm_write_tecplot(dataname);
	prev_iteration = iteration;
}
