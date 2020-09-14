/*
 *  BFBDS.h
 *
 *  Created by Steven Danishevski on SEP 20.
 *
 */

#ifndef GenericTester_H
#define GenericTester_H

#include <climits>
#include <algorithm>
#include "TemplateAStar.h"
#include "AStarOpenClosed.h"
#include "GenericTester.h"
#include "IDAStar.h"
#include "MM.h"
#include "MBBDS.h"
#include "BFBDS.h"
#include "IDBiHS.h"
#include "IDTHSwTrans.h"
#include "MbbdsBloomFilter.h"
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <math.h> 
#include <ctime>
using namespace std;

template <class state, class action, class environment, class hasher>
class GenericTester {
private:
	unsigned long MMstatesQuantityBound;
	unsigned long ASTARstatesQuantityBound;
	unsigned long statesQuantityBound;
	unsigned long statesQuantityBoundDefault = 1000000;
	int secondsLimit = 60*30;

	bool AstarRun=true;
	bool RevAstarRun=true;

	bool IDAstarRun=true;

	bool AstarPIDAstarRun=true;
	bool AstarPIDAstarReverseRun=true;
	bool AstarPIDAstarReverseMinHRun=true;
	bool IDTHSpTrans = true;

	bool BAI=true;
	bool Max_BAI=true;

	bool MMRun=true;

	bool IDBiHSRun=true;
	bool F2Fheuristics=true;

	bool ASTARpIDBiHS=true;
	bool MMpIDBiHS=false;

	bool MBBDSRun=true;
	bool BFBDSRUN=true;
	bool fullMBBDS=true;
	bool revAlgo=false;

	bool threePhase=true;
	bool twoPhase=false;

	bool detectDuplicate=true;
	bool isConsistent=true;
	bool isUpdateByWorkload=true;
	bool printAbstractAlgos = false;
	
public:
	GenericTester() {}
	virtual ~GenericTester() {}
	void genericTest(state original, state goal, environment env, ofstream &myfile);
};

template <class state, class action, class environment, class hasher>
void GenericTester<state, action, environment, hasher>::genericTest(state original, state goal, environment env, ofstream &myfile){
	myfile << "\tStart state: " << original << endl;
	myfile << "\tGoal state: " << goal << endl;
	myfile <<"\tInitial heuristic " << env.HCost(original, goal) << endl;
	Timer timer;
	state start;
	vector<state> astarPath;
	vector<state> mmPath;
	vector<state> fullMbbdsPath;
	vector<state> idaPath;
	state midState;
	statesQuantityBound = ULONG_MAX;
	// A*
	vector<AStarOpenClosedDataWithF<state>> astarOpenList;
	if (AstarRun)
	{
		myfile <<"\t\t_A*_\n";
		TemplateAStar<state, action, environment> astar;
		start = original;
		timer.StartTimer();
		bool solved = astar.GetPathTime(&env, start, goal, astarPath, secondsLimit);
		ASTARstatesQuantityBound = astar.getMemoryStatesUse();
		timer.EndTimer();
		if(solved){
			myfile << boost::format("\t\t\tA* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % env.GetPathLength(astarPath) %
			   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % timer.GetElapsedTime();
			if(printAbstractAlgos)
				myfile << boost::format("\t\t\tI-A* ; %llu expanded;\n") % astar.getIAstarExpansions();	
			statesQuantityBound =  std::min(statesQuantityBound, ASTARstatesQuantityBound);
		}
		else{
			myfile << boost::format("\t\t\tA* failed after %1.4fs\n") % timer.GetElapsedTime();
			myfile << "\t\t\tI-A* failed after because A* failed\n";	
		}
	}
	if (RevAstarRun){
		myfile <<"\t\t_Rev-A*_\n";
		TemplateAStar<state, action, environment> astar;
		start = original;
		timer.StartTimer();
		bool solved = astar.GetPathTime(&env, goal, start, astarPath, secondsLimit);
		ASTARstatesQuantityBound = astar.getMemoryStatesUse();
		timer.EndTimer();
		if(solved){
			myfile << boost::format("\t\t\tRev-A* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % env.GetPathLength(astarPath) %
			   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % timer.GetElapsedTime();
			statesQuantityBound =  std::min(statesQuantityBound, ASTARstatesQuantityBound);
		}
		else{
			myfile << boost::format("\t\t\tRev-A* failed after %1.4fs\n") % timer.GetElapsedTime();

		}
	}
	// MM
	if (MMRun){
		myfile << "\t\t_MM_\n";				
		MM<state, action, environment> mm;
		//goal.Reset();
		start = original;
		timer.StartTimer();
		bool solved = mm.GetPath(&env, start, goal, &env, &env, mmPath, secondsLimit);
		MMstatesQuantityBound = mm.getMemoryStatesUse();
		timer.EndTimer();
		if(solved){
			myfile << boost::format("\t\t\tMM found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % env.GetPathLength(mmPath) % mm.GetNodesExpanded() % mm.GetNecessaryExpansions() % MMstatesQuantityBound % timer.GetElapsedTime();
			if(printAbstractAlgos)
				myfile << boost::format("\t\t\tI-MM ; %llu expanded;\n") % mm.getIMMExpansions();
			statesQuantityBound =  std::min(statesQuantityBound, MMstatesQuantityBound);
		}
		else{
			myfile << boost::format("\t\t\tMM failed after %1.4fs\n") % timer.GetElapsedTime();
			if(printAbstractAlgos)
				myfile << "\t\t\tI-MM failed because MM failed\n";
		}
	}
	if(statesQuantityBound == ULONG_MAX){
		statesQuantityBound = statesQuantityBoundDefault;
	}
	// IDA*
	if (IDAstarRun)
	{
		myfile << "\t\t_IDA*_\n";
		IDAStar<state, action, false> idastar;
		//goal.Reset();
		start = original;
		timer.StartTimer();
		bool solved = idastar.GetPath(&env, start, goal, idaPath, secondsLimit);
		timer.EndTimer();
		if(solved){
			myfile << boost::format("\t\t\tIDA* found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % env.GetPathLength(idaPath) %
			   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % timer.GetElapsedTime();
			if(printAbstractAlgos)
				myfile << boost::format("\t\t\tD-A* ; %llu expanded;\n") % idastar.getDAstarExpansions();
		}
		else{
			myfile << boost::format("\t\t\tIDA* failed after %1.4fs\n") % timer.GetElapsedTime();
			if(printAbstractAlgos)
				myfile << "\t\t\tD-A* failed because IDA* failed\n";
		}					   
	}

	double percentages[3] = {0.5, 0.1, 0.01};
	long stateSize = sizeof(original);
	//FullMBBDS
	if (BFBDSRUN){
		bool threePhase = false;
		myfile << "\t\t_BFBDS_\n";
		long stateSize = sizeof(original);
		bool solved;
		unsigned long nodesExpanded;
		for(double percentage : percentages){
			timer.StartTimer();
			unsigned long statesQuantityBoundforMBBDS = std::max(statesQuantityBound*percentage,2.0);;
			BFBDS<state, action, environment, MbbdsBloomFilter<state, hasher>, false> bfbds(statesQuantityBoundforMBBDS, isUpdateByWorkload, isConsistent, revAlgo, F2Fheuristics) ;
			bool threePhase = true;
			solved = bfbds.GetMidState(&env, start, goal, midState, fullMbbdsPath, secondsLimit, threePhase);
			if(solved){
				myfile << boost::format("\t\t\tBFBDS(k=1, ThreePhase=%d) using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(bfbds.isThreePhase()) % statesQuantityBoundforMBBDS % stateSize % percentage % bfbds.getPathLength() % bfbds.getNodesExpanded() % bfbds.getNecessaryExpansions() % bfbds.getIterationsNum() % timer.GetElapsedTime();
			}
			else{
				myfile << boost::format("\t\t\tBFBDS(k=1, ThreePhase=%d) using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % int(bfbds.isThreePhase()) % statesQuantityBoundforMBBDS % stateSize % percentage % timer.GetElapsedTime();
				break;
			} 
		}
	}
	/*if(MMpIDBiHS){
		myfile << "\t\t_MM+IDBiHS_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforMMpIDBiHS = std::max(statesQuantityBound*percentage,2.0);;
			MM<state, action, environment> mm;
			//goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = mm.GetPath(&env, start, goal, &env, &env, mmPath, secondsLimit, statesQuantityBoundforMMpIDBiHS);
			unsigned long nodesExpanded = mm.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			if(solved){
				timer.EndTimer();
				necessaryNodesExpanded = mm.GetNecessaryExpansions();
				myfile << boost::format("\t\t\tMM+IDBiHS MM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDBiHS % stateSize % percentage % env.GetPathLength(mmPath) %
				   nodesExpanded % necessaryNodesExpanded  % timer.GetElapsedTime();
			}
			else{
				IDBiHS<environment, state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload);
				state midState;
				bool solved = idbihs.GetMidStateFromLists(&env, start, goal, midState, secondsLimit-timer.GetElapsedTime(), mm.getLastBound(), mm.GetForwardItems(), mm.GetBackwardItems());
				nodesExpanded += idbihs.GetNodesExpanded();
				necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tMM+IDBiHS IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDBiHS % stateSize % percentage % idbihs.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tMM+IDBiHS IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforMMpIDBiHS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}*/
	if(ASTARpIDBiHS){
		myfile << "\t\t_A*+IDBiHS_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound*percentage, 2.0);;
			IDBiHS<environment, state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload);
			state midState;
			bool solved = idbihs.Astar_plus_IDBiHS(&env, start, goal, midState, statesQuantityBoundforASPIDBiHS, secondsLimit-timer.GetElapsedTime(), detectDuplicate);
			timer.EndTimer();
			if(solved){
			  myfile << boost::format("\t\t\tA*+IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % idbihs.getPathLength() % idbihs.GetNodesExpanded() % idbihs.GetNecessaryExpansions() % timer.GetElapsedTime();
			}
			else{
			  myfile << boost::format("\t\t\tA*+IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % timer.GetElapsedTime();
			  break;
			}	
		}
	}
	if(IDTHSpTrans){
		{
			myfile << "\t\t_IDTHSpTrans_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound*percentage, 2.0);
				//goal.Reset();
				start = original;
				timer.StartTimer();
				unsigned long nodesExpanded = 0;
				unsigned long necessaryNodesExpanded = 0;
				IDTHSwTrans<state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload, 1 , true);
				state midState;
				bool solved = idbihs.GetPath(&env, start, goal, /*astarPath, */secondsLimit, statesQuantityBoundforASPIDBiHS);
				nodesExpanded += idbihs.GetNodesExpanded();
				necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tIDTHSpTrans IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % idbihs.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tIDTHSpTrans IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
		{
			myfile << "\t\t_IDTHSpTrans_NDD_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound*percentage,2.0);;
				//goal.Reset();
				start = original;
				timer.StartTimer();
				unsigned long nodesExpanded = 0;
				unsigned long necessaryNodesExpanded = 0;
				IDTHSwTrans<state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload, 1 , false);
				state midState;
				bool solved = idbihs.GetPath(&env, start, goal, /*astarPath, */secondsLimit, statesQuantityBoundforASPIDBiHS);
				nodesExpanded += idbihs.GetNodesExpanded();
				necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tIDTHSpTrans_NDD IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % idbihs.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tIDTHSpTrans_NDD IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}
	if (AstarPIDAstarRun){
		myfile << "\t\t_Astar+IDAstar_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforASPIDAS = std::max(statesQuantityBound*percentage,2.0);;
			TemplateAStar<state, action, environment> astar;
			//goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDAS);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			if(solved){
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				myfile << boost::format("\t\t\tA*+IDA* A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % env.GetPathLength(astarPath) %
				   nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
			}
			else{
				IDAStar<state, action, false> idastar;
				solved = idastar.ASpIDA(&env, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime(), detectDuplicate);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tA*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % idastar.getSolLength() % nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tA*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}
	/*if (AstarPIDAstarRun){
		myfile << "\t\t_Astar+IDAstar_\n";
		for(double percentage : percentages){
			start = original;
			unsigned long statesQuantityBoundforASPIDAS = std::max(statesQuantityBound*percentage,2.0);;
			IDAStar<state, action, false, environment> idastar;
			bool solved = idastar.Astar_plus_IDAstar(&env, start, goal, idaPath, statesQuantityBoundforASPIDAS, secondsLimit-timer.GetElapsedTime(), detectDuplicate);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tA*+IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % idastar.getSolLength() % idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % timer.GetElapsedTime();
			}
			else{
				myfile << boost::format("\t\t\tA*+IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % timer.GetElapsedTime();
				break;
			}	
		}
	}*/
	if (AstarPIDAstarReverseRun){
		myfile << "\t\t_Astar+IDAstar+Reverse_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);
			TemplateAStar<state, action, environment> astar;
			//goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS, false);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			if(solved){
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				myfile << boost::format("\t\t\tA*+IDA*_Reverse A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % env.GetPathLength(astarPath) %
				   nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
			}
			else{
				IDAStar<state, action, false> idastar;
				solved = idastar.ASpIDArev(&env, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(),secondsLimit-timer.GetElapsedTime());
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tA*+IDA*_Reverse IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
					   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tA*+IDA*_Reverse IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}
	if (AstarPIDAstarReverseMinHRun){
		myfile << "\t\t_Astar+IDAstar+Reverse+MinH_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);;
			TemplateAStar<state, action, environment> astar;
			//goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS, false);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			if(solved){
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				myfile << boost::format("\t\t\tA*+IDA*_Reverse+MinH A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % env.GetPathLength(astarPath) % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
			}
			else{
				IDAStar<state, action, false> idastar;
				solved = idastar.ASpIDArev(&env, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(),secondsLimit-timer.GetElapsedTime(),isConsistent, true);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tA*+IDA*_Reverse+MinH IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() % nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tA*+IDA*_Reverse+MinH IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}
	if (BAI){
		myfile << "\t\t_BAI_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);;
			TemplateAStar<state, action, environment> astar;
			//goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			if(solved){
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				myfile << boost::format("\t\t\tBAI A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % env.GetPathLength(astarPath) % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
			}
			else{
				IDAStar<state, action, false> idastar;
				solved = idastar.BAI(&env, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime(), false);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tBAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() % nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tBAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}
	if (Max_BAI){
		myfile << "\t\t_Max_BAI_\n";
		for(double percentage : percentages){
			unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);;
			TemplateAStar<state, action, environment> astar;
			//goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&env, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS);
			unsigned long nodesExpanded = astar.GetNodesExpanded();
			unsigned long necessaryNodesExpanded = 0;
			if(solved){
				timer.EndTimer();
				necessaryNodesExpanded = astar.GetNecessaryExpansions();
				myfile << boost::format("\t\t\tMax_BAI A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % env.GetPathLength(astarPath) % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
			}
			else{
				IDAStar<state, action, false> idastar;
				solved = idastar.BAI(&env, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime(), true);
				nodesExpanded += idastar.GetNodesExpanded();
				necessaryNodesExpanded += idastar.GetNecessaryExpansions();
				timer.EndTimer();
				if(solved){
					myfile << boost::format("\t\t\tMax_BAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() % nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tMax_BAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % timer.GetElapsedTime();
					break;
				}	
			}
		}
	}
	//IDBiHS
	if(IDBiHSRun){
		myfile << "\t\t_IDBiHS_\n";
		IDBiHS<environment, state, action, false> idbihs(F2Fheuristics, isConsistent, isUpdateByWorkload);
		//goal.Reset();
		start = original;
		state midState;
		timer.StartTimer();
		bool solved = idbihs.GetMidState(&env, start, goal, midState, secondsLimit);
		timer.EndTimer();
		if(solved){
			myfile << boost::format("\t\t\tIDBiHS found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed; ") % idbihs.getPathLength() % idbihs.GetNodesExpanded() % idbihs.GetNodesTouched() % idbihs.GetNecessaryExpansions() % timer.GetElapsedTime();
			myfile << "Mid state: " << midState << endl;
			if(printAbstractAlgos)
				myfile << boost::format("\t\t\tD-MM ; %llu expanded;\n") % idbihs.getDMMExpansions();
		}
		else{
			myfile << boost::format("\t\t\tIDBiHS failed after %1.4fs\n") % timer.GetElapsedTime();
			if(printAbstractAlgos)
				myfile << "\t\t\tD-MM failed because IDBiHS failed\n";
		}   				
	}	
}

#endif
