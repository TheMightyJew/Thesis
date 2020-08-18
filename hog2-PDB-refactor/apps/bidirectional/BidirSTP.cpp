//
//  BidirSTP.cpp
//  hog2 glut
//
//  Created by Nathan Sturtevant on 2/7/17.
//  Copyright © 2017 University of Denver. All rights reserved.
//
#include "BidirSTP.h"
#include "MNPuzzle.h"
#include "WeightedVertexGraph.h"
#include "STPInstances.h"
#include "LexPermutationPDB.h"
#include "MR1PermutationPDB.h"

#include <algorithm>
#include "STPHasher.h"
#include "TemplateAStar.h"
#include "AStarOpenClosed.h"
#include "IDAStar.h"
#include "MM.h"
#include "MBBDS.h"
#include "FullMBBDS.h"
#include "IDMM.h"
#include "MbbdsBloomFilter.h"
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <math.h> 
#include <ctime>
using namespace std;


static void StevenTest(int problems_num=1, bool randomSTP=true, vector<int> skipVector = vector<int>());

static int walkLength = 120;
static unsigned long MMstatesQuantityBound;
static unsigned long ASTARstatesQuantityBound;
static unsigned long statesQuantityBound = 1000000;
static int secondsLimit = 60*30;
static bool AstarRun=true;
static bool AstarPIDAstarRun=true;
static bool AstarPIDAstarReverseRun=false;
static bool ASTARpIDMM=true;
static bool MMRun=true;
static bool MMpIDMM=false;
static bool IDAstarRun=true;
static bool MBBDSRun=true;
static bool threePhase=true;
static bool twoPhase=false;
static bool IDMMRun=true;
static bool idmmF2fFlag=true;
static bool isConsistent=false;


static string datetime()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,80,"%d-%m-%Y_%H-%M-%S",timeinfo);
    return string(buffer);
}

static ofstream myfile;
static string filename = "Test_Results/STP/results_" + datetime() + ".txt";


void TestSTP(int algorithm)
{
	cout << "running..." << endl;
	myfile.open (filename);
	
	StevenTest(100, false);

	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
	exit(0);
}

void StevenTest(int problems_num, bool randomSTP, vector<int> skipVector)
{
	srandom(2017218);
	MNPuzzleState<4, 4> original, start, goal;
	MNPuzzle<4, 4> mnp;
	MNPuzzle<4, 4> mnp2;
	
	vector<MNPuzzleState<4, 4>> astarPath;
	vector<MNPuzzleState<4, 4>> mmPath;
	vector<MNPuzzleState<4, 4>> fullMbbdsPath;
	vector<MNPuzzleState<4, 4>> idaPath;
	MNPuzzleState<4, 4> midState;
	
	Timer timer;
	myfile << boost::format("TestSTP:(Random: %d, Hard: %d)\n") % randomSTP % (!randomSTP);
	for (int count = 0; count < problems_num; count++)
	{
		goal.Reset();
		original.Reset();
		if(randomSTP){
			original = STP::GetRandomInstance(walkLength);
		}
		else{
			original = STP::GetKorfInstance(count);
		}
		if(std::find(skipVector.begin(), skipVector.end(), count+1) != skipVector.end()) {
			continue;
		}
		myfile << boost::format("\tProblem %d of %d\n") % (count+1) % problems_num;
		myfile << "\tStart state: " << original << endl;
		myfile << "\tGoal state: " << goal << endl;
		myfile <<"\tInitial heuristic " << mnp.HCost(original, goal) << endl;
		// A*
		vector<AStarOpenClosedDataWithF<MNPuzzleState<4, 4>>> astarOpenList;
		if (AstarRun)
		{
			myfile <<"\t\t_A*_\n";
			TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&mnp, start, goal, astarPath, secondsLimit);
			ASTARstatesQuantityBound = astar.getMemoryStatesUse();
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tA* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % mnp.GetPathLength(astarPath) %
				   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % timer.GetElapsedTime();
				myfile << boost::format("\t\t\tI-A* ; %llu expanded;\n") % astar.getIAstarExpansions();	
			}
			else{
				myfile << boost::format("\t\t\tA* failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tI-A* failed after because A* failed\n";	
			}
		}
		// MM
		if (MMRun)
		{
			myfile << "\t\t_MM_\n";				
			MM<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> mm;
			goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = mm.GetPath(&mnp, start, goal, &mnp, &mnp2, mmPath, secondsLimit);
			MMstatesQuantityBound = mm.getMemoryStatesUse();
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tMM found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % mnp.GetPathLength(mmPath) %
					   mm.GetNodesExpanded() % mm.GetNecessaryExpansions() % MMstatesQuantityBound % timer.GetElapsedTime();
				myfile << boost::format("\t\t\tI-MM ; %llu expanded;\n") % mm.getIMMExpansions();				
			}
			else{
				myfile << boost::format("\t\t\tMM failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tI-MM failed because MM failed\n";
			}
		}
		if(MMRun && AstarRun){
			statesQuantityBound = min(ASTARstatesQuantityBound, MMstatesQuantityBound);
		}
		else if(MMRun && !AstarRun){
			statesQuantityBound = MMstatesQuantityBound;
		}
		else if(!MMRun && AstarRun){
			statesQuantityBound = ASTARstatesQuantityBound;
		}
		// IDA*
		if (IDAstarRun)
		{
			myfile << "\t\t_IDA*_\n";
			IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
			goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = idastar.GetPath(&mnp, start, goal, idaPath, secondsLimit);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tIDA* found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % mnp.GetPathLength(idaPath) %
				   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % timer.GetElapsedTime();
				myfile << boost::format("\t\t\tD-A* ; %llu expanded;\n") % idastar.getDAstarExpansions();
			}
			else{
				myfile << boost::format("\t\t\tIDA* failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-A* failed because IDA* failed\n";
			}					   
		}
		//double percentages[6] = {1, 0.9, 0.75, 0.5, 0.25, 0.1};
		double percentages[5] = {0.9, 0.75, 0.5, 0.25, 0.1};
		long stateSize = sizeof(original);
		
		if(MMpIDMM){
			myfile << "\t\t_MM+IDMM_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforMMpIDMM = statesQuantityBound*percentage;
				MM<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> mm;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = mm.GetPath(&mnp, start, goal, &mnp, &mnp2, mmPath, secondsLimit, statesQuantityBoundforMMpIDMM);
				unsigned long nodesExpanded = mm.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = mm.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tMM+IDMM MM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDMM % stateSize % percentage % mnp.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % timer.GetElapsedTime();
				}
				else{
					IDMM<MNPuzzleState<4, 4>, slideDir, false> idmm(idmmF2fFlag,isConsistent);
					MNPuzzleState<4, 4> midState;
					bool solved = idmm.GetMidStateFromLists(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), mm.getLastBound(), mm.GetForwardItems(), mm.GetBackwardItems());
					nodesExpanded += idmm.GetNodesExpanded();
					necessaryNodesExpanded += idmm.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tMM+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDMM % stateSize % percentage % idmm.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tMM+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforMMpIDMM % stateSize % percentage % timer.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if(ASTARpIDMM){
			myfile << "\t\t_A*+IDMM_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDMM = statesQuantityBound*percentage;
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDMM);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDMM A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % mnp.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % timer.GetElapsedTime();
				}
				else{
					IDMM<MNPuzzleState<4, 4>, slideDir, false> idmm(idmmF2fFlag,isConsistent);
					MNPuzzleState<4, 4> midState;
					bool solved = idmm.GetMidStateFromForwardList(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), astar.getStatesList());
					nodesExpanded += idmm.GetNodesExpanded();
					necessaryNodesExpanded += idmm.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % idmm.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % timer.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (AstarPIDAstarRun){
			myfile << "\t\t_Astar+IDAstar_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDAS = statesQuantityBound*percentage;
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDAS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDA* A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % mnp.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
					solved = idastar.ASpIDA(&mnp, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime());
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % timer.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (AstarPIDAstarReverseRun){
			myfile << "\t\t_Astar+IDAstar+Reverse_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDARS = statesQuantityBound*percentage;
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS, false);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDA*_Reverse A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % mnp.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
					solved = idastar.ASpIDArev(&mnp, goal, start, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime());
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
		// MBBDS
		if (MBBDSRun && (threePhase || twoPhase)){
			bool doThree = threePhase;
			bool doTwo = twoPhase;
			for(int i=0;i<2;i++){
				if((i==0 && doThree) || (i==1 && twoPhase)){
					myfile << boost::format("\t\t_MBBDS(ThreePhase=%d)_\n")% int(doThree);
					bool solved;
					unsigned long nodesExpanded;
					for(double percentage : percentages){
						timer.StartTimer();
						unsigned long statesQuantityBoundforMBBDS = statesQuantityBound*percentage;
						solved = false;
						unsigned long nodesExpanded = 0;
						double lastBound = 0;
						MM<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> mm;
						if(doThree){
							goal.Reset();
							start = original;
							solved = mm.GetPath(&mnp, start, goal, &mnp, &mnp2, mmPath, secondsLimit, statesQuantityBoundforMBBDS);
							nodesExpanded += mm.GetNodesExpanded();
							lastBound = mm.getLastBound();
						}
						if(solved){
							timer.EndTimer();
							myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) MM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % mnp.GetPathLength(mmPath) %
							   nodesExpanded % mm.GetNecessaryExpansions() % 0 % timer.GetElapsedTime();
						}
						else{
							MBBDS<MNPuzzleState<4, 4>, slideDir, MbbdsBloomFilter<MNPuzzleState<4, 4>, STPHasher>, false> mbbds(statesQuantityBoundforMBBDS) ;
							goal.Reset();
							start = original;
							solved = mbbds.GetMidState(&mnp, start, goal, midState, secondsLimit - timer.GetElapsedTime(), int(lastBound));
							nodesExpanded += mbbds.GetNodesExpanded();
							lastBound = mbbds.getLastBound();
							if(solved){
								timer.EndTimer();
								myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) MBBDS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % mbbds.getPathLength() % nodesExpanded % mbbds.GetNecessaryExpansions() % mbbds.getIterationNum() % timer.GetElapsedTime();
							}
							else{
								IDMM<MNPuzzleState<4, 4>, slideDir, false> idmm(idmmF2fFlag,isConsistent);
								goal.Reset();
								start = original;
								solved = idmm.GetMidState(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), int(lastBound));
								nodesExpanded += idmm.GetNodesExpanded();
								timer.EndTimer();
								if(solved){
									myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % idmm.getPathLength() % nodesExpanded % mbbds.GetNecessaryExpansions() % mbbds.getIterationNum() % timer.GetElapsedTime();
								}
								else{
									myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs and %d iterations\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % timer.GetElapsedTime() % mbbds.getIterationNum();
									break;
								} 
							}
						}
					}
				}
				doThree = false;
			}
		}
		//IDMM
		if(IDMMRun)
		{
			myfile << "\t\t_IDMM_\n";
			IDMM<MNPuzzleState<4, 4>, slideDir, false> idmm(idmmF2fFlag,isConsistent);
			goal.Reset();
			start = original;
			MNPuzzleState<4, 4> midState;
			timer.StartTimer();
			bool solved = idmm.GetMidState(&mnp, start, goal, midState, secondsLimit);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tIDMM found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed; ") % idmm.getPathLength() %
				   idmm.GetNodesExpanded() % idmm.GetNodesTouched() % idmm.GetNecessaryExpansions() % timer.GetElapsedTime();
				myfile << "Mid state: " << midState << endl;
				myfile << boost::format("\t\t\tD-MM ; %llu expanded;\n") % idmm.getDMMExpansions();
			}
			else{
				myfile << boost::format("\t\t\tIDMM failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-MM failed because IDMM failed\n";
			}   				
		}
	}
}