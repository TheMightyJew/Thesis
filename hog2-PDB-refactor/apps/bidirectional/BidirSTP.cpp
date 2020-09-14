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
#include "STPHasher.h"

#include <climits>
#include <algorithm>
#include "TemplateAStar.h"
#include "IDTHSwTrans.h"
#include "AStarOpenClosed.h"
#include "IDAStar.h"
#include "MM.h"
#include "MBBDS.h"
#include "BFBDS.h"
#include "IDBiHS.h"
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
static unsigned long statesQuantityBound;
static unsigned long statesQuantityBoundDefault = 1000000;
static int secondsLimit = 60*30;

static bool AstarRun=true;
static bool RevAstarRun=true;

static bool IDAstarRun=true;

static bool IDTHSpTrans = true;

static bool AstarPIDAstarRun=true;

static bool isSFBDS=true;


static bool AstarPIDAstarReverseRun=true;
static bool AstarPIDAstarReverseMinHRun=true;

static bool BAI=true;
static bool Max_BAI=true;

static bool MMRun=true;

static bool IDBiHSRun=true;
static bool idbihsF2fFlag=true;

static bool ASTARpIDBiHS=true;

static bool MMpIDBiHS=false;

static bool MBBDSRun=true;
static bool threePhase=true;
static bool twoPhase=false;

static bool isConsistent=true;
static bool isUpdateByWorkload=true;

static bool isDuplicateDetection=true;

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
static string filename;;


void TestSTP(string file)
{
	if(file == "")
		file = "results_" + datetime();
	filename = "Test_Results/STP/" + file + ".txt";
	cout << "running..." << endl;
	myfile.open (filename);
	
	StevenTest(1, true);

	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
	//exit(0);
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
    
    statesQuantityBound = ULONG_MAX;
    
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
        statesQuantityBound =  std::min(statesQuantityBound, ASTARstatesQuantityBound);
			}
			else{
				myfile << boost::format("\t\t\tA* failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tI-A* failed after because A* failed\n";	
			}
		}
    if (RevAstarRun)
		{
			myfile <<"\t\t_A*_\n";
			TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
			start = original;
			timer.StartTimer();
			bool solved = astar.GetPathTime(&mnp, goal, start, astarPath, secondsLimit);
			ASTARstatesQuantityBound = astar.getMemoryStatesUse();
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tRev-A* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % mnp.GetPathLength(astarPath) %
				   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % timer.GetElapsedTime();
        statesQuantityBound =  std::min(statesQuantityBound, ASTARstatesQuantityBound);
			}
			else{
				myfile << boost::format("\t\t\tRev-A* failed after %1.4fs\n") % timer.GetElapsedTime();
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
        statesQuantityBound =  std::min(statesQuantityBound, MMstatesQuantityBound);        
			}
			else{
				myfile << boost::format("\t\t\tMM failed after %1.4fs\n") % timer.GetElapsedTime();
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
		double percentages[3] = {0.5, 0.1, 0.01};
		long stateSize = sizeof(original);
		
		/*if(MMpIDBiHS){
			myfile << "\t\t_MM+IDBiHS_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforMMpIDBiHS = std::max(statesQuantityBound*percentage,2.0);
				MM<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> mm;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = mm.GetPath(&mnp, start, goal, &mnp, &mnp2, mmPath, secondsLimit, statesQuantityBoundforMMpIDBiHS);
				unsigned long nodesExpanded = mm.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = mm.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tMM+IDBiHS MM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDBiHS % stateSize % percentage % mnp.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % timer.GetElapsedTime();
				}
				else{
					IDBiHS<MNPuzzle<4, 4>, MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag,isConsistent);
					MNPuzzleState<4, 4> midState;
					//bool solved = idbihs.GetMidStateFromLists(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), mm.getLastBound(), mm.GetForwardItems(), mm.GetBackwardItems());
					bool solved = false;
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
      if(IDTHSpTrans){
      {
        myfile << "\t\t_IDTHSpTrans_\n";
        for(double percentage : percentages){
          unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound*percentage,2.0);;
          goal.Reset();
          start = original;
          timer.StartTimer();
          unsigned long nodesExpanded = 0;
          unsigned long necessaryNodesExpanded = 0;

          IDTHSwTrans<MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, isConsistent, isUpdateByWorkload, 1 , true);
          MNPuzzleState<4, 4> midState;
          bool solved = idbihs.GetPath(&mnp, start, goal, /*astarPath, */secondsLimit, statesQuantityBoundforASPIDBiHS);
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
          goal.Reset();
          start = original;
          timer.StartTimer();
          unsigned long nodesExpanded = 0;
          unsigned long necessaryNodesExpanded = 0;

          IDTHSwTrans<MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, isConsistent, isUpdateByWorkload, 1 , false);
          MNPuzzleState<4, 4> midState;
          bool solved = idbihs.GetPath(&mnp, start, goal, /*astarPath, */secondsLimit, statesQuantityBoundforASPIDBiHS);
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
		if(ASTARpIDBiHS){
			myfile << "\t\t_A*+IDBiHS_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound*percentage,2.0);
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDBiHS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDBiHS A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % mnp.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % timer.GetElapsedTime();
				}
				else{
					IDBiHS<MNPuzzle<4, 4>, MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, isConsistent, isUpdateByWorkload);
					MNPuzzleState<4, 4> midState;
					bool solved = idbihs.GetMidStateFromForwardList(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), astar.getStatesList(),true);
					nodesExpanded += idbihs.GetNodesExpanded();
					necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDBiHS IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % idbihs.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDBiHS IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % timer.GetElapsedTime();
						break;
					}	
				}
			}
		}
    if(ASTARpIDBiHS){
			myfile << "\t\t_A*+IDBiHS-NDD_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDBiHS = std::max(statesQuantityBound*percentage,2.0);
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDBiHS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDBiHS-NDD A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % mnp.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % timer.GetElapsedTime();
				}
				else{
					IDBiHS<MNPuzzle<4, 4>, MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, isConsistent, isUpdateByWorkload);
					MNPuzzleState<4, 4> midState;
					bool solved = idbihs.GetMidStateFromForwardList(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), astar.getStatesList(),false);
					nodesExpanded += idbihs.GetNodesExpanded();
					necessaryNodesExpanded += idbihs.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDBiHS-NDD IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % idbihs.getPathLength() % nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDBiHS-NDD IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDBiHS % stateSize % percentage % timer.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (AstarPIDAstarRun){
			myfile << "\t\t_Astar+IDAstar_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDAS = std::max(statesQuantityBound*percentage,2.0);
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
					solved = idastar.ASpIDA(&mnp, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime(),isDuplicateDetection);
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
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS, false);
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
					solved = idastar.ASpIDArev(&mnp, start, goal, idaPath, astar.getStatesList(),astar.getPrevF(), secondsLimit-timer.GetElapsedTime());
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
		if (BAI){
			myfile << "\t\t_BAI_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tBAI A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % mnp.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
					solved = idastar.BAI(&mnp, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime(),false);
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tBAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
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
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);
				TemplateAStar<MNPuzzleState<4, 4>, slideDir, MNPuzzle<4, 4>> astar;
				goal.Reset();
				start = original;
				timer.StartTimer();
				bool solved = astar.GetPathTime(&mnp, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					timer.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tMax_BAI A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % mnp.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % timer.GetElapsedTime();
				}
				else{
					IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
					solved = idastar.BAI(&mnp, start, goal, idaPath, astar.getStatesList(), secondsLimit-timer.GetElapsedTime(),true);
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					timer.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tMax_BAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % timer.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tMax_BAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % timer.GetElapsedTime();
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
						unsigned long statesQuantityBoundforMBBDS = std::max(statesQuantityBound*percentage,2.0);
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
							MBBDS<MNPuzzleState<4, 4>, slideDir, MbbdsBloomFilter<MNPuzzleState<4, 4>, STPHasher>, false> mbbds(statesQuantityBoundforMBBDS, isUpdateByWorkload, isConsistent) ;
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
								IDBiHS<MNPuzzle<4, 4>, MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, isConsistent, isUpdateByWorkload);
								goal.Reset();
								start = original;
								solved = idbihs.GetMidState(&mnp, start, goal, midState, secondsLimit-timer.GetElapsedTime(), int(lastBound));
								nodesExpanded += idbihs.GetNodesExpanded();
								timer.EndTimer();
								if(solved){
									myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % idbihs.getPathLength() % nodesExpanded % mbbds.GetNecessaryExpansions() % mbbds.getIterationNum() % timer.GetElapsedTime();
								}
								else{
									myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) IDBiHS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs and %d iterations\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % timer.GetElapsedTime() % mbbds.getIterationNum();
									break;
								} 
							}
						}
					}
				}
				doThree = false;
			}
		}
		//IDBiHS
		if(IDBiHSRun)
		{
			myfile << "\t\t_IDBiHS-0.5_\n";
			IDBiHS<MNPuzzle<4, 4>, MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, false, true);
			goal.Reset();
			start = original;
			MNPuzzleState<4, 4> midState;
			timer.StartTimer();
			bool solved = idbihs.GetMidState(&mnp, start, goal, midState, secondsLimit);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tIDBiHS-0.5 found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed; ") % idbihs.getPathLength() %
				   idbihs.GetNodesExpanded() % idbihs.GetNodesTouched() % idbihs.GetNecessaryExpansions() % timer.GetElapsedTime();
				myfile << "Mid state: " << midState << endl;
				myfile << boost::format("\t\t\tD-MM ; %llu expanded;\n") % idbihs.getDMMExpansions();
			}
			else{
				myfile << boost::format("\t\t\tIDBiHS-0.5 failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-MM failed because IDBiHS failed\n";
			}   				
		}
    /*
    if(IDBiHSRun)
		{
			myfile << "\t\t_IDBiHS-BW_\n";
			IDBiHS<MNPuzzle<4, 4>, MNPuzzleState<4, 4>, slideDir, false> idbihs(idbihsF2fFlag, true, true);
			goal.Reset();
			start = original;
			MNPuzzleState<4, 4> midState;
			timer.StartTimer();
			bool solved = idbihs.GetMidState(&mnp, start, goal, midState, secondsLimit);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tIDBiHS-BW found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed; ") % idbihs.getPathLength() %
				   idbihs.GetNodesExpanded() % idbihs.GetNodesTouched() % idbihs.GetNecessaryExpansions() % timer.GetElapsedTime();
				myfile << "Mid state: " << midState << endl;
				myfile << boost::format("\t\t\tD-MM ; %llu expanded;\n") % idbihs.getDMMExpansions();
			}
			else{
				myfile << boost::format("\t\t\tIDBiHS-BW failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-MM failed because IDBiHS failed\n";
			}   				
		}
    */
if(isSFBDS)
      {
			myfile << "\t\t_SFBDS1_\n";
			IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
			goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = idastar.SFBDS(&mnp, start, goal, idaPath, secondsLimit,1,false,false);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tSFBDS1 found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % mnp.GetPathLength(idaPath) %
				   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % timer.GetElapsedTime();
			}
			else{
				myfile << boost::format("\t\t\tSFBDS1 failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-A* failed because IDA* failed\n";
			}					   
		}
    
    {
			myfile << "\t\t_SFBDS1M_\n";
			IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
			goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = idastar.SFBDS(&mnp, start, goal, idaPath, secondsLimit,1,false,true);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tSFBDS1M found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % mnp.GetPathLength(idaPath) %
				   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % timer.GetElapsedTime();
			}
			else{
				myfile << boost::format("\t\t\tSFBDS1M failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-A* failed because IDA* failed\n";
			}					   
		}
    
		{
			myfile << "\t\t_SFBDSA_\n";
			IDAStar<MNPuzzleState<4, 4>, slideDir, false> idastar;
			goal.Reset();
			start = original;
			timer.StartTimer();
			bool solved = idastar.SFBDS(&mnp, start, goal, idaPath, secondsLimit,0,true,false);
			timer.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tSFBDSA found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % mnp.GetPathLength(idaPath) %
				   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % timer.GetElapsedTime();
			}
			else{
				myfile << boost::format("\t\t\tSFBDSA failed after %1.4fs\n") % timer.GetElapsedTime();
				myfile << "\t\t\tD-A* failed because IDA* failed\n";
			}					   
		}
    
    
	}
}