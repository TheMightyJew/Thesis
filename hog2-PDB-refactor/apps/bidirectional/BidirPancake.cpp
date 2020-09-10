//
//  BidirPancake.cpp
//  hog2 glut
//
//  Created by Nathan Sturtevant on 2/7/17.
//  Copyright © 2017 University of Denver. All rights reserved.
//
#include <climits>
#include <algorithm>
#include "BidirPancake.h"
#include "PancakePuzzle.h"
#include "TemplateAStar.h"
#include "AStarOpenClosed.h"
#include "IDAStar.h"
#include "MM.h"
#include "PancakeInstances.h"
#include "MBBDS.h"
#include "BFBDS.h"
#include "IDMM.h"
#include "IDTHSwTrans.h"
#include "MbbdsBloomFilter.h"
#include "PancakeHasher.h"
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <math.h> 
#include <ctime>
using namespace std;


static void StevenTest(int gap=0, int problems_num=1, bool randomPancake=true, vector<int> skipVector = vector<int>());

static unsigned long MMstatesQuantityBound;
static unsigned long ASTARstatesQuantityBound;
static unsigned long statesQuantityBound;
static unsigned long statesQuantityBoundDefault = 1000000;
static int secondsLimit = 60*30;

static bool AstarRun=true;
static bool RevAstarRun=true;

static bool IDAstarRun=true;

static bool AstarPIDAstarRun=true;
static bool AstarPIDAstarReverseRun=true;
static bool AstarPIDAstarReverseMinHRun=true;
static bool IDTHSpTrans = true;

static bool BAI=true;
static bool Max_BAI=true;

static bool MMRun=true;

static bool IDMMRun=true;
static bool idmmF2fFlag=true;

static bool ASTARpIDMM=true;
static bool MMpIDMM=false;

static bool MBBDSRun=true;
static bool BFBDSRUN=true;
static bool fullMBBDS=true;
static bool revAlgo=true;

static bool threePhase=true;
static bool twoPhase=false;

static bool detectDuplicate=true;
static bool isConsistent=true;
static bool isUpdateByWorkload=true;


string datetime()
{
    time_t rawtime;
    struct tm * timeinfo;
    char buffer[80];

    time (&rawtime);
    timeinfo = localtime(&rawtime);

    strftime(buffer,80,"%d-%m-%Y_%H-%M-%S",timeinfo);
    return string(buffer);
}

ofstream myfile;
string filename;


void TestPancake(string file)
{
	if(file == "")
		file = "results_" + datetime();
	filename = "Test_Results/PancakeSorting/" + file + ".txt";
	cout << "running..." << endl;
	myfile.open (filename);
	

	StevenTest(0, 5, true);
	StevenTest(1, 5, true);
	StevenTest(2, 5, true);

	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
	//exit(0);
}

void StevenTest(int gap, int problems_num, bool randomPancake, vector<int> skipVector)
{
	const int pancakes_num = 8;
	srandom(2017218);
	PancakePuzzleState<pancakes_num> start;
	PancakePuzzleState<pancakes_num> original;
	PancakePuzzleState<pancakes_num> goal;
	PancakePuzzle<pancakes_num> pancake(gap);
	PancakePuzzle<pancakes_num> pancake2(gap);
	pancake.SetUseRealValueEdges(false);
	pancake2.SetUseRealValueEdges(false);
	
	vector<PancakePuzzleState<pancakes_num>> astarPath;
	vector<PancakePuzzleState<pancakes_num>> mmPath;
	vector<PancakePuzzleState<pancakes_num>> fullMbbdsPath;
	vector<PancakePuzzleState<pancakes_num>> idaPath;
	PancakePuzzleState<pancakes_num> midState;
	
	Timer t1, t2, t3, t4, t5, timer, t7, t8;
	myfile << boost::format("TestPancake:(Pancakes: %d, Gap: %d, Random: %d, Hard: %d)\n") % pancakes_num % gap % randomPancake % (!randomPancake);
	for (int count = 0; count < problems_num; count++){
		goal.Reset();
		original.Reset();
		if(randomPancake){
			for (int x = 0; x < pancakes_num; x++)
				swap(original.puzzle[x], original.puzzle[x+random()%(pancakes_num-x)]);
		}
		else{
			GetPancakeInstance(original, count);
		}
		if(std::find(skipVector.begin(), skipVector.end(), count+1) != skipVector.end()) {
			continue;
		}
    
    /*if (count+1 != 17 || gap!=0){
      continue;
    }*/
	cout << "Running: Gap=" << gap << ", ProblemID=" << count+1 << endl;
    
		myfile << boost::format("\tProblem %d of %d\n") % (count+1) % problems_num;
		myfile << "\tStart state: " << original << endl;
		myfile << "\tGoal state: " << goal << endl;
		myfile <<"\tInitial heuristic " << pancake.HCost(original, goal) << endl;
    
    statesQuantityBound = ULONG_MAX;
    
		// A*
		vector<AStarOpenClosedDataWithF<PancakePuzzleState<pancakes_num>>> astarOpenList;
		if (AstarRun)
		{
			myfile <<"\t\t_A*_\n";
			TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
			start = original;
			t1.StartTimer();
			bool solved = astar.GetPathTime(&pancake, start, goal, astarPath, secondsLimit);
			ASTARstatesQuantityBound = astar.getMemoryStatesUse();
			t1.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tA* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % pancake.GetPathLength(astarPath) %
				   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % t1.GetElapsedTime();
				myfile << boost::format("\t\t\tI-A* ; %llu expanded;\n") % astar.getIAstarExpansions();	
        statesQuantityBound =  std::min(statesQuantityBound, ASTARstatesQuantityBound);
			}
			else{
				myfile << boost::format("\t\t\tA* failed after %1.4fs\n") % t1.GetElapsedTime();
				myfile << "\t\t\tI-A* failed after because A* failed\n";	
			}
		}
    if (RevAstarRun)
		{
			myfile <<"\t\t_Rev-A*_\n";
			TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
			start = original;
			t1.StartTimer();
			bool solved = astar.GetPathTime(&pancake, goal, start, astarPath, secondsLimit);
			ASTARstatesQuantityBound = astar.getMemoryStatesUse();
			t1.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tRev-A* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % pancake.GetPathLength(astarPath) %
				   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % t1.GetElapsedTime();
           statesQuantityBound =  std::min(statesQuantityBound, ASTARstatesQuantityBound);
			}
			else{
				myfile << boost::format("\t\t\tRev-A* failed after %1.4fs\n") % t1.GetElapsedTime();

			}
		}
		// MM
		if (MMRun)
		{
			myfile << "\t\t_MM_\n";				
			MM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> mm;
			goal.Reset();
			start = original;
			t4.StartTimer();
			bool solved = mm.GetPath(&pancake, start, goal, &pancake, &pancake2, mmPath, secondsLimit);
			MMstatesQuantityBound = mm.getMemoryStatesUse();
			t4.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tMM found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % pancake.GetPathLength(mmPath) %
					   mm.GetNodesExpanded() % mm.GetNecessaryExpansions() % MMstatesQuantityBound % t4.GetElapsedTime();
				myfile << boost::format("\t\t\tI-MM ; %llu expanded;\n") % mm.getIMMExpansions();
        statesQuantityBound =  std::min(statesQuantityBound, MMstatesQuantityBound);				
			}
			else{
				myfile << boost::format("\t\t\tMM failed after %1.4fs\n") % t4.GetElapsedTime();
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
			IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
			goal.Reset();
			start = original;
			t3.StartTimer();
			bool solved = idastar.GetPath(&pancake, start, goal, idaPath, secondsLimit);
			t3.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tIDA* found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % pancake.GetPathLength(idaPath) %
				   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % idastar.GetNecessaryExpansions() % t3.GetElapsedTime();
				myfile << boost::format("\t\t\tD-A* ; %llu expanded;\n") % idastar.getDAstarExpansions();
			}
			else{
				myfile << boost::format("\t\t\tIDA* failed after %1.4fs\n") % t3.GetElapsedTime();
				myfile << "\t\t\tD-A* failed because IDA* failed\n";
			}					   
		}
    
		//double percentages[6] = {1, 0.9, 0.75, 0.5, 0.25, 0.1};
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
				BFBDS<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>, MbbdsBloomFilter<PancakePuzzleState<pancakes_num>, PancakeHasher<pancakes_num>>, false> bfbds(statesQuantityBoundforMBBDS) ;
				bool threePhase = true;
				solved = bfbds.solve(&pancake, start, goal, midState, fullMbbdsPath, secondsLimit, threePhase);
				if(solved){
					myfile << boost::format("\t\t\tBFBDS(k=1, ThreePhase=%d) using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(bfbds.isThreePhase()) % statesQuantityBoundforMBBDS % stateSize % percentage % bfbds.getPathLength() % bfbds.getNodesExpanded() % bfbds.getNecessaryExpansions() % bfbds.getIterationsNum() % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tHARD-%d GAP-%d BFBDS(k=1, ThreePhase=%d) using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % count % gap % int(bfbds.isThreePhase()) % statesQuantityBoundforMBBDS % stateSize % percentage % timer.GetElapsedTime();
					break;
				} 
			}
		}
		
		/*if(MMpIDMM){
			myfile << "\t\t_MM+IDMM_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforMMpIDMM = std::max(statesQuantityBound*percentage,2.0);;
				MM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> mm;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = mm.GetPath(&pancake, start, goal, &pancake, &pancake2, mmPath, secondsLimit, statesQuantityBoundforMMpIDMM);
				unsigned long nodesExpanded = mm.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = mm.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tMM+IDMM MM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDMM % stateSize % percentage % pancake.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % t4.GetElapsedTime();
				}
				else{
					IDMM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm(idmmF2fFlag, isConsistent, isUpdateByWorkload);
					PancakePuzzleState<pancakes_num> midState;
					bool solved = idmm.GetMidStateFromLists(&pancake, start, goal, midState, secondsLimit-t1.GetElapsedTime(), mm.getLastBound(), mm.GetForwardItems(), mm.GetBackwardItems());
					nodesExpanded += idmm.GetNodesExpanded();
					necessaryNodesExpanded += idmm.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tMM+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforMMpIDMM % stateSize % percentage % idmm.getPathLength() % nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tMM+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforMMpIDMM % stateSize % percentage % t1.GetElapsedTime();
						break;
					}	
				}
			}
		}*/
		if(ASTARpIDMM){
			myfile << "\t\t_A*+IDMM_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDMM = std::max(statesQuantityBound*percentage,2.0);;
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDMM);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDMM A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % pancake.GetPathLength(mmPath) %
					   nodesExpanded % necessaryNodesExpanded  % t1.GetElapsedTime();
				}
				else{
					IDMM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm(idmmF2fFlag, isConsistent, isUpdateByWorkload);
					PancakePuzzleState<pancakes_num> midState;
					bool solved = idmm.GetMidStateFromForwardList(&pancake, start, goal, midState, secondsLimit-t1.GetElapsedTime(), astar.getStatesList(), detectDuplicate);
					nodesExpanded += idmm.GetNodesExpanded();
					necessaryNodesExpanded += idmm.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
					  myfile << boost::format("\t\t\tA*+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % idmm.getPathLength() % nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
					  myfile << boost::format("\t\t\tA*+IDMM IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % t1.GetElapsedTime();
					  break;
					}	
				}
			}
		}
    if(IDTHSpTrans){
      {
        myfile << "\t\t_IDTHSpTrans_\n";
        for(double percentage : percentages){
          unsigned long statesQuantityBoundforASPIDMM = std::max(statesQuantityBound*percentage,2.0);;
          goal.Reset();
          start = original;
          t1.StartTimer();
          unsigned long nodesExpanded = 0;
          unsigned long necessaryNodesExpanded = 0;

          IDTHSwTrans<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm(idmmF2fFlag, isConsistent, isUpdateByWorkload, 1 , true);
          PancakePuzzleState<pancakes_num> midState;
          bool solved = idmm.GetPath(&pancake, start, goal, /*astarPath, */secondsLimit, statesQuantityBoundforASPIDMM);
          nodesExpanded += idmm.GetNodesExpanded();
          necessaryNodesExpanded += idmm.GetNecessaryExpansions();
          t1.EndTimer();
          if(solved){
            myfile << boost::format("\t\t\tIDTHSpTrans IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % idmm.getPathLength() % nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
          }
          else{
            myfile << boost::format("\t\t\tIDTHSpTrans IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % t1.GetElapsedTime();
            break;
          }	
          
        }
      }
            {
        myfile << "\t\t_IDTHSpTrans_NDD_\n";
        for(double percentage : percentages){
          unsigned long statesQuantityBoundforASPIDMM = std::max(statesQuantityBound*percentage,2.0);;
          goal.Reset();
          start = original;
          t1.StartTimer();
          unsigned long nodesExpanded = 0;
          unsigned long necessaryNodesExpanded = 0;

          IDTHSwTrans<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm(idmmF2fFlag, isConsistent, isUpdateByWorkload, 1 , false);
          PancakePuzzleState<pancakes_num> midState;
          bool solved = idmm.GetPath(&pancake, start, goal, /*astarPath, */secondsLimit, statesQuantityBoundforASPIDMM);
          nodesExpanded += idmm.GetNodesExpanded();
          necessaryNodesExpanded += idmm.GetNecessaryExpansions();
          t1.EndTimer();
          if(solved){
            myfile << boost::format("\t\t\tIDTHSpTrans_NDD IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % idmm.getPathLength() % nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
          }
          else{
            myfile << boost::format("\t\t\tIDTHSpTrans_NDD IDMM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDMM % stateSize % percentage % t1.GetElapsedTime();
            break;
          }	
          
        }
      }
		}
		if (AstarPIDAstarRun){
			myfile << "\t\t_Astar+IDAstar_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDAS = std::max(statesQuantityBound*percentage,2.0);;
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDAS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDA* A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % pancake.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
				}
				else{
					IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
					solved = idastar.ASpIDA(&pancake, start, goal, idaPath, astar.getStatesList(), secondsLimit-t1.GetElapsedTime(), detectDuplicate);
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDAS % stateSize % percentage % t1.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (AstarPIDAstarReverseRun){
			myfile << "\t\t_Astar+IDAstar+Reverse_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS, false);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDA*_Reverse A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % pancake.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
				}
				else{
					IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
					solved = idastar.ASpIDArev(&pancake, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(),secondsLimit-t1.GetElapsedTime());
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDA*_Reverse IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDA*_Reverse IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % t1.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (AstarPIDAstarReverseMinHRun){
			myfile << "\t\t_Astar+IDAstar+Reverse+MinH_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);;
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS, false);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tA*+IDA*_Reverse+MinH A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % pancake.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
				}
				else{
					IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
					solved = idastar.ASpIDArev(&pancake, start, goal, idaPath, astar.getStatesList(), astar.getPrevF(),secondsLimit-t1.GetElapsedTime(),isConsistent, true);
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tA*+IDA*_Reverse+MinH IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tA*+IDA*_Reverse+MinH IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % t1.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (BAI){
			myfile << "\t\t_BAI_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);;
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tBAI A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % pancake.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
				}
				else{
					IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
					solved = idastar.BAI(&pancake, start, goal, idaPath, astar.getStatesList(), secondsLimit-t1.GetElapsedTime(), false);
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tBAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tBAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % t1.GetElapsedTime();
						break;
					}	
				}
			}
		}
		if (Max_BAI){
			myfile << "\t\t_Max_BAI_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDARS = std::max(statesQuantityBound*percentage,2.0);;
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, goal, start, astarPath, secondsLimit, true, statesQuantityBoundforASPIDARS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				unsigned long necessaryNodesExpanded = 0;
				if(solved){
					t1.EndTimer();
					necessaryNodesExpanded = astar.GetNecessaryExpansions();
					myfile << boost::format("\t\t\tMax_BAI A* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % pancake.GetPathLength(astarPath) %
					   nodesExpanded % necessaryNodesExpanded % t1.GetElapsedTime();
				}
				else{
					IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
					solved = idastar.BAI(&pancake, start, goal, idaPath, astar.getStatesList(), secondsLimit-t1.GetElapsedTime(), true);
					nodesExpanded += idastar.GetNodesExpanded();
					necessaryNodesExpanded += idastar.GetNecessaryExpansions();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tMax_BAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % necessaryNodesExpanded % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tMax_BAI IDA* using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) failed after %1.4fs\n") % statesQuantityBoundforASPIDARS % stateSize % percentage % t1.GetElapsedTime();
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
						//printf("precentage= %1.2f\n", percentage);
						timer.StartTimer();
						unsigned long statesQuantityBoundforMBBDS = std::max(statesQuantityBound*percentage,2.0);;
						solved = false;
						unsigned long nodesExpanded = 0;
						double lastBound = 0;
						MM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> mm;
						if(doThree){
							goal.Reset();
							start = original;
							solved = mm.GetPath(&pancake, start, goal, &pancake, &pancake2, mmPath, secondsLimit, statesQuantityBoundforMBBDS);
							nodesExpanded += mm.GetNodesExpanded();
							lastBound = mm.getLastBound();
						}
						if(solved){
							timer.EndTimer();
							myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) MM using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % pancake.GetPathLength(mmPath) %
							   nodesExpanded % mm.GetNecessaryExpansions() % 0 % timer.GetElapsedTime();
						}
						else{
							/*MBBDS<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, MbbdsBloomFilter<PancakePuzzleState<pancakes_num>, PancakeHasher<pancakes_num>>, false> mbbds(statesQuantityBoundforMBBDS, isUpdateByWorkload, isConsistent, revAlgo) ;*/
							MBBDS<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, MbbdsBloomFilter<PancakePuzzleState<pancakes_num>, PancakeHasher<pancakes_num>>, false> mbbds(statesQuantityBoundforMBBDS, isUpdateByWorkload, false, false) ;
							goal.Reset();
							start = original;
							solved = mbbds.GetMidState(&pancake, start, goal, midState, secondsLimit - timer.GetElapsedTime(), int(lastBound));
							nodesExpanded += mbbds.GetNodesExpanded();
							lastBound = mbbds.getLastBound();
							if(solved){
								timer.EndTimer();
								myfile << boost::format("\t\t\tMBBDS(k=1,ThreePhase=%d) MBBDS using memory for %1.0llu states(state size: %d bits, Memory_Percentage=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % mbbds.getPathLength() % nodesExpanded % mbbds.GetNecessaryExpansions() % mbbds.getIterationNum() % timer.GetElapsedTime();
							}
							else{
								IDMM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm(idmmF2fFlag, isConsistent, isUpdateByWorkload);
								goal.Reset();
								start = original;
								solved = idmm.GetMidState(&pancake, start, goal, midState, secondsLimit-timer.GetElapsedTime(), int(lastBound));
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
		if(IDMMRun){
			myfile << "\t\t_IDMM_\n";
			IDMM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm(idmmF2fFlag, isConsistent, isUpdateByWorkload);
			goal.Reset();
			start = original;
			PancakePuzzleState<pancakes_num> midState;
			t8.StartTimer();
			bool solved = idmm.GetMidState(&pancake, start, goal, midState, secondsLimit);
			t8.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tIDMM found path length %1.0f; %llu expanded; %llu generated; %llu necessary; %1.4fs elapsed; ") % idmm.getPathLength() %
				   idmm.GetNodesExpanded() % idmm.GetNodesTouched() % idmm.GetNecessaryExpansions() % t8.GetElapsedTime();
				myfile << "Mid state: " << midState << endl;
				myfile << boost::format("\t\t\tD-MM ; %llu expanded;\n") % idmm.getDMMExpansions();
			}
			else{
				myfile << boost::format("\t\t\tIDMM failed after %1.4fs\n") % t8.GetElapsedTime();
				myfile << "\t\t\tD-MM failed because IDMM failed\n";
			}   				
		}
	}
}