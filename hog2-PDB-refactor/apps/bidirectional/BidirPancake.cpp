//
//  BidirPancake.cpp
//  hog2 glut
//
//  Created by Nathan Sturtevant on 2/7/17.
//  Copyright © 2017 University of Denver. All rights reserved.
//
#include <algorithm>
#include "BidirPancake.h"
#include "PancakePuzzle.h"
#include "TemplateAStar.h"
#include "AStarOpenClosed.h"
#include "IDAStar.h"
#include "MM.h"
#include "PancakeInstances.h"
#include "MBBDS.h"
#include "FullMBBDS.h"
#include "IDMM.h"
#include "MbbdsBloomFilter.h"
#include "PancakeHasher.h"
#include <fstream>
#include <iostream>
#include <boost/format.hpp>
#include <math.h> 
#include <ctime>
using namespace std;


void StevenTest(int gap=0, int problems_num=1, bool randomPancake=true, vector<int> skipVector = vector<int>());

int all_problems_num = 100;
unsigned long MMstatesQuantityBound = 1000000;
unsigned long ASTARstatesQuantityBound = 1000000;
int secondsLimit = 60*30;
bool AstarRun=true;
bool AstarPIDAstarRun=true;
bool MMRun=true;
bool IDAstarRun=true;
bool MBBDSRun=true;
bool threePhase=true;
bool twoPhase=true;
bool IDMMRun=true;

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
string filename = "test_results/results_" + datetime() + ".txt";


void TestPancake()
{
	cout << "running..." << endl;
	myfile.open (filename);
	
	StevenTest(3, 100, false, {5, 23, 32, 43, 60, 63, 73});
	
	myfile << "completed!" << endl;
	myfile.close();
	cout << "completed!" << endl;
	exit(0);
}

void StevenTest(int gap, int problems_num, bool randomPancake, vector<int> skipVector)
{
	const int pancakes_num = 16;
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
	myfile << boost::format("TestPancakeHard:(Pancakes: %d, Gap: %d, Random: %d, Hard: %d)\n") % pancakes_num % gap % randomPancake % (!randomPancake);
	for (int count = 0; count < problems_num; count++)
	{
		if(std::find(skipVector.begin(), skipVector.end(), count+1) != skipVector.end()) {
			continue;
		}
		goal.Reset();
		original.Reset();
		if(randomPancake){
			for (int x = 0; x < pancakes_num; x++)
				swap(original.puzzle[x], original.puzzle[x+random()%(pancakes_num-x)]);
		}
		else{
			GetPancakeInstance(original, count);
		}
		myfile << boost::format("\tProblem %d of %d\n") % (count+1) % problems_num;
		myfile << "\tStart state: " << original << endl;
		myfile << "\tGoal state: " << goal << endl;
		myfile <<"\tInitial heuristic " << pancake.HCost(original, goal) << endl;
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
				myfile << boost::format("\t\t\tGAP-%d A* found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % gap % pancake.GetPathLength(astarPath) %
				   astar.GetNodesExpanded() % astar.GetNecessaryExpansions() % ASTARstatesQuantityBound % t1.GetElapsedTime();
				myfile << boost::format("\t\t\tGAP-%d I-A* ; %llu expanded;\n") % gap % astar.getIAstarExpansions();	
			}
			else{
				myfile << boost::format("\t\t\tHard-GAP-%d A* failed after %1.4fs\n") % gap % t1.GetElapsedTime();
				myfile << boost::format("\t\t\tHard-GAP-%d I-A* failed after because A* failed\n") % gap;	
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
				myfile << boost::format("\t\t\tHard-GAP-%d MM found path length %1.0f; %llu expanded; %llu necessary; using %1.0llu states in memory; %1.4fs elapsed\n") % gap % pancake.GetPathLength(mmPath) %
					   mm.GetNodesExpanded() % mm.GetNecessaryExpansions() % MMstatesQuantityBound % t4.GetElapsedTime();
				myfile << boost::format("\t\t\tHard-GAP-%d I-MM ; %llu expanded;\n") % gap % mm.getIMMExpansions();				
			}
			else{
				myfile << boost::format("\t\t\tHard-GAP-%d MM failed after %1.4fs\n") % gap % t4.GetElapsedTime();
				myfile << boost::format("\t\t\tHard-GAP-%d I-MM failed because MM failed\n") % gap;
			}
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
				myfile << boost::format("\t\t\tHard-GAP-%d IDA* found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed\n") % gap % pancake.GetPathLength(idaPath) %
				   idastar.GetNodesExpanded() % idastar.GetNodesTouched() % t3.GetElapsedTime();
				myfile << boost::format("\t\t\tHard-GAP-%d D-A* ; %llu expanded;\n") % gap % idastar.getDAstarExpansions();
			}
			else{
				myfile << boost::format("\t\t\tHard-GAP-%d IDA* failed after %1.4fs\n") % gap % t3.GetElapsedTime();
				myfile << boost::format("\t\t\tHard-GAP-%d D-A* failed because IDA* failed\n") % gap;
			}					   
		}
		//FullMBBDS
		/*if (false){
			bool threePhase = false;
			myfile << "\t\t_MBBDS_\n";
			double percentages[6] = {1, 0.9, 0.75, 0.5, 0.25, 0.1};
			long stateSize = sizeof(original);
			bool solved;
			unsigned long nodesExpanded;
			for(double percentage : percentages){
				timer.StartTimer();
				unsigned long statesQuantityBoundforMBBDS = MMstatesQuantityBound*percentage;
				FullMBBDS<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>, MbbdsBloomFilter<PancakePuzzleState<pancakes_num>, PancakeHasher<pancakes_num>>, pancakes_num, false> fullMbbds(statesQuantityBoundforMBBDS) ;
				bool threePhase = true;
				solved = fullMbbds.solve(&pancake, start, goal, midState, fullMbbdsPath, secondsLimit, threePhase);
				if(solved){
					myfile << boost::format("\t\t\tHARD-%d GAP-%d MBBDS(k=1, ThreePhase=%d) IDMM using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) found path length %1.0f; %llu expanded; %1.4fs elapsed;\n") % count % gap % int(fullMbbds.isThreePhase()) % statesQuantityBoundforMBBDS % stateSize % percentage % fullMbbds.getPathLength() % fullMbbds.GetNodesExpanded() % timer.GetElapsedTime();
				}
				else{
					myfile << boost::format("\t\t\tHARD-%d GAP-%d MBBDS(k=1, ThreePhase=%d) using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) failed after %1.4fs\n") % count % gap % int(fullMbbds.isThreePhase()) % statesQuantityBoundforMBBDS % stateSize % percentage % timer.GetElapsedTime();
					break;
				} 
			}
		}*/
		double percentages[6] = {1, 0.9, 0.75, 0.5, 0.25, 0.1};
		long stateSize = sizeof(original);
		// AstarPIDAstarRun
		if (AstarPIDAstarRun){
			myfile << "\t\t_Astar+IDAstar_\n";
			for(double percentage : percentages){
				unsigned long statesQuantityBoundforASPIDAS = ASTARstatesQuantityBound*percentage;
				TemplateAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, PancakePuzzle<pancakes_num>> astar;
				goal.Reset();
				start = original;
				t1.StartTimer();
				bool solved = astar.GetPathTime(&pancake, start, goal, astarPath, secondsLimit, true, statesQuantityBoundforASPIDAS);
				unsigned long nodesExpanded = astar.GetNodesExpanded();
				if(solved){
					t1.EndTimer();
					myfile << boost::format("\t\t\tGAP-%d A*+IDA* A* using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n") % gap % statesQuantityBoundforASPIDAS % stateSize % percentage % pancake.GetPathLength(astarPath) %
					   nodesExpanded % astar.GetNecessaryExpansions() % t1.GetElapsedTime();
				}
				else{
					IDAStar<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idastar;
					solved = idastar.GetPath(&pancake, start, goal, idaPath, secondsLimit-t1.GetElapsedTime(), true, astar.getOpenList());
					nodesExpanded += idastar.GetNodesExpanded();
					t1.EndTimer();
					if(solved){
						myfile << boost::format("\t\t\tHard-GAP-%d A*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed\n") % gap % statesQuantityBoundforASPIDAS % stateSize % percentage % idastar.getSolLength() %
						   nodesExpanded % idastar.GetNodesTouched() % t1.GetElapsedTime();
					}
					else{
						myfile << boost::format("\t\t\tHard-GAP-%d A*+IDA* IDA* using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) failed after %1.4fs\n") % gap % statesQuantityBoundforASPIDAS % stateSize % percentage % t1.GetElapsedTime();
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
						unsigned long statesQuantityBoundforMBBDS = MMstatesQuantityBound*percentage;
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
							myfile << boost::format("\t\t\tHARD-%d GAP-%d MBBDS(k=1, ThreePhase=%d) MM using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % count % gap % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % pancake.GetPathLength(mmPath) %
							   nodesExpanded % mm.GetNecessaryExpansions() % 0 % timer.GetElapsedTime();
						}
						else{
							MBBDS<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, MbbdsBloomFilter<PancakePuzzleState<pancakes_num>, PancakeHasher<pancakes_num>>, false> mbbds(statesQuantityBoundforMBBDS) ;
							goal.Reset();
							start = original;
							solved = mbbds.GetMidState(&pancake, start, goal, midState, secondsLimit - timer.GetElapsedTime(), int(lastBound));
							nodesExpanded += mbbds.GetNodesExpanded();
							lastBound = mbbds.getLastBound();
							if(solved){
								timer.EndTimer();
								myfile << boost::format("\t\t\tHARD-%d GAP-%d MBBDS(k=1, ThreePhase=%d) MBBDS using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % count % gap % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % mbbds.getPathLength() % nodesExpanded % mbbds.GetNecessaryExpansions() % mbbds.getIterationNum() % timer.GetElapsedTime();
							}
							else{
								IDMM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm;
								goal.Reset();
								start = original;
								solved = idmm.GetMidState(&pancake, start, goal, midState, secondsLimit-timer.GetElapsedTime(), int(lastBound));
								nodesExpanded += idmm.GetNodesExpanded();
								timer.EndTimer();
								if(solved){
									myfile << boost::format("\t\t\tHARD-%d GAP-%d MBBDS(k=1, ThreePhase=%d) IDMM using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) found path length %1.0f; %llu expanded; %llu necessary; %d iterations; %1.4fs elapsed;\n") % count % gap % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % idmm.getPathLength() % nodesExpanded % mbbds.GetNecessaryExpansions() % mbbds.getIterationNum() % timer.GetElapsedTime();
								}
								else{
									myfile << boost::format("\t\t\tHARD-%d GAP-%d MBBDS(k=1, ThreePhase=%d) using memory for %1.0llu states(state size: %d bits, Memory_percentage_from_MM=%1.2f) failed after %1.4fs and %d iterations\n") % count % gap % int(doThree) % statesQuantityBoundforMBBDS % stateSize % percentage % timer.GetElapsedTime() % mbbds.getIterationNum();
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
			IDMM<PancakePuzzleState<pancakes_num>, PancakePuzzleAction, false> idmm;
			goal.Reset();
			start = original;
			PancakePuzzleState<pancakes_num> midState;
			t8.StartTimer();
			bool solved = idmm.GetMidState(&pancake, start, goal, midState, secondsLimit);
			t8.EndTimer();
			if(solved){
				myfile << boost::format("\t\t\tHARD-%d GAP-%d IDMM found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed; ") % count % gap % idmm.getPathLength() %
				   idmm.GetNodesExpanded() % idmm.GetNodesTouched() % t8.GetElapsedTime();
				myfile << "Mid state: " << midState << endl;	
				myfile << boost::format("\t\t\tGAP-%d D-MM ; %llu expanded;\n") % gap % idmm.getDMMExpansions();
			}
			else{
				myfile << boost::format("\t\t\tGAP-%d IDMM failed after %1.4fs\n") % gap % t8.GetElapsedTime();
				myfile << boost::format("\t\t\tGAP-%d D-MM failed because IDMM failed\n") % gap;
			}   				
		}
	}
}