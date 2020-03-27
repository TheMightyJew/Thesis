//
//  BidirPancake.cpp
//  hog2 glut
//
//  Created by Nathan Sturtevant on 2/7/17.
//  Copyright © 2017 University of Denver. All rights reserved.
//

#include "BidirPancake.h"
#include "PancakePuzzle.h"
#include "TemplateAStar.h"
#include "NBS.h"
#include "IDAStar.h"
#include "MM.h"
#include "BSStar.h"
#include "PancakeInstances.h"
#include "WeightedVertexGraph.h"
#include "HeuristicError.h"
#include "MBBDS.h"
#include "MyBloomFilter.h"
#include "PancakeHasher.h"
#include <iostream>
#include <math.h> 



const int S = 10; // must be factor of sizes below

void TestPancakeTR();
void TestPancakeRandom();
void TestPancakeHard(int gap);
void TestRob();
void TestVariants();
void TestError();

void TestPancake()
{
//	TestRob();
	TestPancakeRandom();
	//TestPancakeHard(0); // GAP heuristic #
	//TestPancakeHard(1);
	//TestPancakeHard(2);
	//TestPancakeHard(3);
	//TestPancakeHard(4);
	//TestPancakeHard(pancakes_num); // Heuristic 0
	
//	TestError();
//	TestVariants();
	exit(0);
}

const int pancakes_num = 16;
int all_problems_num = 1;


void TestRob()
{
//	0 3 2 1
	PancakePuzzleState<4> start;
	PancakePuzzleState<4> goal;
	PancakePuzzle<4> cake(1);
	ZeroHeuristic<PancakePuzzleState<4>> z;
	std::vector<PancakePuzzleState<4>> path;
	start.puzzle[0] = 0;
	start.puzzle[1] = 3;
	start.puzzle[2] = 2;
	start.puzzle[3] = 1;
	goal.puzzle[0] = 1;
	goal.puzzle[1] = 3;
	goal.puzzle[2] = 2;
	goal.puzzle[3] = 0;
	NBS<PancakePuzzleState<4>, PancakePuzzleAction, PancakePuzzle<4>> nbs;
	MM<PancakePuzzleState<4>, PancakePuzzleAction, PancakePuzzle<4>> mm;
	mm.GetPath(&cake, start, goal, &cake, &cake, path);
	printf("MM: %lld expansions\n", mm.GetNodesExpanded());
	mm.GetPath(&cake, start, goal, &z, &z, path);
	printf("MM0: %lld expansions\n", mm.GetNodesExpanded());
	
	exit(0);
}

void TestPancakeTR()
{
	// multiples of 5
	int arrangement[] = {0,2,4,1,3,5,7,9,6,8,10,12,14,11,13,15,17,19,16,18,20,22,24,21,23,25,27,29,26,28,30,32,34,31,33,35,37,39,36,38,40,42,44,41,43,45,47,49,46,48,50,52,54,51,53,55,57,59,56,58,60,62,64,61,63,65,67,69,66,68,70,72,74,71,73,75,77,79,76,78,80,82,84,81,83,85,87,89,86,88,90,92,94,91,93,95,97,99,96,98,};
	// multiples of 9
//	const int arrangement[] = {0,4,7,2,5,8,3,6,1,9,13,16,11,14,17,12,15,10,18,22,25,20,23,26,21,24,19,27,31,34,29,32,35,30,33,28,36,40,43,38,41,44,39,42,37,45,49,52,47,50,53,48,51,46,54,58,61,56,59,62,57,60,55,63,67,70,65,68,71,66,69,64,72,76,79,74,77,80,75,78,73,81,85,88,83,86,89,84,87,82,90,94,97,92,95,98,93,96,91};

	for (int gap = 0; gap < 10; gap++)
	{
		
		PancakePuzzleState<S> start;
		PancakePuzzleState<S> goal;
		PancakePuzzle<S> pancake(gap);
		PancakePuzzle<S> pancake2(gap);
		
		NBS<PancakePuzzleState<S>, PancakePuzzleAction, PancakePuzzle<S>> nbs;
		MM<PancakePuzzleState<S>, PancakePuzzleAction, PancakePuzzle<S>> mm;
		TemplateAStar<PancakePuzzleState<S>, PancakePuzzleAction, PancakePuzzle<S>> astar;
		IDAStar<PancakePuzzleState<S>, PancakePuzzleAction, false> idastar;
		
		std::vector<PancakePuzzleState<S>> nbsPath;
		std::vector<PancakePuzzleState<S>> astarPath;
		std::vector<PancakePuzzleState<S>> mmPath;
		std::vector<PancakePuzzleAction> idaPath;
		Timer t1, t2, t3, t4;
		
		
		goal.Reset();
		for (int x = 0; x < S; x++)
			start.puzzle[x] = arrangement[x];
		t1.StartTimer();
		astar.GetPath(&pancake, start, goal, astarPath);
		t1.EndTimer();
		uint64_t necessary = 0;
		double solutionCost = pancake.GetPathLength(astarPath);
		for (unsigned int x = 0; x < astar.GetNumItems(); x++)
		{
			const auto &item = astar.GetItem(x);
			if ((item.where == kClosedList) && (item.g+item.h < solutionCost))
				necessary++;
		}
		printf("A* found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", pancake.GetPathLength(astarPath),
			   astar.GetNodesExpanded(), necessary, t1.GetElapsedTime());
		
		goal.Reset();
		for (int x = 0; x < S; x++)
			start.puzzle[x] = arrangement[x];
		t2.StartTimer();
		nbs.GetPath(&pancake, start, goal, &pancake, &pancake2, nbsPath);
		t2.EndTimer();
		printf("NBS found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", pancake.GetPathLength(nbsPath),
			   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), t2.GetElapsedTime());
		
		goal.Reset();
		for (int x = 0; x < S; x++)
			start.puzzle[x] = arrangement[x];
		t3.StartTimer();
		idastar.GetPath(&pancake, start, goal, idaPath);
		t3.EndTimer();
		printf("IDA* found path length %ld; %llu expanded; %llu generated; %1.4fs elapsed\n", idaPath.size(),
			   idastar.GetNodesExpanded(), idastar.GetNodesTouched(), t3.GetElapsedTime());
		
		
		goal.Reset();
		for (int x = 0; x < S; x++)
			start.puzzle[x] = arrangement[x];
		t4.StartTimer();
		mm.GetPath(&pancake, start, goal, &pancake, &pancake2, mmPath);
		t4.EndTimer();
		printf("MM found path length %1.0f; %llu expanded; %1.4fs elapsed\n", pancake.GetPathLength(mmPath),
			   mm.GetNodesExpanded(), t4.GetElapsedTime());
		
		printf("Problem & IDA* & & A* & & MM & & NBS* & \\\\\n");
		printf("%d G-%d & %llu & %1.4fs & %llu & %1.4fs & %llu & %1.4fs & %llu & %1.4fs \\\\ \n", S, gap,
			   idastar.GetNodesExpanded(), t3.GetElapsedTime(),
			   astar.GetNodesExpanded(), t1.GetElapsedTime(),
			   mm.GetNodesExpanded(), t4.GetElapsedTime(),
			   nbs.GetNodesExpanded(), t2.GetElapsedTime());
	}
	exit(0);
}

const int N = pancakes_num;
void TestPancakeRandom()
{
	int singleGap = 0;
	for (int gap = singleGap; gap <= singleGap; gap = gap + 1)
	{
		printf("\nTestPancakeRandom:(Pancakes: %d, Gap: %d)\n", N, gap);
		srandom(2017218);
		PancakePuzzleState<N> start;
		PancakePuzzleState<N> original;
		PancakePuzzleState<N> goal;
		PancakePuzzle<N> pancake(gap);
		PancakePuzzle<N> pancake2(gap);
		pancake.SetUseRealValueEdges(false);
		
		std::vector<PancakePuzzleState<N>> nbsPath;
		std::vector<PancakePuzzleState<N>> bsPath;
		std::vector<PancakePuzzleState<N>> astarPath;
		std::vector<PancakePuzzleState<N>> mmPath;
		std::vector<PancakePuzzleState<N>> idaPath;
		Timer t1, t2, t3, t4, t5, t6,t7;
		int problems_num = 1;
		for (int count = 0; count < problems_num; count++)
		{
			srandom(random());
			
			goal.Reset();
			original.Reset();
			for (int x = 0; x < N; x++)
				std::swap(original.puzzle[x], original.puzzle[x+random()%(N-x)]);
			
			printf("\tProblem %d of %d\n", count+1, problems_num);
			std::cout << "\tStart state: " << original << std::endl;
			std::cout << "\tGoal state: " << goal << std::endl;
			std::cout <<"\tInitial heuristic " << pancake.HCost(original, goal) << std::endl;
			// A*
			if (1)
			{
				TemplateAStar<PancakePuzzleState<N>, PancakePuzzleAction, PancakePuzzle<N>> astar;
				start = original;
				t1.StartTimer();
				astar.GetPath(&pancake, start, goal, astarPath);
				t1.EndTimer();
				printf("\t\tGAP-%d A* found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", gap, pancake.GetPathLength(astarPath),
					   astar.GetNodesExpanded(), astar.GetNecessaryExpansions(), t1.GetElapsedTime());
			}
			
			// NBS
			if (0)
			{
				NBS<PancakePuzzleState<N>, PancakePuzzleAction, PancakePuzzle<N>> nbs;
				goal.Reset();
				start = original;
				t2.StartTimer();
				nbs.GetPath(&pancake, start, goal, &pancake, &pancake2, nbsPath);
				t2.EndTimer();
				printf("GAP-%d NBS found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed; %f meeting\n", gap, pancake.GetPathLength(nbsPath),
					   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), t2.GetElapsedTime(), nbs.GetMeetingPoint());
			}
			
			// BS*
			if (0)
			{
				BSStar<PancakePuzzleState<N>, PancakePuzzleAction, PancakePuzzle<N>> bs;
				goal.Reset();
				start = original;
				t2.StartTimer();
				bs.GetPath(&pancake, start, goal, &pancake, &pancake2, bsPath);
				t2.EndTimer();
				printf("GAP-%d BS* found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", gap, pancake.GetPathLength(bsPath),
					   bs.GetNodesExpanded(), bs.GetNecessaryExpansions(), t2.GetElapsedTime());
			}
			
			// IDA*
			if (1)
			{
				IDAStar<PancakePuzzleState<N>, PancakePuzzleAction, true> idastar;
				goal.Reset();
				start = original;
				t3.StartTimer();
				idastar.GetPath(&pancake, start, goal, idaPath);
				t3.EndTimer();
				printf("\t\tGAP-%d IDA* found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed\n", gap, pancake.GetPathLength(idaPath),
					   idastar.GetNodesExpanded(), idastar.GetNodesTouched(), t3.GetElapsedTime());
			}
			
			// MBBDS
			const long byte = 8;
			const long kb = pow(2,10) * byte;
			const long mb = pow(2,10) * kb;
			const unsigned long statesQuantityBound = 100;
			const long stateSize = sizeof(original);
			const long memoryBound = statesQuantityBound*stateSize/2;
			if (1)
			{
				//k=1
				MBBDS<PancakePuzzleState<N>, PancakePuzzleAction, MyBloomFilter<PancakePuzzleState<N>, PancakeHasher<N>, memoryBound, 1>, false> mbbds(statesQuantityBound, memoryBound) ;
				goal.Reset();
				start = original;
				PancakePuzzleState<N> midState;
				t6.StartTimer();
				double path_length = mbbds.GetMidState(&pancake, start, goal, midState);
				t6.EndTimer();
				printf("\t\tHARD-%d GAP-%d MBBDS(k=1) using memory for %1.0llu states(state size: %d bits) found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed; \n", count, gap, statesQuantityBound, stateSize, path_length,
					   mbbds.GetNodesExpanded(), mbbds.GetNodesTouched(), t6.GetElapsedTime());
				std::cout << "\t\t\t\Mid state: " << midState << std::endl;
			}
			if(1)
			{
				//k=3				
				MBBDS<PancakePuzzleState<N>, PancakePuzzleAction, MyBloomFilter<PancakePuzzleState<N>, PancakeHasher<N>, memoryBound, 3>, false> mbbds(statesQuantityBound, memoryBound) ;
				goal.Reset();
				start = original;
				PancakePuzzleState<N> midState;
				t7.StartTimer();
				double path_length = mbbds.GetMidState(&pancake, start, goal, midState);
				t7.EndTimer();
				printf("\t\tHARD-%d GAP-%d MBBDS(k=3) using memory for %1.0llu states(state size: %d bits) found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed; \n", count, gap, statesQuantityBound, stateSize, path_length,
					   mbbds.GetNodesExpanded(), mbbds.GetNodesTouched(), t7.GetElapsedTime());
				std::cout << "\t\t\t\Mid state: " << midState << std::endl;			
			}
			
			// MM
			if (1)
			{
				MM<PancakePuzzleState<N>, PancakePuzzleAction, PancakePuzzle<N>> mm;
				goal.Reset();
				start = original;
				t4.StartTimer();
				mm.GetPath(&pancake, start, goal, &pancake, &pancake2, mmPath);
				t4.EndTimer();
				printf("\t\tGAP-%d MM found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", gap, pancake.GetPathLength(mmPath),
					   mm.GetNodesExpanded(), mm.GetNecessaryExpansions(), t4.GetElapsedTime());
			}

			// MM0
			if (0 && gap == 3)
			{
				MM<PancakePuzzleState<N>, PancakePuzzleAction, PancakePuzzle<N>> mm;
				ZeroHeuristic<PancakePuzzleState<N>> z;
				goal.Reset();
				start = original;
				t4.StartTimer();
				mm.GetPath(&pancake, start, goal, &z, &z, mmPath);
				t4.EndTimer();
				printf("GAP-%d MM0 found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", gap, pancake.GetPathLength(mmPath),
					   mm.GetNodesExpanded(), mm.GetNecessaryExpansions(), t1.GetElapsedTime());
			}
		
		}
	}
}
const int CNT = pancakes_num;
void TestPancakeHard(int gap)
{
	srandom(2017218);
	PancakePuzzleState<CNT> start;
	PancakePuzzleState<CNT> original;
	PancakePuzzleState<CNT> goal;
	PancakePuzzle<CNT> pancake(gap);
	PancakePuzzle<CNT> pancake2(gap);
	pancake.SetUseRealValueEdges(false);
	
	std::vector<PancakePuzzleState<CNT>> nbsPath;
	std::vector<PancakePuzzleState<CNT>> bsPath;
	std::vector<PancakePuzzleState<CNT>> astarPath;
	std::vector<PancakePuzzleState<CNT>> mmPath;
	std::vector<PancakePuzzleState<CNT>> idaPath;
	Timer t1, t2, t3, t4, t5, t6;
	
	int problems_num = all_problems_num;
	printf("\nTestPancakeHard:(Pancakes: %d, Gap: %d)\n", CNT, gap);
	for (int count = 0; count < problems_num; count++)
	{
		goal.Reset();
		original.Reset();
		GetPancakeInstance(original, count);
		printf("\tProblem %d of %d\n", count+1, problems_num);
		std::cout << "\tStart state: " << original << std::endl;
		std::cout << "\tGoal state: " << goal << std::endl;
		std::cout <<"\tInitial heuristic " << pancake.HCost(original, goal) << std::endl;
		// A*
		if (1)
		{
			TemplateAStar<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> astar;
			start = original;
			t1.StartTimer();
			astar.GetPath(&pancake, start, goal, astarPath);
			t1.EndTimer();
			printf("\t\tHARD-%d-G%d A* found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", count, gap, pancake.GetPathLength(astarPath),
				   astar.GetNodesExpanded(), astar.GetNecessaryExpansions(), t1.GetElapsedTime());
//			std::unordered_map<int, bool> m;
//			for (int x = 0; x < astar.GetNumItems(); x++)
//			{
//				auto &i = astar.GetItem(x);
//				if (i.where != kClosedList)
//					continue;
//				int tmp = (((int)i.g)<<10)|(int)i.h;
//				if (m.find(tmp) == m.end())
//				{
//					m[tmp] = true;
//					printf("(%d, %d)\n", (int)i.g, (int)i.h);
//				}
//			}
		}

		// Reverse A*
		if (0)
		{
			TemplateAStar<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> astar;
			start = original;
			t1.StartTimer();
			astar.GetPath(&pancake, goal, start, astarPath);
			t1.EndTimer();
			printf("HARD-%d-G%d ReverseA* found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", count, gap, pancake.GetPathLength(astarPath),
				   astar.GetNodesExpanded(), astar.GetNecessaryExpansions(), t1.GetElapsedTime());
//			std::unordered_map<int, bool> m;
//			for (int x = 0; x < astar.GetNumItems(); x++)
//			{
//				auto &i = astar.GetItem(x);
//				if (i.where != kClosedList)
//					continue;
//				int tmp = (((int)i.g)<<10)|(int)i.h;
//				if (m.find(tmp) == m.end())
//				{
//					m[tmp] = true;
//					printf("(%d, %d)\n", (int)i.g, (int)i.h);
//				}
//			}
		}
		
		// Find minimum
		if (0)
		{
			start = original;

			std::string t = "/Users/nathanst/bidir/pancake/pancake_";
			t += std::to_string(count);
			t += "_GAP";
			t += std::to_string(gap);

			BidirectionalProblemAnalyzer<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> p(start, goal, &pancake, &pancake, &pancake);
			p.drawFullGraph = true;
			p.drawProblemInstance = false;
			p.drawMinimumVC = true;
			p.drawAllG = false;
			p.drawStatistics = false;
//			p.SaveSVG((t+"-full.svg").c_str());
			p.drawFullGraph = false;
			p.drawProblemInstance = false;
			p.drawAllG = true;
			p.drawStatistics = false;
//			p.SaveSVG((t+"-min.svg").c_str());
		}

		
		// NBS e=1
		if (0)
		{
			NBS<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>, NBSQueue<PancakePuzzleState<CNT>, 1>> nbs;
			goal.Reset();
			start = original;
			t2.StartTimer();
			nbs.GetPath(&pancake, start, goal, &pancake, &pancake2, nbsPath);
			t2.EndTimer();
			printf("HARD-%d-G%d NBSe1 found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed; %f meeting\n", count, gap, pancake.GetPathLength(nbsPath),
				   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), t2.GetElapsedTime(), nbs.GetMeetingPoint());
		}
		// NBS e=1
		if (0)
		{
			NBS<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>, NBSQueue<PancakePuzzleState<CNT>, 1, true>> nbs;
			goal.Reset();
			start = original;
			t2.StartTimer();
			nbs.GetPath(&pancake, start, goal, &pancake, &pancake2, nbsPath);
			t2.EndTimer();
			printf("HARD-%d-G%d NBSAe1 found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed; %f meeting\n", count, gap, pancake.GetPathLength(nbsPath),
				   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), t2.GetElapsedTime(), nbs.GetMeetingPoint());
		}
		// NBS
		if (0)
		{
			NBS<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>, NBSQueue<PancakePuzzleState<CNT>, 0>> nbs;
			goal.Reset();
			start = original;
			t2.StartTimer();
			nbs.GetPath(&pancake, start, goal, &pancake, &pancake2, nbsPath);
			t2.EndTimer();
			printf("HARD-%d-G%d NBSe0 found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed; %f meeting\n", count, gap, pancake.GetPathLength(nbsPath),
				   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), t2.GetElapsedTime(), nbs.GetMeetingPoint());
		}
		// NBS0
		if (0)
		{
			ZeroHeuristic<PancakePuzzleState<CNT>> z;
			NBS<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>, NBSQueue<PancakePuzzleState<CNT>, 0>> nbs;
			goal.Reset();
			start = original;
			t2.StartTimer();
			nbs.GetPath(&pancake, start, goal, &z, &z, nbsPath);
			t2.EndTimer();
			printf("HARD-%d-G%d NBS0 found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed; %f meeting\n", count, gap, pancake.GetPathLength(nbsPath),
				   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), t2.GetElapsedTime(), nbs.GetMeetingPoint());
		}
		
		// BS*
		if (0)
		{
			BSStar<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> bs;
			goal.Reset();
			start = original;
			t2.StartTimer();
			bs.GetPath(&pancake, start, goal, &pancake, &pancake2, bsPath);
			t2.EndTimer();
			printf("HARD-%d BS* found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", count, pancake.GetPathLength(bsPath),
				   bs.GetNodesExpanded(), bs.GetNecessaryExpansions(), t2.GetElapsedTime());
		}
		
		// MM
		if (1)
		{
			MM<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> mm;
			goal.Reset();
			start = original;
			t4.StartTimer();
			mm.GetPath(&pancake, start, goal, &pancake, &pancake2, mmPath);
			t4.EndTimer();
			printf("\t\tHARD-%d MM found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", count, pancake.GetPathLength(mmPath),
				   mm.GetNodesExpanded(), mm.GetNecessaryExpansions(), t4.GetElapsedTime());
		}
		
		// MM0
		if (0)
		{
			MM<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> mm;
			ZeroHeuristic<PancakePuzzleState<CNT>> z;
			goal.Reset();
			start = original;
			t4.StartTimer();
			mm.GetPath(&pancake, start, goal, &z, &z, mmPath);
			t4.EndTimer();
			printf("HARD-%d MM0 found path length %1.0f; %llu expanded; %llu necessary; %1.4fs elapsed\n", count, pancake.GetPathLength(mmPath),
				   mm.GetNodesExpanded(), mm.GetNecessaryExpansions(), t1.GetElapsedTime());
		}
		
		// IDA*
		if (1)
		{
			IDAStar<PancakePuzzleState<CNT>, PancakePuzzleAction, true> idastar;
			goal.Reset();
			start = original;
			t3.StartTimer();
			idastar.GetPath(&pancake, start, goal, idaPath);
			t3.EndTimer();
			printf("\t\tHARD-%d GAP-%d IDA* found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed\n", count, gap, pancake.GetPathLength(idaPath),
				   idastar.GetNodesExpanded(), idastar.GetNodesTouched(), t3.GetElapsedTime());
		}
		
		// MBBDS
		if (1)
		{
			const long byte = 8;
			const long kb = pow(2,10) * byte;
			const long mb = pow(2,10) * kb;
			const unsigned long statesQuantityBound = pow(10,8);
			const long stateSize = sizeof(original);
			const long memoryBound = statesQuantityBound*stateSize/2;
			MBBDS<PancakePuzzleState<CNT>, PancakePuzzleAction, MyBloomFilter<PancakePuzzleState<CNT>, PancakeHasher<CNT>, memoryBound, 3>, false> mbbds(statesQuantityBound, memoryBound) ;
			goal.Reset();
			start = original;
			PancakePuzzleState<CNT> midState;
			t6.StartTimer();
			double path_length = mbbds.GetMidState(&pancake, start, goal, midState);
			t6.EndTimer();
			printf("\t\tHARD-%d GAP-%d MBBDS using memory for %1.0llu states(state size: %d bits) found path length %1.0f; %llu expanded; %llu generated; %1.4fs elapsed; \n", count, gap, statesQuantityBound, stateSize, path_length,
				   mbbds.GetNodesExpanded(), mbbds.GetNodesTouched(), t6.GetElapsedTime());
			std::cout << "\t\t\t\Mid state: " << midState << std::endl;				
		}
	}
}

void Solve(Heuristic<PancakePuzzleState<CNT>> *h, const char *name)
{
	PancakePuzzle<CNT> pancake;
	PancakePuzzleState<CNT> start;
	PancakePuzzleState<CNT> goal;
	GetPancakeInstance(start, 11);
	goal.Reset();

	if (1)
	{
		BidirectionalProblemAnalyzer<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> p(start, goal, &pancake, h, h);
		p.drawProblemInstance = false;
		p.drawStatistics = false;
		p.drawMinimumVC = true;
		p.drawSumOnEdge = true;
		p.drawShortenedEdges = false;
		{
			p.drawProblemInstance = true;
			std::string s(name);
			s += "-instance.svg";
			p.SaveSVG(s.c_str());
			p.drawProblemInstance = false;
		}

		p.drawFullGraph = true;
		p.drawAllG = false;
		p.flipBackwardsGCost = false;
		{
			std::string s(name);
			s += "-ey-gn-fn.svg";
			p.SaveSVG(s.c_str());
		}

		p.drawFullGraph = true;
		p.drawAllG = false;
		p.flipBackwardsGCost = true;
		p.drawSumOnEdge = false;
		{
			std::string s(name);
			s += "-ey-gn-fy.svg";
			p.SaveSVG(s.c_str());
		}

		p.drawAllG = true;
		{
			std::string s(name);
			s += "-ey-gy-fy.svg";
			p.SaveSVG(s.c_str());
		}

		p.drawFullGraph = false;
		p.drawAllG = true;
		p.flipBackwardsGCost = true;
		p.drawSumOnEdge = true;
		{
			std::string s(name);
			s += "-en-gy-fy.svg";
			p.SaveSVG(s.c_str());
		}

		p.drawFullGraph = false;
		p.drawAllG = true;
		p.flipBackwardsGCost = true;
		p.drawShortenedEdges = true;
		p.drawSumOnEdge = true;
		{
			std::string s(name);
			s += "-es-gy-fy.svg";
			p.SaveSVG(s.c_str());
		}

}
	
	if (0)
	{
		std::vector<PancakePuzzleState<CNT>> nbsPath;
		NBS<PancakePuzzleState<CNT>, PancakePuzzleAction, PancakePuzzle<CNT>> nbs;
		nbs.GetPath(&pancake, start, goal, h, h, nbsPath);
		printf("NBS found path length %1.0f; %llu expanded; %llu necessary; %f meeting\n", pancake.GetPathLength(nbsPath),
			   nbs.GetNodesExpanded(), nbs.GetNecessaryExpansions(), nbs.GetMeetingPoint());
	}
}

void TestError()
{
	srandom(2017218);
	PancakePuzzleState<CNT> goal;
	PancakePuzzle<4> pancake(2);
	PancakePuzzleState<4> smallGoal;
	PancakePuzzle<CNT> pancake0(0);
	PancakePuzzle<CNT> pancake1(1);
	PancakePuzzle<CNT> pancake2(2);
	OffsetHeuristic<PancakePuzzleState<CNT>> o1(&pancake0, 1);
	OffsetHeuristic<PancakePuzzleState<CNT>> o2(&pancake0, 2);
	OffsetHeuristic<PancakePuzzleState<CNT>> o3(&pancake0, 3);
	WeightedHeuristic<PancakePuzzleState<CNT>> w9(&pancake0, 0.9);
	WeightedHeuristic<PancakePuzzleState<CNT>> w8(&pancake0, 0.8);
	WeightedHeuristic<PancakePuzzleState<CNT>> w7(&pancake0, 0.7);

	goal.Reset();
	float f;

//	f = MeasureHeuristicErrors(&pancake, smallGoal, &pancake, 3, 2, [](float i){return i < 1;});
//	printf("GAP\\2 Error percentage: %1.1f (4-pancake)\n", f*100);

	f = MeasureHeuristicErrors(&pancake0, goal, &pancake0, 3, 2, [](float i){return i <1;});
	printf("GAP\\0 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake1, goal, &pancake1, 3, 2, [](float i){return i <1;});
	printf("GAP\\1 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake2, goal, &pancake2, 3, 2, [](float i){return i <1;});
	printf("GAP\\2 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &o1, 3, 2, [](float i){return i <1;});
	printf("GAP-1 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &o2, 3, 2, [](float i){return i <1;});
	printf("GAP-2 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &o3, 3, 2, [](float i){return i <1;});
	printf("GAP-3 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w9, 3, 2, [](float i){return i <1;});
	printf("GAPx.9 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w8, 3, 2, [](float i){return i <1;});
	printf("GAPx.8 Error percentage (3,2): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w7, 3, 2, [](float i){return i <1;});
	printf("GAPx.7 Error percentage (3,2): %1.1f\n", f*100);

	printf("\n----\n\n");
	
	f = MeasureHeuristicErrors(&pancake0, goal, &pancake0, 3, 1, [](float i){return i <2;});
	printf("GAP\\0 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake1, goal, &pancake1, 3, 1, [](float i){return i <2;});
	printf("GAP\\1 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake2, goal, &pancake2, 3, 1, [](float i){return i <2;});
	printf("GAP\\2 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &o1, 3, 1, [](float i){return i <2;});
	printf("GAP-1 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &o2, 3, 1, [](float i){return i <2;});
	printf("GAP-2 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &o3, 3, 1, [](float i){return i <2;});
	printf("GAP-3 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w9, 3, 2, [](float i){return i <2;});
	printf("GAPx.9 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w8, 3, 2, [](float i){return i <2;});
	printf("GAPx.8 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w7, 3, 2, [](float i){return i <2;});
	printf("GAPx.7 Error percentage (3,1): %1.1f\n", f*100);

	
	f = MeasureHeuristicErrors(&pancake0, goal, &w9, 3, 2, [](float i){return i>=1 && i <2;});
	printf("GAPx.9 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w8, 3, 2, [](float i){return i>=1 && i <2;});
	printf("GAPx.8 Error percentage (3,1): %1.1f\n", f*100);
	f = MeasureHeuristicErrors(&pancake0, goal, &w7, 3, 2, [](float i){return i>=1 && i <2;});
	printf("GAPx.7 Error percentage (3,1): %1.1f\n", f*100);

	exit(0);
}


void TestVariants()
{
	srandom(2017218);
	PancakePuzzleState<CNT> start;
	PancakePuzzleState<CNT> original;
	PancakePuzzleState<CNT> goal;
	PancakePuzzle<CNT> pancake;
	PancakePuzzle<CNT> pancake0(0);
	PancakePuzzle<CNT> pancake1(1);
	PancakePuzzle<CNT> pancake2(2);
	OffsetHeuristic<PancakePuzzleState<CNT>> o1(&pancake0, 1);
	OffsetHeuristic<PancakePuzzleState<CNT>> o2(&pancake0, 2);
	OffsetHeuristic<PancakePuzzleState<CNT>> o3(&pancake0, 3);
	WeightedHeuristic<PancakePuzzleState<CNT>> w9(&pancake0, 0.9);
	WeightedHeuristic<PancakePuzzleState<CNT>> w8(&pancake0, 0.8);
	WeightedHeuristic<PancakePuzzleState<CNT>> w7(&pancake0, 0.7);

//	Solve(&pancake0, "/Users/nathanst/bidir/pancake/p11_G0");
//	Solve(&pancake1, "/Users/nathanst/bidir/pancake/p11_G1");
//	Solve(&pancake2, "/Users/nathanst/bidir/pancake/p11_G2");
	Solve(&pancake2, "/Users/nathanst/bidir/pancake/p11_G2-E");
//	Solve(&o1, "/Users/nathanst/bidir/pancake/p11_O1");
//	Solve(&o2, "/Users/nathanst/bidir/pancake/p11_O2");
//	Solve(&o3, "/Users/nathanst/bidir/pancake/p11_O3");
//	Solve(&w9, "/Users/nathanst/bidir/pancake/p11_W9");
//	Solve(&w8, "/Users/nathanst/bidir/pancake/p11_W8");
//	Solve(&w7, "/Users/nathanst/bidir/pancake/p11_W7");
}
