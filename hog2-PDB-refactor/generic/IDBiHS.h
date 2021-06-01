/*
 *  IDBiHS.h
 *
 *  Created by Steven Danishevski on MARCH 20.
 *
 */

#ifndef IDBiHS_H
#define IDBiHS_H

#include <iostream>
#include <unordered_set>
#include "SearchEnvironment.h"
#include "AStarOpenClosed.h"
#include "TemplateAStar.h"
#include "MM.h"
#include <math.h>
#include <typeinfo>

template <class environment, class state, class action, bool verbose = false>
class IDBiHS
{
	typedef AStarOpenClosed<state, AStarCompareWithF<state>, AStarOpenClosedDataWithF<state>> aStarStatesList;

public:
	IDBiHS(bool F2Fheuristics = false, bool isConsistent = false, bool isUpdateByWorkload = false, double smallestEdge = 1)
	{
		this->F2Fheuristics = F2Fheuristics;
		this->isConsistent = isConsistent;
		this->smallestEdge = smallestEdge;
		this->isUpdateByWorkload = isUpdateByWorkload;
	}
	virtual ~IDBiHS() {}
	bool GetMidState(environment *env, state fromState, state toState, state &midState, int secondsLimit = 600, double startingFBound = 0);
	bool Astar_plus_IDBiHS(environment *env, state fromState, state toState, state &midState, unsigned long statesQuantityBound, int secondsLimit = 600, bool detectDuplicates = false);
	bool IDTHSpTrans(environment *env, state fromState, state toState, state &midState, int secondsLimit = 600, unsigned long availableStorage = 0);
	bool GetMidStateFromLists(environment *env, state fromState, state toState, state &midState, int secondsLimit = 600, double startingFBound = 0, aStarStatesList forwardList = aStarStatesList(), aStarStatesList backwardList = aStarStatesList(), bool detectDuplicates = true);
	bool GetMidStateFromForwardList(environment *env, state fromState, state toState, state &midState, int secondsLimit = 600, aStarStatesList forwardList = aStarStatesList(), bool detectDuplicates = true);
	double getPathLength() { return pathLength; }
	uint64_t GetNodesExpanded() { return nodesExpanded; }
	uint64_t GetNecessaryExpansions() { return necessaryExpansions; }
	uint64_t GetNodesTouched() { return nodesTouched; }
	void ResetNodeCount() { nodesExpanded = nodesTouched = 0; }
	unsigned long getDMMExpansions() { return dMMExpansions; }

private:
	bool DoIterationForward(environment *env, state parent, state currState, double g, state &midState);
	bool DoIterationBackward(environment *env, state parent, state currState, double g, state &midState, state possibleMidState, double possibleMidStateG, double otherH = 0, double otherError = 0);
	double updateBoundByFraction(double boundToSplit, double p = 0.5, bool isInteger = true);
	double updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound, uint64_t forwardLoad, uint64_t backwardLoad);
	void UpdateNextBound(double fCost);

	state originGoal, originStart;
	unsigned long nodesExpanded, nodesTouched, dMMExpansions;
	double backwardBound, forwardBound, fBound, nextBound;
	double pathLength = std::numeric_limits<double>::max();
	unsigned long dMMLastIterExpansions = 0;
	unsigned long necessaryExpansions = 0;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	unsigned long forwardExpandedInLastIter = 0;
	unsigned long backwardExpandedInLastIter = 0;
	unsigned long availableStorage = 0;

	aStarStatesList forwardList;
	aStarStatesList backwardList;
	std::vector<AStarOpenClosedDataWithF<state>> forwardOpenList;
	std::vector<AStarOpenClosedDataWithF<state>> backwardOpenList;
	bool readyOpenLists = false;
	bool F2Fheuristics, firstBounds, isConsistent, detectDuplicates, isUpdateByWorkload;
	double minForwardError = 0;
	double minBackwardError = 0;
	double smallestEdge;
	std::chrono::steady_clock::time_point startTimeTest;
};

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::Astar_plus_IDBiHS(environment *env, state start, state goal, state &midState, unsigned long statesQuantityBound, int secondsLimit, bool detectDuplicates)
{
	TemplateAStar<state, action, environment> astar;
	std::vector<state> astarPath;
	Timer timer;
	timer.StartTimer();
	bool solved = astar.GetPathTime(env, start, goal, astarPath, secondsLimit, true, statesQuantityBound);
	nodesExpanded = astar.GetNodesExpanded();
	if (solved)
	{
		timer.EndTimer();
		necessaryExpansions = astar.GetNecessaryExpansions();
		pathLength = env->GetPathLength(astarPath);
	}
	else
	{
		solved = GetMidStateFromForwardList(env, start, goal, midState, secondsLimit - timer.GetElapsedTime(), astar.getStatesList(), detectDuplicates);
		timer.EndTimer();
	}
	nodesExpanded += astar.GetNodesExpanded();
	return solved;
}

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::GetMidStateFromForwardList(environment *env,
																			 state fromState, state toState, state &midState, int secondsLimit, aStarStatesList forwardList, bool detectDuplicates)
{
	aStarStatesList newBackwardList;
	double h = env->HCost(toState, fromState);
	newBackwardList.AddOpenNode(toState, env->GetStateHash(toState), 0 + h, 0, h);

	return GetMidStateFromLists(env, fromState, toState, midState, secondsLimit, 0, forwardList, newBackwardList, detectDuplicates);
}

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::GetMidStateFromLists(environment *env,
																	   state fromState, state toState, state &midState, int secondsLimit, double startingFBound, aStarStatesList forwardList, aStarStatesList backwardList, bool detectDuplicates)
{
	//printf("%d|%d\n", forwardList.getElements().size(), backwardList.getElements().size());
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = true;
	this->forwardList = forwardList;
	this->backwardList = backwardList;
	this->detectDuplicates = detectDuplicates;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double minF = std::numeric_limits<double>::max();
	double minFforward = std::numeric_limits<double>::max();
	double minFbackward = std::numeric_limits<double>::max();
	double maxOpenG = 0;
	double minOpenG = std::numeric_limits<double>::max();
	minForwardError = DBL_MAX;
	minBackwardError = DBL_MAX;

	for (int x = 0; x < forwardList.OpenSize(); x++)
	{
		AStarOpenClosedDataWithF<state> openState = forwardList.getElements()[forwardList.GetOpenItem(x)];
		minOpenG = std::min(minOpenG, openState.g);
		maxOpenG = std::max(maxOpenG, openState.g);
		minFforward = std::min(minFforward, openState.f);
		minForwardError = std::min(minForwardError, openState.g - env->HCost(originStart, openState.data));
		forwardOpenList.push_back(openState);
	}
	for (int x = 0; x < backwardList.OpenSize(); x++)
	{
		AStarOpenClosedDataWithF<state> openState = backwardList.getElements()[backwardList.GetOpenItem(x)];
		minFbackward = std::min(minFbackward, openState.f);
		minBackwardError = std::min(minBackwardError, openState.g - env->HCost(openState.data, originGoal));
		backwardOpenList.push_back(openState);
	}
	if (!isConsistent)
	{
		minBackwardError = minForwardError = 0;
	}
	minF = std::max(minFbackward, minFforward);
	sort(forwardOpenList.begin(), forwardOpenList.end(), [](const AStarOpenClosedDataWithF<state> &lhs, const AStarOpenClosedDataWithF<state> &rhs)
		 { return lhs.h < rhs.h; });
	sort(backwardOpenList.begin(), backwardOpenList.end(), [](const AStarOpenClosedDataWithF<state> &lhs, const AStarOpenClosedDataWithF<state> &rhs)
		 { return lhs.h < rhs.h; });

	double initialHeuristic = env->HCost(fromState, toState);
	fBound = nextBound = std::max(smallestEdge, std::max(startingFBound, std::max(initialHeuristic, minF)));
	forwardBound = std::max(minOpenG, ceil(fBound / 2) - smallestEdge);
	firstBounds = true;
	while (true)
	{
		//printf("bounds: %1.1f|%1.1f\n", fBound, forwardBound);
		dMMLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
		if (elapsed_seconds.count() >= secondsLimit)
		{
			return false;
		}
		if (verbose)
		{
			printf("\t\tBounds: %1.1f and %1.1f: ", forwardBound, backwardBound);
		}
		bool solved = false;
		for (AStarOpenClosedDataWithF<state> openState : forwardOpenList)
		{
			auto currentTime = std::chrono::steady_clock::now();
			std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
			if (elapsed_seconds.count() >= secondsLimit)
			{
				return false;
			}
			solved = DoIterationForward(env, openState.data, openState.data, openState.g, midState);
			if (solved)
			{
				break;
			}
		}
		if (verbose)
		{
			printf("Nodes expanded: %d(%d)\n", nodesExpanded - nodesExpandedSoFar, nodesExpanded);
		}
		if (solved)
		{
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
			for (AStarOpenClosedDataWithF<state> forwardState : forwardList.getElements())
			{
				if (forwardState.where == kClosedList && forwardState.g + forwardState.h < pathLength)
				{
					necessaryExpansions++;
				}
			}
			for (AStarOpenClosedDataWithF<state> backwardState : backwardList.getElements())
			{
				if (backwardState.where == kClosedList && backwardState.g + backwardState.h < pathLength)
				{
					necessaryExpansions++;
				}
			}
			return true;
		}
		else
		{
			if (!isUpdateByWorkload)
			{
				forwardBound = ceil(nextBound / 2) - smallestEdge;
			}
			else
			{
				forwardBound = updateBoundByWorkload(nextBound, fBound, forwardBound, forwardExpandedInLastIter, backwardExpandedInLastIter);
			}
			fBound = nextBound;
			forwardBound = std::max(minOpenG, forwardBound);
		}
	}
	return false;
}

template <class environment, class state, class action, bool verbose>
double IDBiHS<environment, state, action, verbose>::updateBoundByFraction(double boundToSplit, double p, bool isInteger)
{
	if (isInteger)
	{
		return ceil(p * boundToSplit) - smallestEdge;
	}
	else
	{
		return p * boundToSplit - smallestEdge;
	}
}

template <class environment, class state, class action, bool verbose>
double IDBiHS<environment, state, action, verbose>::updateBoundByWorkload(double newbound, double prevBound, double oldForwardBound, uint64_t forwardLoad, uint64_t backwardLoad)
{
	//printf("forwardLoad: %llu,backwardLoad: %llu,oldForwardBound%d\n",forwardLoad,backwardLoad,oldForwardBound);
	if (forwardLoad <= backwardLoad)
	{
		return oldForwardBound + newbound - prevBound;
	}
	else
	{
		return oldForwardBound;
	}
}

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::GetMidState(environment *env,
															  state fromState, state toState, state &midState, int secondsLimit, double startingFBound)
{
	auto startTime = std::chrono::steady_clock::now();
	this->readyOpenLists = false;
	this->forwardList = forwardList;
	this->backwardList = backwardList;
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	double initialHeuristic = env->HCost(fromState, toState);
	fBound = nextBound = std::max(smallestEdge, std::max(startingFBound, initialHeuristic));
	forwardBound = ceil(fBound / 2) - smallestEdge;
	unsigned long nodesExpandedSoFar = 0;
	unsigned long previousIterationExpansions = 0;
	while (true)
	{
		dMMLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
		if (elapsed_seconds.count() >= secondsLimit)
		{
			return false;
		}
		if (verbose)
		{
			printf("\t\tBounds: %1.1f and %1.1f: ", forwardBound, backwardBound);
		}
		//printf("%1.1f|%1.1f|%1.1f\n", fBound, smallestEdge, forwardBound);
		bool solved = DoIterationForward(env, originStart, originStart, 0, midState);
		if (verbose)
		{
			printf("Nodes expanded: %d(%d)\n", nodesExpanded - nodesExpandedSoFar, nodesExpanded);
		}
		if (solved)
		{
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
			//pathLength = std::min(pathLength, fBound);
			return true;
		}
		else
		{

			if (!isUpdateByWorkload)
			{
				forwardBound = ceil(nextBound / 2) - smallestEdge;
			}
			else
			{
				forwardBound = updateBoundByWorkload(nextBound, fBound, forwardBound, forwardExpandedInLastIter, backwardExpandedInLastIter);
			}
			fBound = nextBound;
			forwardExpandedInLastIter = 0;
			backwardExpandedInLastIter = 0;
		}
		previousIterationExpansions = nodesExpanded - nodesExpandedSoFar;
		nodesExpandedSoFar = nodesExpanded;
	}
	return false;
}

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::DoIterationForward(environment *env,
																	 state parent, state currState, double g, state &midState)
{
	double h = env->HCost(currState, originGoal);
	if (currState == originGoal && g <= fBound)
	{
		pathLength = g;
		midState = currState;
		return true;
	}
	if (fgreater(g + h + minBackwardError, fBound))
	{
		UpdateNextBound(g + h + minBackwardError);
		return false;
	}
	else if (g > forwardBound)
	{
		backwardBound = fBound - g - smallestEdge;
		double error = 0;
		if (isConsistent)
		{
			error = g - env->HCost(currState, originStart);
		}
		if (readyOpenLists)
		{
			for (AStarOpenClosedDataWithF<state> openState : backwardOpenList)
			{
				//if(openState.g + openState.h > fBound || openState.g > backwardBound){
				if (openState.g + openState.h > fBound)
				{
					UpdateNextBound(openState.g + openState.h);
					continue;
				}
				if (DoIterationBackward(env, openState.data, openState.data, openState.g, midState, currState, g, h, error))
				{
					pathLength += g;
					midState = currState;
					return true;
				}
			}
			return false;
		}
		else
		{
			if (DoIterationBackward(env, originGoal, originGoal, 0, midState, currState, g, h, error))
			{
				pathLength += g;
				midState = currState;
				return true;
			}
			else
			{
				return false;
			}
		}
	}
	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	forwardExpandedInLastIter++;
	if (g + h == fBound)
	{
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		uint64_t childID;
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent || (detectDuplicates && forwardList.getElements().size() > 1 && forwardList.Lookup(env->GetStateHash(neighbors[x]), childID) != kNotFound && (forwardList.Lookup(childID).where == kClosedList || (forwardList.Lookup(childID).where == kOpenList && forwardList.Lookup(childID).g <= g + edgeCost))))
		{
			continue;
		}
		if (DoIterationForward(env, currState, neighbors[x], g + edgeCost, midState))
		{
			return true;
		}
	}
	return false;
}

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::DoIterationBackward(environment *env,
																	  state parent, state currState, double g, state &midState, state possibleMidState, double possibleMidStateG, double otherH, double otherError)
{
	if (g > backwardBound || (smallestEdge == 0 && g == backwardBound))
	{
		if (currState == possibleMidState && g + possibleMidStateG <= fBound)
		{
			pathLength = g;
			return true;
		}
		else if (g > backwardBound)
		{
			UpdateNextBound(g + possibleMidStateG + smallestEdge);
			return false;
		}
	}

	double fPossibleBound = 0;
	if (F2Fheuristics)
	{
		double h = env->HCost(currState, possibleMidState);
		fPossibleBound = std::max(fPossibleBound, g + h + possibleMidStateG);
	}
	if (isConsistent)
	{
		double error = g - env->HCost(currState, originGoal);
		fPossibleBound = std::max(possibleMidStateG + otherH + error, fPossibleBound);
	}

	double originalH = env->HCost(currState, originStart);
	fPossibleBound = std::max(g + originalH + otherError, fPossibleBound);

	if (fgreater(fPossibleBound, fBound))
	{
		UpdateNextBound(fPossibleBound);
		return false;
	}

	std::vector<state> neighbors;
	env->GetSuccessors(currState, neighbors);
	nodesTouched += neighbors.size();
	nodesExpanded++;
	backwardExpandedInLastIter++;
	if (fPossibleBound == fBound)
	{
		dMMLastIterExpansions++;
	}
	for (unsigned int x = 0; x < neighbors.size(); x++)
	{
		uint64_t childID;
		double edgeCost = env->GCost(currState, neighbors[x]);
		if (neighbors[x] == parent || (detectDuplicates && backwardList.getElements().size() > 1 &&
									   backwardList.Lookup(env->GetStateHash(neighbors[x]), childID) != kNotFound && (backwardList.Lookup(childID).where == kClosedList || (backwardList.Lookup(childID).where == kOpenList && backwardList.Lookup(childID).g <= g + edgeCost))))
		{
			continue;
		}
		if (DoIterationBackward(env, currState, neighbors[x], g + edgeCost, midState, possibleMidState, possibleMidStateG, otherH, otherError))
		{
			return true;
		}
	}
	return false;
}

template <class environment, class state, class action, bool verbose>
void IDBiHS<environment, state, action, verbose>::UpdateNextBound(double fCost)
{
	if (!fgreater(nextBound, fBound))
	{
		nextBound = std::max(fBound, fCost);
	}
	else if (fgreater(fCost, fBound) && fless(fCost, nextBound))
	{
		nextBound = fCost;
	}
}

template <class environment, class state, class action, bool verbose>
bool IDBiHS<environment, state, action, verbose>::IDTHSpTrans(environment *env,
															  state fromState, state toState, state &midState, int secondsLimit, unsigned long availableStorage)
{
	auto startTime = std::chrono::steady_clock::now();
	nodesExpanded = nodesTouched = 0;
	originStart = fromState;
	originGoal = toState;
	this->availableStorage = availableStorage;

	double initialHeuristic = env->HCost(fromState, toState);
	fBound = nextBound = std::max(smallestEdge, initialHeuristic);
	forwardBound = ceil(fBound / 2) - smallestEdge;
	while (true)
	{
		dMMLastIterExpansions = 0;
		auto currentTime = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = currentTime - startTime;
		if (elapsed_seconds.count() >= secondsLimit)
		{
			return false;
		}
		if (verbose)
		{
			printf("\t\tBounds: %1.1f and %1.1f: ", forwardBound, backwardBound);
		}
		bool solved = DoIterationForward(env, originStart, originStart, 0, midState);
		if (verbose)
		{
			printf("Nodes expanded: %d(%d)\n", nodesExpanded - nodesExpandedSoFar, nodesExpanded);
		}
		if (solved)
		{
			dMMExpansions = previousIterationExpansions + dMMLastIterExpansions;
			necessaryExpansions += nodesExpandedSoFar;
			return true;
		}
		if (!isUpdateByWorkload)
		{
			forwardBound = ceil(nextBound / 2) - smallestEdge;
		}
		else
		{
			forwardBound = updateBoundByWorkload(nextBound, fBound, forwardBound, forwardExpandedInLastIter, backwardExpandedInLastIter);
		}
		fBound = nextBound;
	}
	return false;
}

#endif