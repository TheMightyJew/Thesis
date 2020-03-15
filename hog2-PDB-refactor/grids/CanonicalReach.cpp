//
//  CanonicalReach.cpp
//  HOG2 Demos
//
//  Created by Nathan Sturtevant on 5/15/18.
//  Copyright © 2018 NS Software. All rights reserved.
//

#include "CanonicalReach.h"

CanonicalReach::CanonicalReach(CanonicalGrid::CanonicalGrid *m, bool computeNow)
{
	this->me = m;
	x = y = 0;
	done = false;

	reach.clear();
	reach.resize(me->GetMap()->GetMapWidth()*me->GetMap()->GetMapHeight());
	search.SetStopAfterGoal(false);

	if (!computeNow)
		return;
	while (!DoneComputing())
		IncrementalCompute();
}

void CanonicalReach::Draw(Graphics::Display &display)
{
	double maxVal = *std::max_element(reach.begin(), reach.end());
	for (int xx = 0; xx < me->GetMap()->GetMapWidth(); xx++)
	{
		for (int yy = 0; yy < me->GetMap()->GetMapHeight(); yy++)
		{
			if (me->GetMap()->GetTerrainType(xx, yy) == kGround)
			{
				CanonicalGrid::xyLoc tmp(xx, yy);
				me->SetColor(reach[me->GetStateHash(tmp)]/maxVal, 0, 0);
				me->Draw(display, tmp);
			}
		}
	}
}

bool CanonicalReach::DoneComputing()
{
	return done;
}

double CanonicalReach::IncrementalCompute()
{
	if (done)
		return 1.0;
	
	Map *m = me->GetMap();
	
	if (m->GetTerrainType(x, y) == kGround)
	{
		CanonicalGrid::xyLoc l(x, y);
		ComputeCanonicalReach(l);
	}
	x++;
	if (x >= m->GetMapWidth())
	{
		y++;
		x = 0;
	}
	if (y >= m->GetMapHeight())
	{
		x = y = 0;
		done = true;
//		for (auto &f : reach)
//			printf("%1.1f ", f);
		return 1.0;
	}
	
	// Not yet complete
	return double(y*m->GetMapWidth()+x)/double(m->GetMapWidth()*m->GetMapHeight());
}

void CanonicalReach::ComputeCanonicalReach(CanonicalGrid::xyLoc start)
{
	search.GetPath(me, start, start, result);

	for (int i = 0; i < search.GetNumItems(); i++)
	{
		const AStarOpenClosedDataWithF<CanonicalGrid::xyLoc> &d = search.GetItem(i);
		me->GetSuccessors(d.data, neighbors);
		bool foundParent = false;
		// check if we have a neighbor that has this state as a parent
		for (int x = 0; x < neighbors.size(); x++)
		{
			AStarOpenClosedDataWithF<CanonicalGrid::xyLoc> n;
			if (search.GetClosedItem(neighbors[x], n))
			{
				if (n.parentID == i)
				{
					foundParent = true;
					break;
				}
			}
		}
		if (foundParent)
			continue;
		
		// trace from i back to the start state updating the reach of each state
		double perfectHCost = 0;
		AStarOpenClosedDataWithF<CanonicalGrid::xyLoc> v = search.GetItem(i);
		int currParent = i;
		while (true)
		{
			int64_t nextParent = v.parentID;
			if (nextParent == kTAStarNoNode || nextParent == currParent)
				break;
			CanonicalGrid::xyLoc next = search.GetItem(nextParent).data;
			perfectHCost += me->GCost(next, v.data);
			double r = std::min(perfectHCost, search.GetItem(nextParent).g);
			reach[me->GetStateHash(next)] = std::max(r, reach[me->GetStateHash(next)]);
			v = search.GetItem(nextParent);
			currParent = nextParent;
		}
	}}

bool CanonicalReach::ShouldNotGenerate(const xyLoc &start, const xyLoc &parent, const xyLoc &current, double gCost, const xyLoc &goal) const
{
	double r = reach[me->GetStateHash({current.x, current.y})];
	if (gCost > r  && me->HCost({current.x, current.y}, {goal.x, goal.y}) > r)
		return true;
	return false;
}


