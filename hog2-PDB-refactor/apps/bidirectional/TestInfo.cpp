#include "TestInfo.h"

TestInfo::TestInfo(string testDescription, unsigned int testId, string startState, string goalState) : m_testDescription(testDescription), m_testId(testId), m_startState(startState), m_goalState(goalState) {}

string TestInfo::csvSerializeHeaders(){
    return "Description, ID, Start State, Goal State";
}

string TestInfo::csvSerialize(){
    return "\"" + m_testDescription + "\",\"" + to_string(m_testId) + "\",\"" + m_startState + "\",\"" + m_goalState + "\"";
}