#ifndef TestInfo_H
#define TestInfo_H

#include <iostream>
#include <string>
#include "CsvSerializable.h"

using namespace std;

class TestInfo
{
public:
    string m_testDescription;
    unsigned int m_testId;
    string m_startState;
    string m_goalState;

    TestInfo(string testDescription, unsigned int testId, string startState, string goalState);
    static string csvSerializeHeaders();
    string csvSerialize();
};

#endif