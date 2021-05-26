#ifndef CsvSerializable_H
#define CsvSerializable_H
 
#include <iostream>
#include <string>

using namespace std;

class CsvSerializable
{
public:
    virtual string csvSerialize() = 0;
};
 
#endif