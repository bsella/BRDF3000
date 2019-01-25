#ifndef BASETEST_H
#define BASETEST_H

#include <iostream>
#include <functional>
#include <vector>
#include <string>

class BaseTest{
public:
    BaseTest();
    virtual ~BaseTest();
    bool doAllTests(std::ostream&);
protected:
    std::vector<bool(*)(std::string&)> tests;
};

#endif // BASETEST_H