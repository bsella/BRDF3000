#ifndef BASETEST_H
#define BASETEST_H

#include <iostream>
#include <vector>
#include <string>

class BaseTest{
public:
    BaseTest();
    virtual ~BaseTest();
    bool doAllTests(std::ostream&);
protected:
    void addTest(std::istream&(*)(std::istream&),
        const std::string&,
        const std::string&);
private:
    struct testSet{
        testSet(std::istream&(*)(std::istream&),
            const std::string&,
            const std::string&);

        std::istream&(*procedure)(std::istream&);
        const std::string dataPath;
        const std::string gtPath;
        // template<typename Scalar>
        bool execute(std::string&)const;
    };
    std::vector<testSet> tests;
};

#endif // BASETEST_H