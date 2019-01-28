#ifndef BASETEST_H
#define BASETEST_H

#include <iostream>
#include <vector>
#include <string>
#include <sstream>

class BaseTest{
public:
    explicit BaseTest(const std::string&);
    virtual ~BaseTest();
    bool doAllTests(std::ostream&);
protected:
    void addTest(std::istringstream(*)(std::istream&),
        const std::string&,
        const std::string&,
        const std::string&);
private:
    const std::string title;
    struct testSet{
        testSet(std::istringstream(*)(std::istream&),
            const std::string&,
            const std::string&,
            const std::string&);

        std::istringstream(*procedure)(std::istream&);
        const std::string title;
        const std::string dataPath;
        const std::string gtPath;
        // template<typename Scalar>
        bool execute(std::string&)const;
    };
    std::vector<testSet> tests;
};

#endif // BASETEST_H