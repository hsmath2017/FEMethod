/**
 * @file RealFunc.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief Real-valued function
 * @version 0.1
 * @date 2022-09-20
 */
#ifndef _RealFunc
#define _RealFunc
#include<math.h>
#include<map>
class Func{
public:
    virtual double operator()(double x)=0;
};

class constFunc:public Func{
public:
    double operator()(double x){return 1;};
};

class testhomogenious:public Func{
public:
    double operator()(double x){return (1+M_PI*M_PI)*sin(M_PI*x);}
};

class simpletest:public Func{
public:
    double operator()(double x){return 2+x-x*x;}
};

class lineartest:public Func{
public:
    double operator()(double x){return x;}
};

class func5:public Func{
public:
    double operator()(double x){return x*x-1;}
};
std::map<int,Func*> mp={{1,new constFunc()},{2,new testhomogenious()},{3,new simpletest()},{4,new lineartest()}, {5, new func5()}};
#else
#endif