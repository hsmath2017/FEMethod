/**
 * @file RealFunc.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief Real-valued function
 * @version 0.1
 * @date 2022-09-20
 */
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
    double operator()(double x){return 2*sin(x);}
};

std::map<int,Func*> mp={{1,new constFunc()},{2,new testhomogenious()}};