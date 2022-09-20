/**
 * @file RealFunc.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief Real-valued function
 * @version 0.1
 * @date 2022-09-20
 */
#include<math.h>
class Func{
public:
    virtual double operator()(double x)=0;
};
