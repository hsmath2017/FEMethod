/**
 * @file BVP.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief 1D BVP for ODE
 * @version 0.1
 * @date 2022-09-20
 */
#include "FuncBasis.hpp"
#include "RealFunc.hpp"
#include "eigen3/Eigen/Core"
using namespace Eigen;
class BVP1D{
private:
    FuncBasis base;
    Func* RightHand;
public:
    BVP1D(Func* rhs,int nsize,double Left,double Right,bool isuniform=true,const std::vector<double>& mesh={}):base(nsize,Left,Right,isuniform,mesh){
        RightHand=rhs;
    }
};