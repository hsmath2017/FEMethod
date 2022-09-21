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
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/Dense"
using namespace Eigen;
class BVP1D{
public:
    friend class FuncBasis;
private:
    FuncBasis base;
    Func* RightHand;
public:
    BVP1D(Func* rhs,int nsize,double Left,double Right,bool isuniform=true,const std::vector<double>& mesh={}):base(nsize,Left,Right,isuniform,mesh){
        RightHand=rhs;
    }
    Eigen::VectorXd generateOverload(){
        Eigen::VectorXd ans;
        int n=base.size;
        ans.resize(n-1);
        std::vector<double> grid;
        grid.push_back(base.LeftSide);
        for(auto c:base.meshgrid){
            grid.push_back(c);
        }
        grid.push_back(base.RightSide);
        for(int i=0;i<n-1;i++){
            ans(i)=0.5*(grid[i+2]-grid[i])*RightHand->operator()(grid[i+1]);
        }
        return ans;
    }
};