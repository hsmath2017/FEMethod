/**
 * @file uniformtest.cpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2022-09-26
 */
#include "BVP.hpp"
#include "FuncBasis.hpp"
#include "MeshGrid.hpp"
#include "PiecewiseLinear.hpp"
#include "RealFunc.hpp"
int main(){
    std::string filename="../src/Dirichletinput.json";
    BVP1D Problem(filename);
    auto PL=Problem.solve();
    PL.draw();
}