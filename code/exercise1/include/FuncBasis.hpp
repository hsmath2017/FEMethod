/**
 * @file FuncBasis.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief Store the B-spline basis for piecewise-linear vector space
 * @version 0.1
 * @date 2022-09-20
 */
#include "MeshGrid.hpp"
#include "PiecewiseLinear.hpp"
class FuncBasis:public Mesh1D{
private:
    std::vector<PiecewiseLinear> BasisFunc;
public:
    using Mesh1D::meshgrid;
    using Mesh1D::LeftSide;
    using Mesh1D::RightSide;
    using Mesh1D::size;
    FuncBasis(int nsize,double Left,double Right,bool isuniform=true,const std::vector<double>& mesh={}):Mesh1D(nsize,Left,Right,isuniform,mesh){
        for(int i=0;i<mesh.size();i++){
            std::vector<double> nodes;
            std::vector<double> val={0,1,0};
            if(i==0){
                nodes={Left,mesh[i],mesh[i+1]};
            }else if(i==mesh.size()-1){
                nodes={mesh[i-1],mesh[i],Right};
            }else{
                nodes={mesh[i-1],mesh[i],mesh[i+1]};
            }
            PiecewiseLinear P(nodes,val);
            BasisFunc.push_back(P);
        }
    }
};