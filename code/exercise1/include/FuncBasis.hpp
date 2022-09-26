/**
 * @file FuncBasis.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief Store the B-spline basis for piecewise-linear vector space
 * @version 0.1
 * @date 2022-09-20
 */
#include "MeshGrid.hpp"
#include "PiecewiseLinear.hpp"
class FuncBasis{
private:
    std::vector<double> nodes;
    std::vector<PiecewiseLinear> BasisFunc;
public:
    FuncBasis() = default;
    FuncBasis(const std::vector<double>& _nodes){
        for(auto v:_nodes){
            nodes.push_back(v);
        }
        int n=_nodes.size();
        for(int i=0;i<=n-1;i++){
            std::vector<double> vals;
            vals.resize(n);
            vals[i]=1;
            PiecewiseLinear PL(nodes,vals);
            BasisFunc.push_back(PL);
        }
    }
};