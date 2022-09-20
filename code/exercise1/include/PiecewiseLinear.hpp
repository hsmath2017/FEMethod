/**
 * @file PiecewiseLinear.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief PiecewiseLinear function for 1D FEM.
 * @version 0.1
 * @date 2022-09-20
 */
#include "MeshGrid.hpp"
/**
 * @brief Find the position of node nearest to x (when smaller than x)
 * @param [in] x input pos
 * @param [in] vec position vector
 * @param [in] left left-side index
 * @param [in] right right-side index
 * @return int return index
 * 
 * @details 
 */
int binary_find(double x,const std::vector<double>& vec,int left,int right){
    if(left==right){
        return left;
    }
    if(left==right-1){
        return left;
    }
    int mid=(left+right)/2;
    if(x<vec[mid]){
        return binary_find(x,vec,left,mid-1);
    }else{
        return binary_find(x,vec,mid,right);
    }
}
/**
 * @brief Store Piecewise Linear Functions.
 * 
 * @details 
 */
class PiecewiseLinear{
private:
    std::vector<double> nodes;
    std::vector<double> value;
public:
    PiecewiseLinear(const std::vector<double>& anodes,const std::vector<double> avalue):nodes(anodes),value(avalue){};
    PiecewiseLinear operator*(double r){
        for(auto& v:value){
            v=v*r;
        }
        return *this;
    }
    friend PiecewiseLinear operator+(const PiecewiseLinear& P1,const PiecewiseLinear& P2){
        std::vector<double> newnds;
        std::vector<double> newval;
        int i,j=0;
        auto nd1=P1.nodes;
        auto nd2=P2.nodes;
        auto val1=P1.value;
        auto val2=P2.value;
        while(i<nd1.size()&&j<nd2.size()){
            if(nd1[i]<nd2[j]){
                newnds.push_back(nd1[i]);
                newval.push_back(val1[i]);
                i++;
            }else if(nd1[i]>nd2[j]){
                newnds.push_back(nd2[j]);
                newval.push_back(val2[j]);
                j++;
            }else{
                newnds.push_back(nd1[i]);
                newval.push_back(val1[i]+val2[j]);
                i++;
                j++;
            }
        }
        while(i<nd1.size()){
            newnds.push_back(nd1[i]);
            newval.push_back(val1[i]);
            i++;
        }
        while(j<nd2.size()){
            newnds.push_back(nd2[j]);
            newval.push_back(val2[j]);
            j++;
        }
        PiecewiseLinear PL(newnds,newval);
        return PL;
    }
    double operator()(double x){
        int n=nodes.size();
        if(x<=nodes[0]||x>=nodes[n-1]){
            return 0;
        }
        int index=binary_find(x,nodes,0,n-1);
        double k=(value[index+1]-value[index])/(nodes[index+1]-nodes[index]);
        double ans=value[index]+k*(x-nodes[index]);
        return ans;
    }
};