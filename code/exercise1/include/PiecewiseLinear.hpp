/**
 * @file PiecewiseLinear.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief PiecewiseLinear function for 1D FEM.
 * @version 0.1
 * @date 2022-09-20
 */
#ifndef _PiecewisePoly
#define _PiecewisePoly
#include "MeshGrid.hpp"
#include <iostream>
#include <ostream>
#include <sstream>
#include <fstream>

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
    PiecewiseLinear() = default;
    PiecewiseLinear(const std::vector<double>& anodes,const std::vector<double> avalue):nodes(anodes),value(avalue){};
    PiecewiseLinear& operator=(const PiecewiseLinear& p)=default;
    friend PiecewiseLinear operator*(double r,const PiecewiseLinear& p){
        std::vector<double> newval;
        std::vector<double> newnds=p.nodes;
        for(auto v:p.value){
            newval.push_back(v*r);
        }
        return PiecewiseLinear(newnds,newval);
    }
    friend PiecewiseLinear operator+(const PiecewiseLinear& P1,const PiecewiseLinear& P2){
        std::vector<double> newnds;
        std::vector<double> newval;
        int i=0,j=0;
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
    double operator()(double x) const{
        int n=nodes.size();
        if(x==nodes[n-1]){
            return value[n-1];
        }
        if(x==nodes[0]){
            return value[0];
        }
        if(x<nodes[0]||x>nodes[n-1]){
            return 0;
        }
        int index=0;
        for(index;nodes[index]<=x;index++){
        }
        index--;
        //int index=binary_find(x,nodes,0,n-1);
        double k=(value[index+1]-value[index])/(nodes[index+1]-nodes[index]);
        double ans=value[index]+k*(x-nodes[index]);
        return ans;
    }
    double dir(double x) const{
        int n=nodes.size();
        if(x<=nodes[0]||x>=nodes[n-1]){
            return 0;
        }
        int index=0;
        for(index;nodes[index]<=x;index++){
        }
        index--;
        double k=(value[index+1]-value[index])/(nodes[index+1]-nodes[index]);
        return k;        
    }
    void draw() const{
        std::ofstream fout("output.m");
        fout<<"x = [ 0 ";
        for(int i=1;i<=100;i++){
            fout<<", "<<0.01*i;
        }
        fout<<"];"<<std::endl;
        fout<<"y = [";
        fout<<this->operator()(0);
        for(int i=1;i<=100;i++){
            fout<<", "<<this->operator()(0.01*i);
        }
        fout<<"];"<<std::endl;
        fout<<"plot (x,y)"<<std::endl;
        fout.close();
    }
};
#else
#endif