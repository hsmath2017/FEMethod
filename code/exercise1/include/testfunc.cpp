#include "BVP.hpp"
#include "FuncBasis.hpp"
#include "MeshGrid.hpp"
#include "PiecewiseLinear.hpp"
#include "RealFunc.hpp"
using namespace std;
int main(){
    PiecewiseLinear PL{{0,0.5,1},{0,1,0}};
    for(int i=0;i<=10;i++){
        cout<<PL(i*0.1)<<endl;
    }
}