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
#include "jsoncpp/json/json.h"
#include "sstream"
#include "fstream"
using namespace Eigen;
enum BdryType{Dirichlet = 0, Neumann = 1, Robin = 2,};
/**
 * @brief class BVP1D
 * @tparam BCType: Boundary Condition Type.
 * 
 * @details store the details for Boundary-valued Problem (with 1 dimension)
 */
template<BdryType BCType>
class BVP1D{
public:
    friend class FuncBasis;
    friend class Mesh1D;
    typedef Eigen::SparseMatrix<double> SpMat;
private:
    FuncBasis basis;
    Func* k,q,f;
    double LeftBC;
    double RightBC;
    double RobinBeta1;
    double RobinBeta2;
    Mesh1D mesh;
    /**
     * @brief Numerical Integral with algebraic accurate = 2.
     * @tparam Func function type (normally lambda expression)
     * @param [in] Left left_side
     * @param [in] Right right_side
     * @param [in] f integrated function
     * @return double integral value
     * 
     * @details 
     */
    template<class Func>
    double Numerical_Integral(double Left,double Right, Func f) const{
        return (Right-Left)*(f(Left)+f(Right)+4*f((Left+Right)/2))/6;
    };
    /**
     * @brief calculate the inner product of two basis function \phi_{i} with \phi_{j}
     * @param [in] i 
     * @param [in] j 
     * @return double 
     * 
     * @details 
     */
    double inner_prod(int i,int j) const{
        if(abs(i-j)>=2){
            return 0;
        }
        auto funcbase=basis.BasisFunc;
        auto func=[i,j,funcbase,k,q](double x){
            double ans=q->operator()(x)*funcbase[i](x)*funcbase[j](x)+k->opearator()(x)*funcbase[i].dir(x)*funcbase[j].dir(x);
            return ans;
        };
        double res=0;
        if(i==j){
            if(i==0){
                res=Numerical_Integral(LeftBC,mesh.meshgrid[0],func);
            }else if(i==mesh.size){
                res=Numerical_Integral(mesh.meshgrid[i-1],RightBC,func);
            }else{
                res=Numerical_Integral(mesh.meshgrid[i-1],mesh.meshgrid[i],func)+Numerical_Integral(mesh.meshgrid[i],mesh.meshgrid[i+1],func);
            }
        }else{
            int d=min(i,j);
            res=Numerical_Integral(mesh.meshgrid[d],mesh.meshgrid[d+1],func);
        }
        return res;
    };
    /**
     * @brief Calculate the overload vector
     * @return Eigen::VectorXd overload vector b
     * 
     * @details 
     */
    Eigen::VectorXd overload() const;

    Eigen::SparseMatrix<double> coef_matrix() const{
        int size=mesh.size+1;
        Eigen::SparseMatrix<double> ans(size,size);
        switch (BCType)
        {
        case Dirichlet:
            ans.coeffRef(0,0)=1
            ans.coeffRef(size-1,size-1)=1;
            for(int i=1;i<size-1;i++){
                for(int j=max(i-1,1);j<=min(size-2,i+1);j++){
                    ans.coeffRef(i,j)=inner_prod(i,j);
                }
            }
            break;
        case Neumann:
            for(int i=0;i<size;i++){
                for(int j=max(i-1,0);j<=min(size-1,i+1);j++){
                    ans.coeffRef(i,j)=inner_prod(i,j);
                }
            }
            break;
        case Robin:
            for(int i=0;i<size;i++){
                for(int j=max(i-1,0);j<=min(size-1,i+1;j++)){
                    ans.coeffRef(i,j)=inner_prod(i,j);
                }
            }
            ans.coeffRef(0,0)-=RobinBeta1*LeftBC;
            ans.coeffRef(size-1,size-1)-=RobinBeta2*RightBC;
            break;
        }
        return ans;
    }
public:
    /**
     * @brief Construct a new BVP1D object
     * @param [in] jsonfile ：jsonfile act as the input.
     * 
     * @details 
     */
    BVP1D(std::string jsonfile){
        std::ifstream infile(jsonfile);
        assert(infile.good());
        Json::Reader reader;
        Json::Value root;
        if(!reader.parse(root)){
            std::cout<<"read error!"<<std::endl;
            abort();
        }
        bool isuniform;
        double Left,Right;
        int size;
        std::vector<double> meshgrid;
        for(auto it:root["mesh"]){
            isuniform=it["uniform"].asBool();
            Left=it["Leftside"].asDouble();
            Right=it["Rightside"].asDouble();
            size=it["size"].asInt();
            for(auto val:it["meshgrid"]){
                meshgrid.push_back(val.asDouble());
            }
        }
        LeftBC=Left;
        RightBC=Right;
        k=mp[root["Functionk"].asInt()];
        q=mp[root["Functionq"].asInt()];
        f=mp[root["Functionf"].asInt()];
        RobinBeta1=root["RobinBeta1"].asDouble();
        RobinBeta2=root["RobinBeta2"].asDouble();
        std::string bdry=root["BdryType"].asString();
        switch (bdry)
        {
        case "Dirichlet":
            assert(BCType==Dirichlet);
            break;
        
        case "Neumann":
            assert(BCType==Neumann);
            break;
        
        case "Robin":
            assert(BCType==Robin);
            break;
        default:
            break;
        }
        mesh=Mesh1D(size,Left,Right,isuniform,meshgrid);
        basis=FuncBasis(mesh.meshgrid);
    }

    PiecewiseLinear solve() const{
        Eigen::SparseMatrix<double> mat=coef_matrix();
        Eigen::VectorXd load_vec=overload();
        Eigen::SimplicialCholesky<SpMat> chol(mat);
        Eigen::VectorXd coefs=chol.solve(load_vec);
        PiecewiseLinear P;
        for(int i=0;i<mesh.size+1;i++){
            P=P+load_vec(i)*basis.BasisFunc[i];
        }
        return P;
    }
};
template<BdryType BCType>
Eigen::VectorXd BVP1D<BCType>::overload() const{
    Eigen::VectorXd ans;
    int size=mesh.size+1;
    ans.resize(size);
    std::vector<double> nowgrid=mesh.meshgrid;
    auto funcbase=basis.BasisFunc;
    switch (BCType)
    {
    case Dirichlet:
        ans(0)=LeftBC;
        ans(size-1)=RightBC;
        for(int i=1;i<size-1;i++){
            auto func=[f,funcbase](double x){
                return f->operator()(x)*funcbase[i](x);
            };
            ans(i)=ans(i)+Numerical_Integral(nowgrid[i-1],nowgrid[i],func);
            ans(i)=ans(i)+Numerical_Integral(nowgrid[i],nowgrid[i+1],func);
        }
        ans(1)=ans(1)-inner_prod(0,1);
        ans(size-2)=ans(size-2)-inner_prod(size-2,size-1);
        break;
    default:
        ans(0)=LeftBC;
        ans(size-1)=RightBC;
        for(int i=0;i<size;i++){
            auto func=[f,funcbase](double x){
                return f->operator()(x)*funcbase[i](x);
            };
            if(i>=1){
                ans(i)+=Numerical_Integral(nowgrid[i-1],nowgrid[i],func);
            }
            if(i<size-1){
                ans(i)+=Numerical_Integral(nowgrid[i],nowgrid[i+1],func);
            }
        }
        break;
    }
    return ans;
};