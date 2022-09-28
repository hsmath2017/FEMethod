/**
 * @file BVP.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief 1D BVP for ODE
 * @version 0.1
 * @date 2022-09-20
 */
#ifndef _BVP
#define _BVP
#include <iostream>
#include "FuncBasis.hpp"
#include "RealFunc.hpp"
#include "eigen3/Eigen/Core"
#include "eigen3/Eigen/Sparse"
#include "eigen3/Eigen/Dense"
#include "jsoncpp/json/json.h"
#include "sstream"
#include "fstream"
using namespace Eigen;
enum BdryType
{
    Dirichlet = 0,
    Neumann = 1,
    Robin = 2,
};
/**
 * @brief class BVP1D
 * @tparam BCType: Boundary Condition Type.
 *
 * @details store the details for Boundary-valued Problem (with 1 dimension)
 */
class BVP1D
{
public:
    typedef Eigen::SparseMatrix<double> SpMat;

private:
    FuncBasis basis;
    Func *k;
    Func *q;
    Func *f;
    double LeftBC;
    double RightBC;
    double RobinBeta1;
    double RobinBeta2;
    BdryType BCType;
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
    template <class Func>
    double Numerical_Integral(double Left, double Right, Func f) const
    {
        double coe=(Right-Left)/2;
        double x1=(Left+Right)/2+std::sqrt(3)/3*((Right-Left)/2);
        double x2=(Left+Right)/2-std::sqrt(3)/3*((Right-Left)/2);
        return coe*(f(x1)+f(x2));
    };
    /**
     * @brief calculate the inner product of two basis function \phi_{i} with \phi_{j}
     * @param [in] i
     * @param [in] j
     * @return double
     *
     * @details
     */
    double inner_prod(int i, int j) const
    {
        if (abs(i - j) >= 2)
        {
            return 0;
        }
        auto funcbase = basis.BasisFunc;
        auto func = [i, j, funcbase, this](double x)
        {
            double ans = q->operator()(x) * funcbase[i](x) * funcbase[j](x)
             + k->operator()(x) * funcbase[i].dir(x) * funcbase[j].dir(x);
            return ans;
        };
        double res = 0;
        if (i == j)
        {
            if (i == 0)
            {
                res = Numerical_Integral(mesh.meshgrid[0], mesh.meshgrid[1], func);
            }
            else if (i == mesh.size)
            {
                res = Numerical_Integral(mesh.meshgrid[i-1], mesh.meshgrid[i], func);
            }
            else
            {
                res = Numerical_Integral(mesh.meshgrid[i - 1], mesh.meshgrid[i], func) + Numerical_Integral(mesh.meshgrid[i], mesh.meshgrid[i + 1], func);
            }
        }
        else
        {
            int d = std::min(i, j);
            res = Numerical_Integral(mesh.meshgrid[d], mesh.meshgrid[d + 1], func);
        }
        return res;
    };
    /**
     * @brief Calculate the overload vector
     * @return Eigen::VectorXd overload vector b
     *
     * @details
     */
    void overload(Eigen::VectorXd &overload) const;

    void coef_matrix(Eigen::SparseMatrix<double> &mat) const
    {
        int size = mesh.size + 1;
        mat.resize(size, size);
        switch (BCType)
        {
        case Dirichlet:
            mat.coeffRef(0, 0) = 1;
            mat.coeffRef(size - 1, size - 1) = 1;
            for (int i = 1; i < size - 1; i++)
            {
                for (int j = std::max(i - 1, 1); j <= std::min(size - 2, i + 1); j++)
                {
                    mat.coeffRef(i, j) += inner_prod(i, j);
                }
            }
            break;
        case Neumann:
            for (int i = 0; i < size; i++)
            {
                for (int j = std::max(i - 1, 0); j <= std::min(size - 1, i + 1); j++)
                {
                    mat.coeffRef(i, j) = inner_prod(i, j);
                }
            }
            inner_prod(size-1,size-1);
            break;
        case Robin:
            for (int i = 0; i < size; i++)
            {
                for (int j = std::max(i - 1, 0); j <= std::min(size - 1, i + 1); j++)
                {
                    mat.coeffRef(i, j) = inner_prod(i, j);
                }
            }
            mat.coeffRef(0, 0) += RobinBeta1 ;
            mat.coeffRef(size - 1, size - 1) += RobinBeta2 ;
            break;
        }
    }

public:
    /**
     * @brief Construct a new BVP1D object
     * @param [in] jsonfile ：jsonfile act as the input.
     *
     * @details
     */
    BVP1D(std::string jsonfile)
    {
        std::ifstream infile(jsonfile);
        assert(infile.good());
        Json::Reader reader;
        Json::Value root;
        if (!reader.parse(infile, root))
        {
            std::cout << "read error!" << std::endl;
            abort();
        }
        bool isuniform;
        double Left, Right;
        int size;
        std::vector<double> meshgrid;
        auto it = root["mesh"];
        isuniform = it["uniform"].asBool();
        Left = it["leftside"].asDouble();
        Right = it["rightside"].asDouble();
        size = it["size"].asInt();
        for (auto val : it["meshgrid"])
        {
            meshgrid.push_back(val.asDouble());
        }
        LeftBC = root["LeftBC"].asDouble();
        RightBC = root["RightBC"].asDouble();
        k = mp[root["Functionk"].asInt()];
        q = mp[root["Functionq"].asInt()];
        f = mp[root["Functionf"].asInt()];
        RobinBeta1 = root["RobinBeta1"].asDouble();
        RobinBeta2 = root["RobinBeta2"].asDouble();
        std::string bdry = root["BdryType"].asString();
        mesh = Mesh1D(size, Left, Right, isuniform, meshgrid);
        basis = FuncBasis(mesh.meshgrid);
        if (bdry == "Dirichlet")
        {
            BCType = Dirichlet;
        }
        if (bdry == "Neumann")
        {
            BCType = Neumann;
        }
        if (bdry == "Robin")
        {
            BCType = Robin;
        }
    }

    PiecewiseLinear solve() const
    {
        Eigen::SparseMatrix<double> mat;
        coef_matrix(mat);
        Eigen::VectorXd load_vec;
        overload(load_vec);
        //std::cout<<mat<<std::endl;
        //std::cout<<load_vec<<std::endl;
        Eigen::SimplicialCholesky<SpMat> chol(mat);
        Eigen::VectorXd coefs = chol.solve(load_vec);
        //std::cout<<coefs<<std::endl;
        PiecewiseLinear P;
        for (int i = 0; i < mesh.size + 1; i++)
        {
            auto np=coefs(i)*basis.BasisFunc[i];
            P = P + np;
        }
        return P;
    }
};
void BVP1D::overload(Eigen::VectorXd &ans) const
{
    int size = mesh.size + 1;
    ans.resize(size);
    std::vector<double> nowgrid = mesh.meshgrid;
    auto funcbase = basis.BasisFunc;
    switch (BCType)
    {
    case Dirichlet:
        ans(0) = LeftBC;
        ans(size - 1) = RightBC;
        for (int i = 1; i < size - 1; i++)
        {
            auto func = [i, this, funcbase](double x)
            {
                return f->operator()(x) * funcbase[i](x);
            };
            ans(i) = Numerical_Integral(nowgrid[i - 1], nowgrid[i], func);
            ans(i) = ans(i) + Numerical_Integral(nowgrid[i], nowgrid[i + 1], func);
        }
        ans(1) = ans(1) - LeftBC*inner_prod(0, 1);
        ans(size - 2) = ans(size - 2) - RightBC*inner_prod(size - 2, size - 1);
        break;
    default:
        ans(0) = LeftBC;
        ans(size - 1) = RightBC;
        for (int i = 0; i < size; i++)
        {
            if(i!=0&&i!=(size-1)){
                ans(i)=0;
            }
            auto func = [i, this, funcbase](double x)
            {
                return f->operator()(x) * funcbase[i](x);
            };
            if (i >= 1)
            {
                ans(i) += Numerical_Integral(nowgrid[i - 1], nowgrid[i], func);
            }
            if (i < size - 1)
            {
                ans(i) += Numerical_Integral(nowgrid[i], nowgrid[i + 1], func);
            }
        }
        break;
    }
};
#else
#endif