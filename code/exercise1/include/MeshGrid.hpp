/**
 * @file meshgrid.hpp
 * @author Shuang Hu (22135017@mail.zju.edu.cn)
 * @brief Mesh on 1-dimensional space [0,1]
 * @version 0.1
 * @date 2022-09-20
 */
#include<vector>
#include<assert.h>
/**
 * @brief Mesh grids in 1-dimensional space.
 * 
 * @details 
 */
class Mesh1D{
protected:
    bool uniform;//uniform grid or not
    int size;//number of small interval
    double LeftSide;
    double RightSide;
    std::vector<double> meshgrid;//\{x_{i}\},except for the Left-side and Right-side.
public:
    /**
     * @brief Construct a new Mesh 1 D object
     * 
     * @details Default constructor
     */
    Mesh1D() = default;
    /**
     * @brief Construct a new Mesh 1 D object
     * @param [in] nsize number of interval pieces
     * @param [in] Left left side
     * @param [in] Right right side
     * @param [in] isuniform uniform grid or not
     * @param [in] mesh if the grid isn't uniform, input the grids.
     * 
     * @details 
     */
    Mesh1D(int nsize,double Left,double Right,bool isuniform=true,const std::vector<double>& mesh={}){
        size=nsize;
        assert(isuniform||(mesh.size()==nsize+1));
        LeftSide=Left;
        RightSide=Right;
        uniform=isuniform;
        if(isuniform){
            double total_len=Right-Left;
            double dx=total_len/size;
            for(int i=0;i<=size;i++){
                double pos=Left+dx*i;
                meshgrid.push_back(pos);
            }            
        }
        else{
            for(auto c:mesh){
                meshgrid.push_back(c);
            }
        }
    }
};