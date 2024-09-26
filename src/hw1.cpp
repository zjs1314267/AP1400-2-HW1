#include "hw1.h"
using namespace std;
Matrix algebra::zeros(size_t n, size_t m)
{
    return Matrix(n, vector<double>(m, 0));
}
Matrix algebra::ones(size_t n,size_t m)
{
    return Matrix(n, vector<double>(m, 1));
}

Matrix algebra::random(size_t n, size_t m, double min, double max){
    if (min > max){
        throw logic_error("min cannot be greater than max");
    }   
    random_device rd; //作种
    default_random_engine e(rd());//初始化随机引擎
    uniform_real_distribution<double> distr(min, max);//约束范围
    Matrix matrix;
    for (size_t i = 0; i < n; ++i) {
        vector<double> vec;
        for (size_t j = 0; j < m; ++j) {
            vec.push_back(distr(e));
        }
        matrix.push_back(vec);  
        }
    return matrix;
}
void algebra::show(const Matrix& matrix)
{
    for(const auto row :matrix){
        for(const double element : row){
                cout<<setw(10)<< fixed<<setprecision(3)<<element;
        }
        cout <<'\n';
    }
}
Matrix algebra::multiply(const Matrix& matrix, double c)
{   
    Matrix A;
    for(const auto row :matrix){
        vector<double> Arow;
        for(const double element : row ){
            Arow.push_back(element*c);
        }
        A.push_back(Arow);
    }
    return A;
}
Matrix algebra::multiply(const Matrix& matrix1,const Matrix& matrix2)
{
    if(matrix1.empty()|matrix2.empty())
    return Matrix{};
    if(matrix1[0].size() != matrix2.size())
    throw logic_error("matrices with wrong dimensions cannot be multiplied");
   
    Matrix A=zeros(matrix1.size(),matrix2[0].size());

    for(int i=0;i<matrix1.size();i++){
        for(int j=0;j<matrix2[0].size();j++){
            double B=0;
            for(int m=0;m<matrix1[0].size();m++)
            {
                B+=matrix1[i][m]*matrix2[m][j];
            }
            A[i][j]=B;
        }
    }
    return A;
}
Matrix algebra::sum(const Matrix& matrix, double c)
{
    Matrix A=matrix;
    for(auto &row:A){
        for(auto &element :row)
        {
            element+=c;
        }
    }
    return A;
}
Matrix algebra::sum(const Matrix& matrix1, const Matrix& matrix2)
{

    if(matrix1.empty()&&matrix2.empty())
    return Matrix{};
    if(matrix1.size() != matrix2.size()||matrix1[0].size()!=matrix2[0].size()||matrix1.empty()||matrix2.empty())
    throw logic_error("matrices with wrong dimensions cannot be added");

    Matrix sum=move(matrix1);
    for(int i=0;i<matrix2.size();++i)
    for(int j=0; j<matrix2[0].size();++j)
    sum[i][j]+=matrix2[i][j];

    return sum;
}
Matrix algebra::transpose(const Matrix& matrix)
{
    if(matrix.empty())
    return Matrix{};

    Matrix trans=algebra::zeros(matrix[0].size(), matrix.size());
    for(int i=0; i<matrix[0].size();++i)
    for(int j=0; j<matrix.size();j++)
    trans[i][j]=matrix[j][i];

    return trans;
}
Matrix algebra::minor(const Matrix& matrix, size_t n, size_t m)
{
    Matrix minor;
    int rownum=0;
    for(auto row:matrix)
    {
        int colnum=0;
        vector<double> rowvec;
         if(rownum ==n) {
            rownum++;
            continue;
            }
        else 
        {
            for(auto element:row)
            {
                if(colnum==m)
                {
                    colnum++;
                    continue;
                }
                else {
                    rowvec.push_back(element);
                    colnum++;
                }
            }
            rownum++;
            minor.push_back(rowvec);
        }
    }  
    return minor;
}
double algebra::determinant(const Matrix& matrix)
{
    if(matrix.empty())
    return 1;
    if(matrix.size()!=matrix[0].size())
    throw logic_error("bad math");
    if(matrix.size()==1)
    return matrix[0][0];
    double count=0;
    for(int i=0;i<matrix[0].size();i++)
    {
        count+=pow(-1,i)*matrix[0][i]*determinant(minor(matrix,0,i));
    }
    return count;
}
Matrix algebra::inverse(const Matrix& matrix)
{   
    if(matrix.empty())
    return Matrix{};
    if(matrix.size()!=matrix[0].size())
    throw logic_error("bad math");
    double deter=algebra::determinant(matrix);
    if(deter==0)
    throw logic_error("bad math");
    Matrix inverse(matrix);
    for(int i=0;i<matrix.size();i++){
        for(int j=0 ; j<matrix[0].size();j++){
            inverse[i][j]=determinant(minor(matrix,j,i))/deter;
        }
    }
    return inverse;
}
Matrix algebra::concatenate(const Matrix& matrix1, const Matrix& matrix2, int axis)
{
    Matrix conmatrix(matrix1);
    if(axis==0){
        if(matrix1[0].size()!=matrix2[0].size())
        throw logic_error("bad math");
        else{
            conmatrix.insert(conmatrix.end(),matrix2.begin(),matrix2.end());
        }
    }
    if(axis==1){
        if(matrix1.size()!=matrix2.size())
        throw logic_error("bad math");
        else{
            for(int i=0;i<matrix1.size();i++)
            conmatrix[i].insert(conmatrix[i].end(),matrix2[i].begin(),matrix2[i].end());
        }
    }
    return conmatrix;
}
Matrix algebra::ero_swap(const Matrix& matrix, size_t r1, size_t r2)
{
    if(r1>matrix.size()-1||r2>matrix.size()-1)
        throw logic_error("are you kidding");
    Matrix swapmatrix=move(matrix);
    swap(swapmatrix[r1],swapmatrix[r2]);
    return swapmatrix;
}
Matrix algebra::ero_multiply(const Matrix& matrix, size_t r, double c)
{
    Matrix eromul =move(matrix);
    for(int i=0;i<eromul[0].size();i++)
    {
        eromul[r][i]*=c;
    }
    return eromul;
}
Matrix algebra::ero_sum(const Matrix& matrix,size_t r1, double c,size_t r2)
{
    Matrix erosum = move(matrix);
    for(int i=0;i<erosum[0].size();i++)
    {
        erosum[r2][i]+=erosum[r1][i]*c;
    }
    return erosum;
}
Matrix algebra::upper_triangular(const Matrix& matrix)
{
    if(matrix.empty())
    return Matrix{};
    if(matrix.size()!=matrix[0].size())
    throw logic_error("bad math");
    Matrix tri = move(matrix);
    for(int i=0;i<tri.size();i++)
    {
        while (tri[i][i]==0)
        {
            for(int n=i+1;n<tri.size();n++)
            {
                if(tri[n][i]!=0)
                {
                    tri=ero_swap(tri,i,n);
                    break;
                }
            }
        }
        for(int m=i+1;m<tri.size();m++)
        tri=ero_sum(tri,i, -tri[m][i]/tri[i][i],m);
    }
    return tri;
}
