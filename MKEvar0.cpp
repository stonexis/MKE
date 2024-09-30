#include <iostream>
#include <fstream>
#include <sstream>
#include <array>
#include <vector>
#include <string>
#include <cmath>
#include <math.h>
#include <algorithm>
#include <cctype>
#include <functional> // std::plus // isinf()

#include <stdlib.h>

#include "Header.h"

using namespace std;

const int edim = 3; //количество узлов
const int ndim = 2; //степень свободы

const vector<int> elem =
{
  0,  1,  3,
  1,  4,  3,
  0,  3,  2,
};
const int esize = elem.size() / edim;// sizeof(elem) / sizeof(int) / edim;

vector<double>  node = // [cm]
{
  -4.0,   0.0,
  -2.0,   0.0,
  -5.0,   1.7,
  -3.0,   1.7,
  -1.0,   1.7,
};
const int nsize = node.size() / ndim;

const vector<double> elastic = { 20.6, 0.3 }; // { [MN/cm2], [1]}

const double thickness = 1.0; // [cm]

vector<double> load = // [MN]
{
    0.0,
    0.0,
    0.0,
    0.0,
    0.0,
   -3.0e-3,
    0.0,
   -3.0e-3,
    0.0,
    0.0,
};

vector<double> bound = { // [cm] 
    0.0,
    0.0,
    0.0,
    0.0,
    INFINITY,
    INFINITY,
    INFINITY,
    INFINITY,
    INFINITY,
    INFINITY,
};
//bound[ndim * 2 + 0] = ;

void gradientA(vector<int> P, vector<double>& result)
{ // ai = eijm xj ym; a0 = x1 y2 - x2 y1; a1 = x2 y0 - x0 y2; a2 = x0 y1 - x1 y0;
    result[0] = node[P[1] * ndim + 0] * node[P[2] * ndim + 1] - node[P[2] * ndim + 0] * node[P[1] * ndim + 1];
    result[1] = node[P[2] * ndim + 0] * node[P[0] * ndim + 1] - node[P[0] * ndim + 0] * node[P[2] * ndim + 1];
    result[2] = node[P[0] * ndim + 0] * node[P[1] * ndim + 1] - node[P[1] * ndim + 0] * node[P[0] * ndim + 1];

}

void gradientB(vector<int> P, vector<double>& result)
{ // bj = eijm 1i ym; b0 = y1 - y2; b1 = y2 - y0; b2 = y0 - y1;
    result[0] = node[P[1] * ndim + 1] - node[P[2] * ndim + 1];
    result[1] = node[P[2] * ndim + 1] - node[P[0] * ndim + 1];
    result[2] = node[P[0] * ndim + 1] - node[P[1] * ndim + 1];
}

void gradientC(vector<int> P, vector<double>& result)
{ // cm = eijm 1i xj; c0 = x2 - x1; c1 = x0 - x2; c2 = x1 - x0;

    result[0] = node[P[2] * ndim + 0] - node[P[1] * ndim + 0];
    result[1] = node[P[0] * ndim + 0] - node[P[2] * ndim + 0];
    result[2] = node[P[1] * ndim + 0] - node[P[0] * ndim + 0];
}

void getDet(vector<double> a, double& det)
{      
       for (size_t i=0;i<a.size();i++)
           det+= a[i];
         
}

void blockStressStiffnessPE(double det, double bp, double bq, double cp, double cq, vector<double>& block)
{

    double nu = elastic[1];
    double E = elastic[0];
    double lambda = (nu*E)/((1+nu)*(1-2*nu));
    
    double mu = E / (2 + 2 * nu);
    block[0 * ndim + 0] = ((lambda+2*mu)*bp*bq + mu*cp*cq ) / (2 * det);
    
    block[0 * ndim + 1] = (lambda*bp*cq + mu*cp*bq) / (2*det);

    block[1 * ndim + 0] = (lambda*cp*bq + mu*bp*cq) / (2 * det);
    block[1 * ndim + 1] = ((lambda + 2 * mu) * cp * cq + mu * bp * bq) / (2 * det);
}

void blockStressStiffnessPS(double det, double bp, double bq, double cp, double cq, vector<double>& block) //??
{
    double v = elastic[1];
    double E = elastic[0];

    //   block[0 * ndim + 0] = ;
    //   block[0 * ndim + 1] = ;

    //   block[1 * ndim + 0] = ;
    //   block[1 * ndim + 1] = ;
}

void applyKinematic(int N, vector<double> u, vector<double>& f, vector<double>& stiffness)
{   
    vector<int> known;
    int pos_on_diag;

    for (size_t i = 0; i < N; i++)
    {
        if (u[i] != INFINITY)
          known.push_back(i);
    }
    for (size_t i = 0; i < known.size(); i++)
    {
        pos_on_diag = N * known[i] + known[i];

        for (size_t i = 1; i < N - pos_on_diag % N; i++)
        {
            stiffness[pos_on_diag + i] = 0;
        }
        for (size_t i = 1; i < pos_on_diag % N + 1; i++)
        {
            stiffness[pos_on_diag - i] = 0;
        }
        f[known[i]] = stiffness[pos_on_diag] * u[known[i]];
    }

    for (size_t i = 0; i < known.size(); i++)
    {
        pos_on_diag = N * known[i] + known[i];
        int countiteratiodown = N - pos_on_diag / N - 1;
        int diagpos_appropriate_f_down, diagpos_appropriate_f_app;
        if (pos_on_diag != stiffness.size() - 1)
        {
            for (size_t iterationdown = 1; iterationdown < countiteratiodown+1; iterationdown++) //итерация по позициям под диагональю
            {
                diagpos_appropriate_f_down = (pos_on_diag + iterationdown * N) / N;
                f[diagpos_appropriate_f_down] = f[diagpos_appropriate_f_down] - stiffness[pos_on_diag + iterationdown * N] * u[known[i]];
                stiffness[pos_on_diag + iterationdown * N] = 0;
            }
        }

        int countiteratioapp = pos_on_diag % N;
        if (pos_on_diag != 0) {
            for (size_t iterationapp = 1; iterationapp < countiteratioapp + 1; iterationapp++) //итерация по позициям над диагональю
            {
                diagpos_appropriate_f_app = (pos_on_diag - iterationapp * N) / N;
                f[diagpos_appropriate_f_app] = f[diagpos_appropriate_f_app] - stiffness[pos_on_diag - iterationapp * N] * u[known[i]];
                stiffness[pos_on_diag - iterationapp * N] = 0;
            }
        }

    }




  
    printf("\nKinematic boundary condition is applyed\n\n");
}

int kinematicTest()
{ // Л.Сегерлинд. Применение метода конечных элементов стр. 110-112
    vector<double> u = {
        150,
        INFINITY,
        INFINITY,
        INFINITY,
        40,
    };
    vector<double> f = {
        500,
        2000,
        1000,
        2000,
        900,
    };
    vector<double> stiffness = {
      55, -46,   4,   0,   0,
      -46, 140, -46,   0,   0,
        4, -46, 110, -46,   4,
        0,   0, -46, 142, -46,
        0,   0,   4, -46,  65,
    };
    const int N = u.size();
    printStiffnessF(N, stiffness, f);
    applyKinematic(N, u, f, stiffness);
    printStiffnessF(N, stiffness, f);

    double sum =
        stiffness[0 * N + 0] - 55.0 + f[0] - 8250 +
        stiffness[1 * N + 1] - 140.0 + stiffness[1 * N + 2] + 46.0 + f[1] - 8900 +
        stiffness[2 * N + 1] + 46.0 + stiffness[2 * N + 2] - 110.0 + stiffness[2 * N + 3] + 46.0 + f[2] - 240 +
        stiffness[3 * N + 2] + 46.0 + stiffness[3 * N + 3] - 142.0 + f[3] - 3840 +
        stiffness[4 * N + 4] - 65.0 + f[4] - 2600;
    return sum == 0.0 ? 0 : 1;
}

int sumStiffnessTest(int N, vector<double> stiffness, vector<double> f)
{
    const double tol = 1e-9;
    double sum = 0;

    for (size_t i = 0; i < stiffness.size(); i++)
        sum += stiffness[i];
    cout << "(" << sum << ")";
    printf("\n\nSum stiffness is %5.9e\n\n", sum);
    return abs(sum) < tol ? 0 : 1;
}

// Преобразовать вектор-матрицу жёсткости в формат CSR https://docs.nvidia.com/cuda/cusparse/index.html#csr-format
void transformMatrixToCsr(int M, int N, int& nz, vector<double>& stiffness, vector<double>& _val, vector<int>& _I, vector<int>& _J)
{

    vector<double> value;
    vector<int> column;
    vector<int> row;
    nz = 0; //количество ненулевых элементов
    bool flag;
    for (size_t i = 0; i < stiffness.size(); i++)
    {
        if (i % N == 0)
        {
            flag = true;
        }
        if (stiffness[i] != 0)
        {
            value.push_back(stiffness[i]);
            column.push_back(i % N);
            if (flag == true)
            {
                row.push_back(nz);
                flag = false;
            }
            nz++;
        }

    }
    row.push_back(nz);
    _I.resize(M + 1);
    _J.resize(nz);
    _val.resize(nz);

    for (int i = 0; i < nz; i++)
    {
        if (i < M + 1)
        {
            _I[i] = row[i];
        }
        _val[i] = value[i];
        _J[i] = column[i];
    }
}

void transformMatrixTest()
{
    const int M = 4;
    const int N = 5;

    vector<double> stiffness = {
      1.0, 4.0, 0.0, 0.0, 0.0,
      0.0, 2.0, 3.0, 0.0, 0.0,
      5.0, 0.0, 0.0, 7.0, 8.0,
      0.0, 0.0, 9.0, 0.0, 6.0,
    };

    vector<double> val;
    vector<int> I;
    vector<int> J;
    int nz;
    transformMatrixToCsr(M, N, nz, stiffness, val, I, J);
    printSCR(M, N, nz, val, I, J);
}

int main(int argc, char** argv)
{
    int N = nsize * ndim;

    vector<double> globalStiffness(N * N);
    for (int i = 0; i < N * N; i++)
    {
        globalStiffness[i] = 0;
    }

    printStiffness(N, globalStiffness);

    vector<double> a{ 0,0,0 }, b{ 0,0,0 }, c{ 0,0,0 };
    vector<int> P{ 0,0,0 }, Q{ 0,0,0 };
    vector<double> block(ndim * ndim);

    for (int ie = 0; ie < esize; ie++)
    {
        P[0] = elem[ie * edim + 0];
       
        P[1] = elem[ie * edim + 1];
        
        P[2] = elem[ie * edim + 2];
        
        Q[0] = P[0];
        
         Q[1] = P[1];
         
         Q[2] = P[2];
        

         double det = 0;
         gradientA(P, a);
         gradientB(P, b);
         gradientC(P, c);
         getDet(a, det);

         printABCDet(ie, edim, det, a, b, c);

         //vector<double> elemStiffness[ndim * edim * ndim * edim];

        for (int p = 0; p < edim; p++)
        {

            for (int q = 0; q < edim; q++)
            {
                 blockStressStiffnessPE(det, b[p], b[q], c[p], c[q], block);

                 printf("\nPQ[%d%d] pq[%d%d]\n", P[p], Q[q], p, q);
                 printStiffness(ndim, block);

                for (int i = 0; i < ndim; i++)
                {
                    for (int j = 0; j < ndim; j++)
                    {
                         //elemStiffness[(p * ndim + i) * edim * ndim + (q * ndim + j)] = block[i * ndim + j]; 
                        globalStiffness[(P[p] * ndim + i) * N + (Q[q] * ndim + j)] += block[i * ndim + j];
                    }
                }

            }
        }

         //printStiffness(ndim * edim, elemStiffness);
    }

    int status = sumStiffnessTest(N, globalStiffness, load);

   if (status != 0)
    {
        printf("\nSum stiffness error\n");
        int exit(status);
    }
   else
    {
        printf("\nSum stiffness ok\n");
    }

       printStiffnessF(N, globalStiffness, load);
       applyKinematic(N, bound, load, globalStiffness);
       printStiffnessF(N, globalStiffness, load);

       status = kinematicTest();

    if (status != 0)
    {
        printf("\nKinematic error\n");
        int exit(status);
    }
    else
   {
       printf("\nKinematic test ok\n");
    }

       transformMatrixTest();

    vector<double> val;
    vector<int> I;
    vector<int> J;
    int nz;

       transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
       printSCR(N, N, nz, val, I, J);

    vector<double> result(N);
    for (int i = 0; i < N; i++)
    {
        result[i] = 0.0;
    }

    status = cgTest();

     status = cg_host(N, nz, I, J, val, load, result, 1);

    resultPrint(nsize, ndim, node, result);

    std::cin.get();
    int exit(status);
}