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

const double P = 0.01;
const double b = 1.0;   
const int edim = 3;     // Количество узлов на элемент
const int ndim = 2;     // Степень свободы

vector<int> node_down_index;
vector<int> node_left_index;
vector<double> node;
vector<int> elem;
vector<int> elements_on_down;
vector<int> elements_on_left;
vector<double> bound;
vector<double> load;
vector<double> elastic = { 20.6, 0.3 }; // [MN/cm^2], [1]
double thickness = 1.0;

void loadNodes(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading nodes");
    }

    std::string line;
    int count_nodes = 0;
    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), ','), line.end());
        std::istringstream line_stream(line);
        double value0, value1;
        line_stream >> value0 >> value1;
        node.push_back(value0);
        node.push_back(value1);
        if (value1 == 0)
            node_down_index.push_back(count_nodes);
        if (value0 == 0)
            node_left_index.push_back(count_nodes);
        count_nodes++;
    }
}

void loadElements(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open file for reading nodes");
    }

    std::string line;
    int count_elements = 0;
    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(), line.end(), ','), line.end());
        std::istringstream line_stream(line);
        double node0, node1, node2;
        line_stream >> node0 >> node1 >> node2;

        bool found_down_0 = std::binary_search(node_down_index.begin(), node_down_index.end(), node0);
        bool found_down_1 = std::binary_search(node_down_index.begin(), node_down_index.end(), node1);
        bool found_down_2 = std::binary_search(node_down_index.begin(), node_down_index.end(), node2);

        if (found_down_0 == true || found_down_1 == true || found_down_2 == true)
            elements_on_down.push_back(count_elements);

        bool found_left_0 = std::binary_search(node_left_index.begin(), node_left_index.end(), node0);
        bool found_left_1 = std::binary_search(node_left_index.begin(), node_left_index.end(), node1);
        bool found_left_2 = std::binary_search(node_left_index.begin(), node_left_index.end(), node2);

        if (found_left_0 == true || found_left_1 == true || found_left_2 == true)
            elements_on_left.push_back(count_elements);

        elem.push_back(node0);
        elem.push_back(node1);
        elem.push_back(node2);
        count_elements++;
    }
}

void setBoundaryConditions() {
    int nsize = node.size() / ndim;
    bound.resize(nsize * ndim, std::numeric_limits<double>::infinity());
    load.resize(nsize * ndim, 0.0);

    int numRightNodes = 0;

    for (int i = 0; i < nsize; ++i) {
        double x = node[i * ndim];
        if (x == b) {
            numRightNodes++;
        }
    }

    double loadPerNode = P / numRightNodes;
    for (int i = 0; i < nsize; ++i) {
        double x = node[i * ndim];
        double y = node[i * ndim + 1];

        if (x == 0.0) bound[i * ndim] = 0.0;
        if (y == 0.0) bound[i * ndim + 1] = 0.0;

        if (x == b) {
            load[i * ndim] = loadPerNode; 
        }
    }
}

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
    for (size_t i = 0; i < a.size(); i++)
        det += a[i];

}

void blockStressStiffnessPE(double det, double bp, double bq, double cp, double cq, vector<double>& block)
{

    double nu = elastic[1];
    double E = elastic[0];
    double lambda = (nu * E) / ((1 + nu) * (1 - 2 * nu));

    double mu = E / (2 + 2 * nu);
    block[0 * ndim + 0] = ((lambda + 2 * mu) * bp * bq + mu * cp * cq) / (2 * det);

    block[0 * ndim + 1] = (lambda * bp * cq + mu * cp * bq) / (2 * det);

    block[1 * ndim + 0] = (lambda * cp * bq + mu * bp * cq) / (2 * det);
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
            for (size_t iterationdown = 1; iterationdown < countiteratiodown + 1; iterationdown++) //итерация по позициям под диагональю
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

void calculateStrains(int esize, const vector<double>& displacements, const vector<int>& elem, vector<long double>& strains) {
    for (int ie = 0; ie < esize; ie++) {
        
        vector<int> nodes = { elem[ie * edim], elem[ie * edim + 1], elem[ie * edim + 2] };

        vector<double> b(3), c(3), area(3);
        gradientB(nodes, b);
        gradientC(nodes, c);
        gradientA(nodes, area);

        double totalArea = abs(area[0] + area[1] + area[2]) / 2.0;

        long double exx = (b[0] * displacements[nodes[0] * ndim] + b[1] * displacements[nodes[1] * ndim] + b[2] * displacements[nodes[2] * ndim]) / (2 * totalArea);
        long double eyy = (c[0] * displacements[nodes[0] * ndim + 1] + c[1] * displacements[nodes[1] * ndim + 1] + c[2] * displacements[nodes[2] * ndim + 1]) / (2 * totalArea);
        long double exy = ((b[0] * displacements[nodes[0] * ndim + 1] + b[1] * displacements[nodes[1] * ndim + 1] + b[2] * displacements[nodes[2] * ndim + 1]) +
            (c[0] * displacements[nodes[0] * ndim] + c[1] * displacements[nodes[1] * ndim] + c[2] * displacements[nodes[2] * ndim])) / (4 * totalArea);

        strains.push_back(exx);
        strains.push_back(eyy);
        strains.push_back(exy);
    }
}

void calculateStresses(int esize, const vector<long double>& strains, vector<long double>& stresses) {

    double E = elastic[0];
    double nu = elastic[1];
    double factor = E / ((1 + nu) * (1 - 2 * nu));

    for (int i = 0; i < esize; i++) {
        
        long double exx = strains[i * 3];
        long double eyy = strains[i * 3 + 1];
        long double exy = strains[i * 3 + 2];

        long double sxx = factor * ((1 - nu) * exx + nu * eyy);
        long double syy = factor * (nu * exx + (1 - nu) * eyy);
        long double sxy = factor * (1 - 2 * nu) / 2.0 * exy;

        stresses.push_back(sxx);
        stresses.push_back(syy);
        stresses.push_back(sxy);
    }
}

void StressesDownPrintFile(const std::vector<long double>& stresses) {
    std::ofstream outfile("output_down.txt");
    for (const auto& element : elements_on_down) {
        
        int stress_index = element * edim;
        if (stress_index >= stresses.size()) {
            continue;
        }

        int element_index = element * edim;
        if (element_index >= elem.size()) {
            continue;
        }

        for (int i = 0; i < edim; ++i) {
            int node_index = elem[element_index + i];
            if (std::find(node_down_index.begin(), node_down_index.end(), node_index) != node_down_index.end()) {
                
                long double x = node[node_index * 2];

                outfile << std::scientific << std::setprecision(3)
                    << x << " " << stresses[stress_index] << '\n';
            }
        }
    }

    outfile.close();
}

void StressesLeftPrintFile(const std::vector<long double>& stresses) {
    std::ofstream outfile("output_left.txt");

    for (const auto& element : elements_on_left) {
        
        int stress_index = element * edim;
        if (stress_index >= stresses.size()) {
            continue; 
        }

        int element_index = element * edim;
        if (element_index >= elem.size()) {
            continue;
        }

        for (int i = 0; i < edim; ++i) {
            int node_index = elem[element_index + i];
            
            if (node[node_index * 2] == 0) {
                long double y = node[node_index * 2 + 1];

                outfile << std::scientific << std::setprecision(3)
                    << y << " " << stresses[stress_index] << '\n';
            }
        }
    }

    outfile.close();
}

int main(int argc, char** argv)
{
    const std::string filenodes = "nodes.txt";
    const std::string fileelem = "elements.txt";
    loadNodes(filenodes);
    loadElements(fileelem);
    int esize = elem.size() / edim;// sizeof(elem) / sizeof(int) / edim;
    int nsize = node.size() / ndim;
    int N = nsize * ndim;
    setBoundaryConditions();
    vector<double> globalStiffness(N * N);
    for (int i = 0; i < N * N; i++)
    {
        globalStiffness[i] = 0;
    }

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

        for (int p = 0; p < edim; p++)
        {

            for (int q = 0; q < edim; q++)
            {
                blockStressStiffnessPE(det, b[p], b[q], c[p], c[q], block);             
                for (int i = 0; i < ndim; i++)
                {
                    for (int j = 0; j < ndim; j++)
                    {                       
                       globalStiffness[(P[p] * ndim + i) * N + (Q[q] * ndim + j)] += block[i * ndim + j];
                    }
                }

            }
        }

        
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

    applyKinematic(N, bound, load, globalStiffness);
    
    vector<double> val;
    vector<int> I;
    vector<int> J;
    int nz;

    transformMatrixToCsr(N, N, nz, globalStiffness, val, I, J);
    printSCR(N, N, nz, val, I, J);

    vector<double> displacement(N);
    for (int i = 0; i < N; i++)
    {
        displacement[i] = 0.0;
    }

    status = cgTest();

    status = cg_host(N, nz, I, J, val, load, displacement, 1);
    resultPrintFile(nsize, ndim, node, displacement);

    vector<long double> strains;
    calculateStrains(esize, displacement, elem, strains);

    vector<long double> stresses;
    calculateStresses(esize, strains, stresses);

    StressesDownPrintFile(stresses);
    StressesLeftPrintFile(stresses);

    std::cin.get();
    int exit(status);
}