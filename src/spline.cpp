#include <iostream>
#include <math.h>
#include <vector>
#include <kdl/velocityprofile_trap.hpp>
#include <kdl/rotational_interpolation_sa.hpp>
#include <kdl/path_roundedcomposite.hpp>
#include <kdl/trajectory_segment.hpp>
#include <Eigen/Dense>
#include <fstream>

#define ROBOT_DOF_SIZE 7
#define MAX_P 4

void BasisFuns(int i, double u, int p, double U[], double B[]);
int WhichSpan(double u, double U[], int n_knot, int p);
void BSplinePoint( double u, double U[], int n_knot, int p, std::vector<Eigen::Vector3f> P,int d, double s[]);

void solve(double* a, double* b, double* c, double* d, int n);

using namespace std;
int main()
{

    for(uint i = 0; i < _n; i++){


        // at the endpoint previous span intervall will be used again by setting u − eps for u == u_max

        int span = WhichSpan(_U_[i], _U, _k, _p);
        cout << "span: " << span << endl;


        double _B[4];
        // BasisFuns(span, _U[4], _p, _U, _B); // wrong
        BasisFuns(span, _U_[i], _p, _U, _B);

        for (uint i = 0; i < 4; i++)
            cout << _B[i] << " ";

        cout << endl;
    }

    return 0;

    }

    int p = 3, dim = 3; // p:polynomal degree, dim:cart space dimension

    // cartesian trajectory in x,y,z
    std::vector<Eigen::Vector3f> q;
    q.push_back(Eigen::Vector3f(83, -54, 119));
    q.push_back(Eigen::Vector3f(-64, 10, 124));
    q.push_back(Eigen::Vector3f(42, 79, 226));
    q.push_back(Eigen::Vector3f(-98, 23, 222));
    q.push_back(Eigen::Vector3f(-13, 125, 102));
    q.push_back(Eigen::Vector3f(140, 81, 92));
    q.push_back(Eigen::Vector3f(43, 32, 92));
    q.push_back(Eigen::Vector3f(-65, -17, 143));
    q.push_back(Eigen::Vector3f(-45, -89 , 182));
    q.push_back(Eigen::Vector3f(71, 90 , 192));

    int n = q.size()-1; // to use as index for q
    int n_knot = n + 6;
    cout << "n " << n << endl;


    // knot vector u
    double U[n_knot] = {0, 0, 0, 0, 0.11, 0.23, 0.35, 0.48, 0.60, 0.68, 0.77, 0.84, 1.0, 1.0, 1.0, 1.0};
    // int nr_of_u = sizeof(_U)/sizeof(_U[0]);

    double _U[n]; // remove first and last 3 values from knot vector u
    for (uint i = 0; i < n; i++)
        _U[i] = U[i+3];



//    double P[16*3] = {83, 34, -168, 146, -182, -45, 207, 31, -89, -29, 32, 71,
//                      -5, 128, -41, 177, 88, 21, 14, -172, 30, 90,-54, -32,
//                      119, 120,
//                      88, 252, 245, 68, 98, 83, 121, 218, 188, 192};

    //waypoints q
//    double q[10*d] =    {83, -64 ,42 ,-98, -13 ,140, 43, -65 ,-45, 71,
//                      -54 ,10 ,79 ,23 ,125 ,81 ,32 ,-17 ,-89, 90,
//                     119 ,124 ,226 ,222 ,102 ,92 ,92 ,134 ,182 ,192};



    // derivatives at the endpoints
    Eigen::Vector3f t_0 = Eigen::Vector3f(-1236, 538, 42);
    Eigen::Vector3f t_n = Eigen::Vector3f(732, 1130, 63);

    // control points
    std::vector<Eigen::Vector3f> P;

    P.push_back(q.at(0));
    //Eigen::Vector3f p_1 = q.at(0) + ((U[4]/3.0)*t_0);
    P.push_back(q.at(0) + ((U[4]/3.0)*t_0));

    // vectors a, b, c, d for tridiagonal matrix PB = R
    double a[n-1], b[n-1], c[n-1], d_x[n-1], d_y[n-1], d_z[n-1];
    std::vector<Eigen::Vector3f> d;

    a[0] = 0;
    c[n-2] = 0;
    for (int i = 1; i <= n-1; i++){

        // Eigen::VectorXf B_a(4);
        double B[MAX_P];

        int intervall = WhichSpan(_U[i], U, n_knot, p);
        BasisFuns(intervall, _U[i], p, U, B); //a

        cout << "intervall " << intervall << ", B1: " << B[0]<< ", B2: " << B[1]<< ", B3: " << B[2]<< ", B4: " << B[3] << endl;

        if (i > 1)
            a[i-1] = B[0];

        b[i-1] = B[1];

        if (i <= n-1)
            c[i-1] = B[2];


        if (i == 1){

            d_x[i-1] = (q.at(1) - B[0]* P.at(1))(0);
            d_y[i-1] = (q.at(1) - B[0]* P.at(1))(1);
            d_z[i-1] = (q.at(1) - B[0]* P.at(1))(2);

        }
        else if (i == n-1){

            d_x[i-1] = (q.at(i) - B[2] * (q.at(n) - ((1.0-U[n+2])/3.0)*t_n))(0); // TODO double check [n+2]
            d_y[i-1] = (q.at(i) - B[2] * (q.at(n) - ((1.0-U[n+2])/3.0)*t_n))(1);
            d_z[i-1] = (q.at(i) - B[2] * (q.at(n) - ((1.0-U[n+2])/3.0)*t_n))(2);

        }

        else if (i < n-1){

            d_x[i-1] = (q.at(i))(0);
            d_y[i-1] = (q.at(i))(1);
            d_z[i-1] = (q.at(i))(2);
        }


    solve(a,b,c,d_x,n-1);
    solve(a,b,c,d_y,n-1);
    solve(a,b,c,d_z,n-1);

    for (int i = 0; i < n-1; i++) {
        cout << "[TEST] " <<  d_x[i] << endl;
     }


    for (uint i = 0; i < n-1; i++){
        P.push_back(Eigen::Vector3f(d_x[i], d_y[i], d_z[i]));
    }

    // add only after P_2 ... P_n were added
    P.push_back(q.at(n) - ((1.0-U[n+2])/3.0)*t_n); // TODO double check [n+2]
    P.push_back(q.at(n));
    //Eigen::Vector3f p_1 = q.at(0) + ((U[4]/3.0)*t_0);


    for (double u = 0; u <=1; u+=0.01)
    {
        double s[p];
        BSplinePoint(u, U, n_knot, p, P, dim, s);
        //for (uint i = 0; i < p; i++)
             cout << s[2] <<  ", " ;
        cout << endl;

    }
    return 0;
}


void solve(double* a, double* b, double* c, double* d, int n) {

    // n is the number of unknowns
    double c_tmp[n];
    for (uint i = 0; i < n; i++)
        c_tmp[i] = c[i];

    n--; // since we start from x0 (not x1)

    // make a copy of c, so c gets not changed by the algorithm
    c_tmp[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c_tmp[i] /= b[i] - a[i]*c_tmp[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c_tmp[i-1]);
    } 

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c_tmp[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c_tmp[i]*d[i+1];
    }
}



void BSplinePoint( double u, double U[], int n_knot, int p,
std::vector<Eigen::Vector3f> P,int d, double s[])
/*
Inputs:u
- value of the independent variable
U[]
- Knot vector
n_knot - length of U[] -1
p
- degree of the spline
P[]
- Control point vector
d
- dimensions of a control point (2 in 2D, 3 in 3D, etc.)
Output:s[]
- value of the B-spline at uthe
*/
{
    double B[MAX_P];
    int i, k, j;
    i = WhichSpan(u, U, n_knot, p);
    BasisFuns(i, u, p, U, B);
    for (k = 0; k < d; k++) /* For each components of the B-spline*/
    {
        s[k] = 0;
        for (j = 0; j <= p; j++)
        {
            // s[k] = s[k] + P[k * (n_knot - p) + i - p + j] * B[j];
            s[k] = s[k] + P.at(i-p+j)(k) * B[j];
        }
    }
}



void BasisFuns(int i, double u, int p, double U[], double B[])
/*
Input:
i
- knot span including u
u
- value of the independent variable
p
- degree of the spline
U[] - Knot vector
Output:
B[] - value of the nonvanishing basis function at u
*/
{
    int j, r;
    double temp, acc;
    double DR[MAX_P], DL[MAX_P];
    B[0] = 1;
    for (j = 1; j <= p; j++)
    {
        DL[j] = u - U[i + 1 - j];
        DR[j] = U[i + j] - u;
        acc = 0.0;
        for (r = 0; r <= j - 1; r++)
        {
            temp = B[r] / (DR[r + 1] + DL[j - r]);
            B[r] = acc + DR[r + 1] * temp;
            acc = DL[j - r] * temp;
        }
        B[j] = acc;
    }
}


int WhichSpan(double u, double U[], int n_knot, int p)
/*
Input:
u
- value of the independent variable
U[]
- Knot vector
n_knot - length of U[] -1
p
- degree of the spline
Output: mid
- index of the knot span including u
*/
{

    int high, low, mid;
    high = n_knot - p;

    if (u == U[high])// TODO check if correct (p. 470 in Trajectory_Planning_for_Automatic_Machines_and_Robots.pdf)
        u -= 0.0000000000001;


    low = p;
    if (u == U[high])
        mid = high;
    else
    {
        mid = (high + low) / 2;
        while ((u < U[mid]) || (u >= U[mid + 1]))
        {
            if (u == U[mid + 1])
                mid = mid + 1;
            /* knot with multiplicity >1 */
            else
            {
                if (u > U[mid])
                    low = mid;
                else
                    high = mid;
                mid = (high + low) / 2;
            }
        }
    }
    return mid;
}




