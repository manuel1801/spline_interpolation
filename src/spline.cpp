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
void BSplinePoint(double u, double U[], int n_knot, int p, std::vector<Eigen::Vector3f> P, int d, double s[]);
void solve(double *a, double *b, double *c, double *d, int n);

using namespace std;

int main()
{
    int p = 3, dim = 3; // p:polynomal degree, dim:cart space dimension

    // cartesian trajectory in x,y,z
    std::vector<Eigen::Vector3f> q;
    q.push_back(Eigen::Vector3f(83.0, -54.0, 119.0));
    q.push_back(Eigen::Vector3f(-64.0, 10.0, 124.0));
    q.push_back(Eigen::Vector3f(42.0, 79.0, 226.0));
    q.push_back(Eigen::Vector3f(-98.0, 23.0, 222.0));
    q.push_back(Eigen::Vector3f(-13.0, 125.0, 102.0));
    q.push_back(Eigen::Vector3f(140.0, 81.0, 92.0));
    q.push_back(Eigen::Vector3f(43.0, 32.0, 92.0));
    q.push_back(Eigen::Vector3f(-65.0, -17.0, 143.0));
    q.push_back(Eigen::Vector3f(-45.0, -89.0, 182.0));
    q.push_back(Eigen::Vector3f(71.0, 90.0, 192.0));

    int n = q.size() - 1; // to use as index for q
    int n_knot = n + 6;
    cout << "n " << n << endl;

    // knot vector u
    double U[n_knot] = {0, 0, 0, 0, 0.11, 0.23, 0.35, 0.48, 0.60, 0.68, 0.77, 0.84, 1.0, 1.0, 1.0, 1.0};
    // int nr_of_u = sizeof(_U)/sizeof(_U[0]);

    double _U[n]; // remove first and last 3 values from knot vector u
                  //    for (uint i = 0; i < n; i++)
                  //        _U[i] = U[i+3];

    //    double P[16*3] = {83, 34, -168, 146, -182, -45, 207, 31, -89, -29, 32, 71,
    //                      -5, 128, -41, 177, 88, 21, 14, -172, 30, 90,
    //                 -54, -32, 119, 120,88, 252, 245, 68, 98, 83, 121, 218, 188, 192};

    //controll points from example
    //    {83.0,  34.0,   -168.0,   146.0, -182.0, -45.0,  207.0, 31.0, -89.0, -29.0,  32.0, 71.0,
    //    -54.0, -32.0,  -5.0,      128.0, -41.0,  177.0,  88.0,  21.0,  14.0, -172.0, 30.0,  90.0,
    //    ,119.0, 120.0, 88.0,      252.0, 245.0,  68.0,   98.0,  83.0,  121.0, 218.0, 188.0, 192.0};

    std::vector<Eigen::Vector3f> __P;
    __P.push_back(Eigen::Vector3f(83.0, -54.0, 119.0));
    __P.push_back(Eigen::Vector3f(34.0, -32.0, 120.0));
    __P.push_back(Eigen::Vector3f(-168.0, -5.0, 88.0));
    __P.push_back(Eigen::Vector3f(146.0, 128.0, 252.0));
    __P.push_back(Eigen::Vector3f(-182.0, - 41.0, 245.0));
    __P.push_back(Eigen::Vector3f(-45.0, 177.0, 68.0));
    __P.push_back(Eigen::Vector3f(207.0, 88.0, 98.0));
    __P.push_back(Eigen::Vector3f(31.0, 21.0, 83.0));
    __P.push_back(Eigen::Vector3f(-89.0, 14.0, 121.0));
    __P.push_back(Eigen::Vector3f(-29.0, -172.0, 218.0));
    __P.push_back(Eigen::Vector3f(32.0, 30.0, 188.0));
    __P.push_back(Eigen::Vector3f(71.0, 90.0, 192.0));

    //waypoints q
    //    double q[10*d] =    {83, -64 ,42 ,-98, -13 ,140, 43, -65 ,-45, 71,
    //                      -54 ,10 ,79 ,23 ,125 ,81 ,32 ,-17 ,-89, 90,
    //                     119 ,124 ,226 ,222 ,102 ,92 ,92 ,134 ,182 ,192};

    double __s[p];
    BSplinePoint(1.0, U, n_knot, p, __P, dim, __s);
    // vec.push_back(Eigen::Vector3f(s[0], s[1], s[2]));
    cout << "[point at u=1] " << __s[0] << " " << __s[1] << " " << __s[2] << endl;

    // return 0;

    // derivatives at the endpoints
    Eigen::Vector3f t_0 = Eigen::Vector3f(-1236, 538, 42);
    Eigen::Vector3f t_n = Eigen::Vector3f(732, 1130, 63);

    // control points
    std::vector<Eigen::Vector3f> P;

    P.push_back(q.at(0));
    P.push_back(q.at(0) + ((U[4] / 3.0) * t_0));

    // vectors a, b, c, d for tridiagonal matrix PB = R
    double a[n - 1], b[n - 1], c[n - 1], d_x[n - 1], d_y[n - 1], d_z[n - 1];
    std::vector<Eigen::Vector3f> d;

    a[0] = 0;
    c[n - 2] = 0;
    for (int i = 1; i <= n - 1; i++)
    {

        double B[MAX_P];
        int intervall = WhichSpan(U[i + 3], U, n_knot, p);
        BasisFuns(intervall, U[i + 3], p, U, B); //a
                                                 //        int intervall = WhichSpan(_U[i], U, n_knot, p);
                                                 //        BasisFuns(intervall, _U[i], p, U, B); //a

        cout << "intervall " << intervall << ", B1: " << B[0] << ", B2: " << B[1] << ", B3: " << B[2] << ", B4: " << B[3] << endl;

        if (i > 1)
            a[i - 1] = B[0];

        b[i - 1] = B[1];

        if (i <= n - 1)
            c[i - 1] = B[2];

        if (i == 1)
        {

            d_x[i - 1] = (q.at(1) - B[0] * P.at(1))(0);
            d_y[i - 1] = (q.at(1) - B[0] * P.at(1))(1);
            d_z[i - 1] = (q.at(1) - B[0] * P.at(1))(2);
        }
        else if (i == n - 1)
        {

            d_x[i - 1] = (q.at(i) - B[2] * (q.at(n) - ((1.0 - U[n + 2]) / 3.0) * t_n))(0); // TODO double check [n+2]
            d_y[i - 1] = (q.at(i) - B[2] * (q.at(n) - ((1.0 - U[n + 2]) / 3.0) * t_n))(1);
            d_z[i - 1] = (q.at(i) - B[2] * (q.at(n) - ((1.0 - U[n + 2]) / 3.0) * t_n))(2);
        }

        else if (i < n - 1)
        {

            d_x[i - 1] = (q.at(i))(0);
            d_y[i - 1] = (q.at(i))(1);
            d_z[i - 1] = (q.at(i))(2);
        }
    }

    // solve PB = R for P with tridiagonal matrix algorthm
    solve(a, b, c, d_x, n - 1);
    solve(a, b, c, d_y, n - 1);
    solve(a, b, c, d_z, n - 1);

    for (uint i = 0; i < n - 1; i++)
    {
        P.push_back(Eigen::Vector3f(d_x[i], d_y[i], d_z[i]));
    }

    P.push_back(q.at(n) - ((1.0 - U[n + 2]) / 3.0) * t_n); // TODO check if [n+2] correct
    P.push_back(q.at(n));

    for (auto p : P)
        cout << p(0) << " " << p(1) << " " << p(2) << endl;

    //    std::vector<Eigen::Vector3f> vec;
    std::ofstream file;
    file.open("/home/manuel/hrg/spline_interpolation/output/spline.txt");

    // output only last point
    // double s[p];
    // BSplinePoint(1.0, U, n_knot, p, P, dim, s);
    // vec.push_back(Eigen::Vector3f(s[0], s[1], s[2]));
    // cout << "[last point own spline: ]" << s[0] << " " << s[1] << " " << s[2] << endl;

    for (double u = 0; u <= 1.0; u += 0.01)
    {
        double s[p];
        // cout << "[loop] u=" << u_it /*<< " " <<s[0] << " " << s[1] << " " << s[2] */<< endl;
       BSplinePoint(u, U, n_knot, p, P, dim, s);
       cout << "[loop] u=" << u << " " <<s[0] << " " << s[1] << " " << s[2] << endl<< endl;
       file << s[0] << " " << s[1] << " " << s[2] << endl;

        // vec.push_back(Eigen::Vector3f(s[0], s[1], s[2]));
    }

    double s[p];
    BSplinePoint(1.0, U, n_knot, p, P, dim, s);
    cout << "u=1"<<s[0] << " " << s[1] << " " << s[2] << endl<< endl;
    file << s[0] << " " << s[1] << " " << s[2] << endl;


    // cout << "uit " << u_it << endl;

//    if (u_it == 1.0){

//        cout << "is 1" << endl<< endl;

//        u_it+=0.01;
//        double s[p];
//        BSplinePoint(u_it, U, n_knot, p, P, dim, s);

//        cout << "u=" << u_it << " " <<s[0] << " " << s[1] << " " << s[2] << endl<< endl;
//    }

    file.close();

    //    for (uint i = 0; i < vec.size(); i++)
    //        cout << vec.at(i)(0) << ", " << vec.at(i)(1) << ", " <<vec.at(i)(2) << endl;

    return 0;
}

//bool write(double){
//    std::ofstream file;
//        file.open("/home/manuel/hrg/spline_interpolation/output/spline.txt", std::ofstream::out | std::ofstream::app);
//        for (int i = 0; i < final_list_joint.size(); ++i)
//        {
//                q_tmp = final_list_joint.at(i);
//                p_tmp = final_list_cart.at(i);
//                if (ROB_DOF_SIZE == 7){
//                    file << q_tmp(0) << " " << q_tmp(1) << " " << q_tmp(2) << " " << q_tmp(3) << " " << q_tmp(4) << " " << q_tmp(5) << " " << q_tmp(6) << " ";
//                }else{
//                    file << q_tmp(0) << " " << q_tmp(1) << " " << q_tmp(2) << " " << q_tmp(3) << " " << q_tmp(4) << " " << q_tmp(5) << " ";
//                }
//                file << p_tmp(0) << " " << p_tmp(1) << " " << p_tmp(2) << " " << p_tmp(3) << " " << p_tmp(4) << " " << p_tmp(5) << std::endl;
//        }

//        file.close();
//}

void solve(double *a, double *b, double *c, double *d, int n)
{

    // n is the number of unknowns
    double c_tmp[n];
    for (uint i = 0; i < n; i++)
        c_tmp[i] = c[i];

    n--; // since we start from x0 (not x1)

    // make a copy of c, so c gets not changed by the algorithm
    c_tmp[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++)
    {
        c_tmp[i] /= b[i] - a[i] * c_tmp[i - 1];
        d[i] = (d[i] - a[i] * d[i - 1]) / (b[i] - a[i] * c_tmp[i - 1]);
    }

    d[n] = (d[n] - a[n] * d[n - 1]) / (b[n] - a[n] * c_tmp[n - 1]);

    for (int i = n; i-- > 0;)
    {
        d[i] -= c_tmp[i] * d[i + 1];
    }
}

void BSplinePoint(double u, double U[], int n_knot, int p,
                  std::vector<Eigen::Vector3f> P, int d, double s[])
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
            s[k] = s[k] + P.at(i - p + j)(k) * B[j];
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

    if (u == U[high]) // TODO check if correct (p. 470 in Trajectory_Planning_for_Automatic_Machines_and_Robots.pdf)
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
