#include <iostream>
#include <math.h>
#include <vector>
#include <kdl/velocityprofile_trap.hpp>
#include <kdl/rotational_interpolation_sa.hpp>
#include <kdl/path_roundedcomposite.hpp>
#include <kdl/trajectory_segment.hpp>
#include <Eigen/Dense>

#define ROBOT_DOF_SIZE 7
#define MAX_P 4

void BasisFuns(int i, double u, int p, double U[], double B[]);
void BasisFunEigen(int i, double u, int p, double U[], Eigen::VectorXf &stdB);
int WhichSpan(double u, double U[], int n_knot, int p);
void BSplinePoint( double u, double U[], int n_knot, int p, double P[],int d, double s[]);
void solve_eigen(Eigen::VectorXf a, Eigen::VectorXf b, Eigen::VectorXf c, Eigen::VectorXf d, int n);
void solve(double* a, double* b, double* c, double* d, int n);

using namespace std;
int main()
{

    int p = 3, dim = 3;
    double U[16] = {0, 0, 0, 0, 0.11, 0.23, 0.35, 0.48, 0.60, 0.68, 0.77, 0.84, 1.0, 1.0, 1.0, 1.0};
    int n_knot = sizeof(U)/sizeof(U[0]) - 1;
    cout << "n_knot " << n_knot << endl;


//    double P[16*3] = {83, 34, -168, 146, -182, -45, 207, 31, -89, -29, 32, 71,
//                      -5, 128, -41, 177, 88, 21, 14, -172, 30, 90,-54, -32,
//                      119, 120,
//                      88, 252, 245, 68, 98, 83, 121, 218, 188, 192};

    //waypoints q
//    double q[10*d] =    {83, -64 ,42 ,-98, -13 ,140, 43, -65 ,-45, 71,
//                      -54 ,10 ,79 ,23 ,125 ,81 ,32 ,-17 ,-89, 90,
//                     119 ,124 ,226 ,222 ,102 ,92 ,92 ,134 ,182 ,192};

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
    cout << "n " << n << endl;

    // derivatives at the endpoints
    Eigen::Vector3f t_0 = Eigen::Vector3f(-1236, 538, 42);
    Eigen::Vector3f t_n = Eigen::Vector3f(732, 1130, 63);

    // control points
    std::vector<Eigen::Vector3f> P;

    P.push_back(q.at(0));
    //Eigen::Vector3f p_1 = q.at(0) + ((U[4]/3.0)*t_0);
    P.push_back(q.at(0) + ((U[4]/3.0)*t_0));

    // oder mit u -> which span i

    double a[n-1], b[n-1], c[n-1], d_x[n-1], d_y[n-1], d_z[n-1];
    std::vector<Eigen::Vector3f> d;


    a[0] = 0;
    c[n-2] = 0;
    for (int i = 0; i < n-1; i++){ //1..n-1

        // Eigen::VectorXf B_a(4);
        double B[MAX_P];


        int intervall = WhichSpan(U[i+1], U, n_knot, p);
        BasisFuns(intervall, U[i+1], p, U, B); //a

        cout << "intervall " << intervall << ", B1: " << B[0]<< ", B2: " << B[1]<< ", B3: " << B[2]<< ", B4: " << B[3] << endl;


        // for a: j == i
        // for b: j == i+1
        // for c: j == i + 2

        if (i > 0)
            a[i] = B[1];


        b[i] = B[0];

        if (i < (n-1))
            c[i-1] = B[2];


        // std::cout << "B_3 " << B[3] << endl;


//        for(uint j = 0; j< (sizeof(B)/sizeof(B[0])); j++)
//        {





//            a.push_back();

//            // fill a, b, c, d

//        }
//            cout << B[j] << endl;

//        cout << endl;

           // d:

        if (i == 1){

            d_x[i-1] = (q.at(1) - B[1]* P.at(1))(0);
            d_y[i-1] = (q.at(1) - B[1]* P.at(1))(1);
            d_z[i-1] = (q.at(1) - B[1]* P.at(1))(2);

           // d.push_back(q.at(1) - B[0]* P.at(1));

        }
        else if (i == n-1){

            d_x[i-1] = (q.at(i) - B[1] * (q.at(n) - ((1.0-U[n+3])/3.0)*t_n))(0);

            d_y[i-1] = (q.at(i) - B[1] * (q.at(n) - ((1.0-U[n+3])/3.0)*t_n))(1);

            d_z[i-1] = (q.at(i) - B[1] * (q.at(n) - ((1.0-U[n+3])/3.0)*t_n))(2);

           // d.push_back(q.at(i) - B[1] * (q.at(n) - ((1.0-U[n+3])/3.0)*t_n)); // q_n-1 - B_n+1 * p_n+1
        }

        else if (i < n-1){

            d_x[i-1] = (q.at(i))(0);
            d_y[i-1] = (q.at(i))(1);
            d_z[i-1] = (q.at(i))(2);

          //  d.push_back(q.at(i));

        }






        // now take from all B's the B with indce i+1 for the case of a
        // somit scalar wert -> nur std::vector<double > verwenden
//        a_mat.push_back(B_a);
//        Eigen::VectorXf B_b(4);
//        BasisFunEigen(i, U[i-1], p, U, B_b); //b

//        b_mat.push_back(B_b);

//        Eigen::VectorXf B_c(4);

//        if (i < (n-1)){
//            BasisFunEigen(i+1, U[i-1], p, U, B_c); //c
//            c_mat.push_back(B_c);
//        }

//        Eigen::VectorXf R(4);

//        if (i-1 == 1){
//            //R[1] = q[1] - B_1(u_1)p_1

//        }
//        else if ((i-1) > 1 && (i-1) < (n+1)){

//            //R[i-1] <<  q[i-1]

//        }else if ((i-1) == (n+1)){
//            //R[n-1] = B_1(u_-1)pn_+1
//        }

//        d_mat.push_back(B_d);

    }


    cout << "a: "  <<endl;

    for (uint i = 0; i < n-1; i++)
        cout << a[i] << ", ";

    cout << endl;
    cout << "b: " <<endl;


    for (uint i = 0; i < n-1; i++)
        cout << b[i] << ", ";
    cout << endl;
    cout << "c: " <<endl;


    for (uint i = 0; i < n-1; i++)
        cout << c[i] << ", ";

    cout << endl;
 //   cout << "d: " <<endl;

//    for (uint i = 0; i < d.size(); i++)
//        cout << d.at(i).transpose() << endl;
//    cout << endl;


    cout << "d_x: " <<endl;
    for (uint i = 0; i < n-1; i++)
        cout << d_x[i] << ", ";

    cout << endl;

    cout << "d_y: " <<endl;
    for (uint i = 0; i < n-1; i++)
        cout << d_y[i] << ", ";

    cout << endl;

    cout << "d_z: " <<endl;
    for (uint i = 0; i < n-1; i++)
        cout << d_z[i] << ", ";


    // now calculate Pn [n=2..n] with tridiagonal alg


    solve(a,b,c,d_x,n-1);
//    solve(a,b,c,d_y,n-1);
//    solve(a,b,c,d_z,n-1);

            for (int i = 0; i < n-1; i++) {
                cout << "[TEST] " <<  d_x[i] << endl;
            }


    // testing the tridiag alg.
//    int  n_ = 4;
//    double a_[4] = {0, -1, -1, -1 };
//    double b_[4] = { 4,  4,  4,  4 };
//    double c_[4] = {-1, -1, -1 ,0};
//    double d_[4] = { 5,  5, 10, 23 };
//    // results    { 2,  3,  5, 7  }

//    solve(a_,b_,c_,d_,n_);
//        for (int i = 0; i < n_; i++) {
//            cout << "[TEST] " <<  d_[i] << endl;
//        }



    for (uint i = 0; i < n-1; i++){
        P.push_back(Eigen::Vector3f(d_x[i], d_y[i], d_z[i]));
    }

    // add only after P_2 ... P_n were added
    P.push_back(q.at(n) - ((1.0-U[n+3])/3.0)*t_n);
    P.push_back(q.at(n));
    //Eigen::Vector3f p_1 = q.at(0) + ((U[4]/3.0)*t_0);






    cout << "FINAL RESULT OF P" << endl;
    for (auto item : P)
        cout << item.transpose() << endl;


      cout << endl;



    // derivitives of q at the endpoints
//    double t_0[3] = {-1236, 538, 42};
//    double  t_9[3] = {732, 1130, 63};

    // generate controll points:
    // p0
//    idx = 0;
//    double P[12*3];
//    for (uint i = 0; i < d; i++){
//        P[i*(n+2+1)+idx] = q[i*(n+1)+idx];
//        cout << P[i*(n+2+1)+idx] << ", ";
//    }

//    cout << endl;

//    // p1:
//    idx = 1;
//    for (uint i = 0; i < d; i++){
//        P[i*(n+2+1)+idx] = q[i*(n+1)+0]+(U[4]/3.0)*t_0[i];
//        cout << P[i*(n+2+1)+idx] << ", ";
//    }
//    cout << endl;

//    // p n+2:
//    idx = n+2;
//    for (uint i = 0; i < d; i++){
//        P[i*(n+2+1)+(n+2)] = q[i*(n+1)+n];
//        cout << P[i*(n+2+1)+(n+2)] << ", ";
//    }
//    cout << endl;


//    // p n+1:
//    idx = n+1;
//    idx = n+2;
//    for (uint i = 0; i < d; i++){
//        P[i*(n+2+1)+(n+1)] = q[i*(n+1)+n] - ((1-U[n+3])/3.0)*t_9[i];
//        cout << P[i*(n+2+1)+(n+1)] << ", ";
//    }
//    cout << endl;

//    std::vector<Eigen::VectorXf> a_mat, b_mat, c_mat, d_mat;
//    Eigen::VectorXf B_a(4), B_b(4),B_c(4),R(4);

//    for (int i = 2; i <= n; i++){
//        Eigen::VectorXf B_a(4);
//        if (i < n){
//            BasisFunEigen(i, U[i], p, U, B_a); //a
//            a_mat.push_back(B_a);
//        }
//        Eigen::VectorXf B_b(4);
//        BasisFunEigen(i, U[i-1], p, U, B_b); //b

//        b_mat.push_back(B_b);

//        Eigen::VectorXf B_c(4);

//        if (i < (n-1)){
//            BasisFunEigen(i+1, U[i-1], p, U, B_c); //c
//            c_mat.push_back(B_c);
//        }

//        Eigen::VectorXf R(4);

//        if (i-1 == 1){
//            //R[1] = q[1] - B_1(u_1)p_1

//        }
//        else if ((i-1) > 1 && (i-1) < (n+1)){

//            //R[i-1] <<  q[i-1]

//        }else if ((i-1) == (n+1)){
//            //R[n-1] = B_1(u_-1)pn_+1
//        }

//        d_mat.push_back(R);


//    }


    // oder mit u -> which span i

//    for (int i = 1; i < n; i++){ //1..n-1
//        Eigen::VectorXf B_a(4);
//        int i = WhichSpan(U[i], U, n_knot, p);
//        BasisFunEigen(i, U[i], p, U, B_a); //a
//        // now take from all B's the B with indce i+1 for the case of a
//        // somit scalar wert -> nur std::vector<double > verwenden
//        a_mat.push_back(B_a);
//        Eigen::VectorXf B_b(4);
//        BasisFunEigen(i, U[i-1], p, U, B_b); //b

//        b_mat.push_back(B_b);

//        Eigen::VectorXf B_c(4);

//        if (i < (n-1)){
//            BasisFunEigen(i+1, U[i-1], p, U, B_c); //c
//            c_mat.push_back(B_c);
//        }

//        Eigen::VectorXf R(4);

//        if (i-1 == 1){
//            //R[1] = q[1] - B_1(u_1)p_1

//        }
//        else if ((i-1) > 1 && (i-1) < (n+1)){

//            //R[i-1] <<  q[i-1]

//        }else if ((i-1) == (n+1)){
//            //R[n-1] = B_1(u_-1)pn_+1
//        }

//        d_mat.push_back(B_d);


//    }


 //    solve_eigen(a_mat,b_mat,c_mat,d_mat, n); // n ?

//    cout << "a" << endl;


//    for (uint i = 0; i < a_mat.size(); i++)
//        cout << a_mat.at(i).transpose()<< endl;
//    cout << "b" << endl;


//    for (uint i = 0; i < b_mat.size(); i++)
//        cout << b_mat.at(i).transpose() << endl;
//    cout << "c" << endl;


//    for (uint i = 0; i < c_mat.size(); i++)
//        cout << c_mat.at(i).transpose()<< endl;

    // then the remaining with tridiagonal matrix
    // a: for i = 2..n   -> B_i(u_i-1)
    // b: for i = 2..n-1 -> B_i(u_i)
    // c: for i = 3..n   -> B_i(u_i-2)
    // d: for r

    // P = calcMatrix(a,b,c,d)

    // then
    // BSplinePoint(u, U, n_knot, p, P, d, s)

return 0;


//    int n_knot = sizeof(U)/sizeof(U[0]) - 1;


//    for (double u = 0; u <=1; u+=0.01)

//    {

//        double s[p];
//        BSplinePoint(u, U, n_knot, p, P, d, s);

//        //for (uint i = 0; i < p; i++)
//       //      cout << s[2] <<  ", " ;
//        //cout << endl;

//    }



    return 0;
}


void solve(double* a, double* b, double* c, double* d, int n) {

    // n is the number of unknowns
/*
    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */

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


void solve_eigen(Eigen::VectorXf a, Eigen::VectorXf b, Eigen::VectorXf c, Eigen::VectorXf d, int n) {
    /*
    // n is the number of unknowns

    |b0 c0 0 ||x0| |d0|
    |a1 b1 c1||x1|=|d1|
    |0  a2 b2||x2| |d2|

    1st iteration: b0x0 + c0x1 = d0 -> x0 + (c0/b0)x1 = d0/b0 ->

        x0 + g0x1 = r0               where g0 = c0/b0        , r0 = d0/b0

    2nd iteration:     | a1x0 + b1x1   + c1x2 = d1
        from 1st it.: -| a1x0 + a1g0x1        = a1r0
                    -----------------------------
                          (b1 - a1g0)x1 + c1x2 = d1 - a1r0

        x1 + g1x2 = r1               where g1=c1/(b1 - a1g0) , r1 = (d1 - a1r0)/(b1 - a1g0)

    3rd iteration:      | a2x1 + b2x2   = d2
        from 2nd it. : -| a2x1 + a2g1x2 = a2r2
                       -----------------------
                       (b2 - a2g1)x2 = d2 - a2r2
        x2 = r2                      where                     r2 = (d2 - a2r2)/(b2 - a2g1)
    Finally we have a triangular matrix:
    |1  g0 0 ||x0| |r0|
    |0  1  g1||x1|=|r1|
    |0  0  1 ||x2| |r2|

    Condition: ||bi|| > ||ai|| + ||ci||

    in this version the c matrix reused instead of g
    and             the d matrix reused instead of r and x matrices to report results
    Written by Keivan Moradi, 2014
    */
    n--; // since we start from x0 (not x1)

    for (uint j = 0; j < 4; j++)
    c[0] /= b[0];
    d[0] /= b[0];

    for (int i = 1; i < n; i++) {
        c[i] /= b[i] - a[i]*c[i-1];
        d[i] = (d[i] - a[i]*d[i-1]) / (b[i] - a[i]*c[i-1]);
    }

    d[n] = (d[n] - a[n]*d[n-1]) / (b[n] - a[n]*c[n-1]);

    for (int i = n; i-- > 0;) {
        d[i] -= c[i]*d[i+1];
    }
}



void BSplinePoint( double u, double U[], int n_knot, int p,
double P[],int d, double s[])
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
            s[k] = s[k] + P[k * (n_knot - p) + i - p + j] * B[j];
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

void BasisFunEigen(int i, double u, int p, double U[], Eigen::VectorXf &B)
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
    B(0) = 1;
    for (j = 1; j <= p; j++)
    {
        DL[j] = u - U[i + 1 - j];
        DR[j] = U[i + j] - u;
        acc = 0.0;
        for (r = 0; r <= j - 1; r++)
        {
            temp = B[r] / (DR[r + 1] + DL[j - r]);
            B(r) = acc + DR[r + 1] * temp;
            acc = DL[j - r] * temp;
        }
        B(j) = acc;
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




