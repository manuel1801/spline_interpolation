#include "spline_interpolation.hpp"

using namespace std;

SplineInterpolation::SplineInterpolation()
{
}

SplineInterpolation::~SplineInterpolation()
{
}

bool SplineInterpolation::generate_cubic_b_spline(std::vector<Eigen::VectorXf> q, Eigen::VectorXf v_start, Eigen::VectorXf v_end)
{
    if (q.size() < 4){
        cout << "not enough points given!" << endl;
    }

    p = 3.0;   // polynomal degree
    dim = 3.0; // dimension of cartesian space
    n = q.size() - 1;
    n_knot = n + 6; // as index for knot vector


    // 1. create knot vector u
    for (uint i = 0; i < 4; i++)
        u_knots.push_back(0);

    double d;
    for (uint i = 1; i < q.size(); i++)
        d += (q.at(i) - q.at(i - 1)).norm();

    for (uint k = 1; k < n; k++)
        u_knots.push_back(u_knots[k + 2] + ((q.at(k) - q.at(k - 1)).norm()) / d);

    for (uint i = 0; i < 4; i++)
        u_knots.push_back(1.0);

//    for (uint k = 0; k < u_knots.size(); k++)
//        cout << u_knots[k] << ", ";

    // 2. create controll points p

    // add first two pints
    P.push_back(q.at(0));
    P.push_back(q.at(0) + ((u_knots[4] / 3.0) * v_start));

    //calculate p_2 .. p_n-3
    // vectors a, b, c, d for tridiagonal matrix PB = R
    double a[n - 1], b[n - 1], c[n - 1], d_x[n - 1], d_y[n - 1], d_z[n - 1];

    a[0] = 0;
    c[n - 2] = 0;

    // get Basis Function for u_1 ... u_n-1
    for (int i = 1; i <= n - 1; i++)
    {
        double B[p + 1];
        BasisFuns(WhichSpan(u_knots[i + 3], u_knots, n_knot, p), u_knots[i + 3], p, u_knots, B);

        // Fill in vectors for tridiagonal matrix of basis functions
        // a;
        if (i > 1)
            a[i - 1] = B[0];

        // b:
        b[i - 1] = B[1];

        // c:
        if (i <= n - 1)
            c[i - 1] = B[2];

        // d:
        if (i == 1)
        {
            d_x[i - 1] = (q.at(1) - B[0] * P.at(1))(0);
            d_y[i - 1] = (q.at(1) - B[0] * P.at(1))(1);
            d_z[i - 1] = (q.at(1) - B[0] * P.at(1))(2);
        }
        else if (i == n - 1)
        {
            d_x[i - 1] = (q.at(i) - B[2] * (q.at(n) - ((1.0 - u_knots[n + 2]) / 3.0) * v_end))(0); // [n+2] ?
            d_y[i - 1] = (q.at(i) - B[2] * (q.at(n) - ((1.0 - u_knots[n + 2]) / 3.0) * v_end))(1);
            d_z[i - 1] = (q.at(i) - B[2] * (q.at(n) - ((1.0 - u_knots[n + 2]) / 3.0) * v_end))(2);
        }
        else if (i < n - 1)
        {
            d_x[i - 1] = (q.at(i))(0);
            d_y[i - 1] = (q.at(i))(1);
            d_z[i - 1] = (q.at(i))(2);
        }
    }

    // solve PB = R for P with tridiagonal matrix algorthm, res in d
    solveMatrix(a, b, c, d_x, n - 1);
    solveMatrix(a, b, c, d_y, n - 1);
    solveMatrix(a, b, c, d_z, n - 1);

    for (uint i = 0; i < n - 1; i++)
    {
        P.push_back(Eigen::Vector3f(d_x[i], d_y[i], d_z[i]));
    }

    // add last two pints
    P.push_back(q.at(n) - ((1.0 - u_knots[n + 2]) / 3.0) * v_end); // TODO check if [n+2] correct
    P.push_back(q.at(n));

    // for (auto p : P)
    //     cout << p(0) << " " << p(1) << " " << p(2) << endl;
}



bool SplineInterpolation::get_waypoint_at(double u, Eigen::VectorXf &s)
{
    if (P.size() < 4){
        cout << "not enough controll points given" << endl;
    }

    double B[p + 1];
    int i, k, j;
    i = WhichSpan(u, u_knots, n_knot, p);
    BasisFuns(i, u, p, u_knots, B);
    for (k = 0; k < dim; k++) /* For each components of the B-spline*/
    {
        s(k) = 0;
        for (j = 0; j <= p; j++)
        {
            // s[k] = s[k] + P[k * (n_knot - p) + i - p + j] * B[j];
            s(k)  = s(k)  + P.at(i - p + j)(k) * B[j];
        }
    }
    return 0;
}



void SplineInterpolation::solveMatrix(double *a, double *b, double *c, double *d, int n)
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



void SplineInterpolation::BasisFuns(int i, double u, int p, std::vector<double> U, double B[])
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
    double DR[p + 1], DL[p + 1];
    B[0] = 1;
    for (j = 1; j <= p; j++)
    {
//        DL[j] = u - U[i + 1 - j];
//        DR[j] = U[i + j] - u;

        DL[j] = u - U.at(i + 1 - j);
        DR[j] = U.at(i + j) - u;


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

int SplineInterpolation::WhichSpan(double u, std::vector<double> U, int n_knot, int p)
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

    if (u == U[high])
        u -= 0.0000000000001; // check if correct (p. 470 in Trajectory_Planning_for_Automatic_Machines_and_Robots.pdf)

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

int main()
{

    SplineInterpolation splInterp;

    std::vector<Eigen::VectorXf> q;
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

    Eigen::VectorXf t_0 = Eigen::Vector3f(-1236, 538, 42);
    Eigen::VectorXf t_n = Eigen::Vector3f(732, 1130, 63);



    splInterp.generate_cubic_b_spline(q, t_0, t_n);

    Eigen::VectorXf p_out;
    std::ofstream file;
    file.open("/home/manuel/hrg/spline_interpolation/output/spline.txt");
    for (double k = 0; k <= 1.0; k += 0.01)
    {
        p_out.setZero(3);
        splInterp.get_waypoint_at(k, p_out);
        cout << "k. " << k << endl;
        file <<p_out(0) << " " << p_out(1) << " " << p_out(2) << endl;

     }
    p_out.setZero(3);
    splInterp.get_waypoint_at(1.0, p_out);
    file <<p_out(0) << " " << p_out(1) << " " << p_out(2) << endl;



    // knot vector u
    // double u_knots[n_knot] = {0, 0, 0, 0, 0.11, 0.23, 0.35, 0.48, 0.60, 0.68, 0.77, 0.84, 1.0, 1.0, 1.0, 1.0};

    //    double P[16*3] = {83, 34, -168, 146, -182, -45, 207, 31, -89, -29, 32, 71,
    //                      -5, 128, -41, 177, 88, 21, 14, -172, 30, 90,
    //                 -54, -32, 119, 120,88, 252, 245, 68, 98, 83, 121, 218, 188, 192};

    //  controll points from example
    //    {83.0,  34.0,   -168.0,   146.0, -182.0, -45.0,  207.0, 31.0, -89.0, -29.0,  32.0, 71.0,
    //    -54.0, -32.0,  -5.0,      128.0, -41.0,  177.0,  88.0,  21.0,  14.0, -172.0, 30.0,  90.0,
    //    ,119.0, 120.0, 88.0,      252.0, 245.0,  68.0,   98.0,  83.0,  121.0, 218.0, 188.0, 192.0};

    //    std::vector<Eigen::Vector3f> _P;
    //    _P.push_back(Eigen::Vector3f(83.0, -54.0, 119.0));
    //    _P.push_back(Eigen::Vector3f(34.0, -32.0, 120.0));
    //    _P.push_back(Eigen::Vector3f(-168.0, -5.0, 88.0));
    //    _P.push_back(Eigen::Vector3f(146.0, 128.0, 252.0));
    //    _P.push_back(Eigen::Vector3f(-182.0, - 41.0, 245.0));
    //    _P.push_back(Eigen::Vector3f(-45.0, 177.0, 68.0));
    //    _P.push_back(Eigen::Vector3f(207.0, 88.0, 98.0));
    //    _P.push_back(Eigen::Vector3f(31.0, 21.0, 83.0));
    //    _P.push_back(Eigen::Vector3f(-89.0, 14.0, 121.0));
    //    _P.push_back(Eigen::Vector3f(-29.0, -172.0, 218.0));
    //    _P.push_back(Eigen::Vector3f(32.0, 30.0, 188.0));
    //    _P.push_back(Eigen::Vector3f(71.0, 90.0, 192.0));

    //waypoints q
    //    double q[10*d] =    {83, -64 ,42 ,-98, -13 ,140, 43, -65 ,-45, 71,
    //                      -54 ,10 ,79 ,23 ,125 ,81 ,32 ,-17 ,-89, 90,
    //                     119 ,124 ,226 ,222 ,102 ,92 ,92 ,134 ,182 ,192};

    // derivatives at the endpoints

    // std::vector<Eigen::Vector3f> vec;
//    std::ofstream file;
//    file.open("/home/manuel/hrg/spline_interpolation/output/spline.txt");
//    for (double k = 0; k <= 1.0; k += 0.01)
//    {

//        splInterp.get_waypoint_at() double s[p];
//        splInterp.BSplinePoint(k, u_knots, n_knot, p, P, dim, s);
//        cout << "Spline Points at u: " << k << " x=" << s[0] << " y=" << s[1] << " Z=" << s[2] << endl;
//        file << s[0] << " " << s[1] << " " << s[2] << endl;
//        // vec.push_back(Eigen::Vector3f(s[0], s[1], s[2]));
//    }

//    double s[p];
//    splInterp.BSplinePoint(1.0, u_knots, n_knot, p, P, dim, s);
//    cout << "Spline Points at u: " << 1.0 << " x=" << s[0] << " y=" << s[1] << " Z=" << s[2] << endl
//         << endl;
//    file << s[0] << " " << s[1] << " " << s[2] << endl;
//    file.close();
//    return 0;
}






//void SplineInterpolation::BSplinePoint(double u, double U[], int n_knot, int p,
//                                       std::vector<Eigen::Vector3f> P, int d, double s[])
//
//Inputs:u
//- value of the independent variable
//U[]
//- Knot vector
//n_knot - length of U[] -1
//p
//- degree of the spline
//P[]
//- Control point vector
//d
//- dimensions of a control point (2 in 2D, 3 in 3D, etc.)
//Output:s[]
//- value of the B-spline at uthe
//
//{
//    double B[p + 1];
//    int i, k, j;
//    i = WhichSpan(u, U, n_knot, p);
//    BasisFuns(i, u, p, U, B);
//    for (k = 0; k < d; k++) /* For each components of the B-spline*/
//    {
//        s[k] = 0;
//        for (j = 0; j <= p; j++)
//        {
//            // s[k] = s[k] + P[k * (n_knot - p) + i - p + j] * B[j];
//            s[k] = s[k] + P.at(i - p + j)(k) * B[j];
//        }
//    }
//}

//void SplineInterpolation::BSplinePoint(double u, double U[], int n_knot, int p,
//                                       std::vector<Eigen::Vector3f> P, int d, double s[])
//
//Inputs:u
//- value of the independent variable
//U[]
//- Knot vector
//n_knot - length of U[] -1
//p
//- degree of the spline
//P[]
//- Control point vector
//d
//- dimensions of a control point (2 in 2D, 3 in 3D, etc.)
//Output:s[]
//- value of the B-spline at uthe
//
//{
//    double B[p + 1];
//    int i, k, j;
//    i = WhichSpan(u, U, n_knot, p);
//    BasisFuns(i, u, p, U, B);
//    for (k = 0; k < d; k++) /* For each components of the B-spline*/
//    {
//        s[k] = 0;
//        for (j = 0; j <= p; j++)
//        {
//            // s[k] = s[k] + P[k * (n_knot - p) + i - p + j] * B[j];
//            s[k] = s[k] + P.at(i - p + j)(k) * B[j];
//        }
//    }
//}



//bool SplineInterpolation::get_waypoint_at(double u, std::vector<Eigen::Vector3f> &p)
//{

//    cout << u_knots[4] << endl;

//    // here make a check if control points were created
//    // 3. evaluate spline for given u (using BSplinePoint())
//    //TODO: translate into t (time)

//    // for (auto p : P)
//    //     cout << p(0) << " " << p(1) << " " << p(2) << endl;

//    // std::vector<Eigen::Vector3f> vec;
//    std::ofstream file;
//    file.open("/home/manuel/hrg/spline_interpolation/output/spline.txt");
//    for (double k = 0; k <= 1.0; k += 0.01)
//    {
//        double s[p];
//        BSplinePoint(k, u_knots, n_knot, p, P, dim, s);
//        cout << "Spline Points at u: " << k << " x=" << s[0] << " y=" << s[1] << " Z=" << s[2] << endl;
//        file << s[0] << " " << s[1] << " " << s[2] << endl;
//        // vec.push_back(Eigen::Vector3f(s[0], s[1], s[2]));
//    }

//    double s[p];
//    BSplinePoint(1.0, u_knots, n_knot, p, P, dim, s);
//    cout << "Spline Points at u: " << 1.0 << " x=" << s[0] << " y=" << s[1] << " Z=" << s[2] << endl
//         << endl;
//    file << s[0] << " " << s[1] << " " << s[2] << endl;
//    file.close();
//    return 0;
//}

