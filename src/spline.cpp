#include <iostream>
#include <math.h>
#include <vector>
#include <kdl/velocityprofile_trap.hpp>
#include <kdl/rotational_interpolation_sa.hpp>
#include <kdl/path_roundedcomposite.hpp>
#include <kdl/trajectory_segment.hpp>

#define ROBOT_DOF_SIZE 7
#define MAX_P 100

void BasisFuns(int i, double u, int p, double U[], double B[]);
int WhichSpan(double u, double U[], int n_knot, int p);
void BSplinePoint( double u, double U[], int n_knot, int p, double P[],int d, double s[]);


using namespace std;
int main()
{

    double U[16] = {0, 0, 0, 0, 0.11, 0.23, 0.35, 0.48, 0.60, 0.68, 0.77, 0.84, 1, 1, 1, 1};
    double P[16*3] = {83, 34, -168, 146, -182, -45, 207, 31, -89, -29, 32, 71,
                      -5, 128, -41, 177, 88, 21, 14, -172, 30, 90,-54, -32,
                      119, 120,
                      88, 252, 245, 68, 98, 83, 121, 218, 188, 192};

    int n_knot = sizeof(U)/sizeof(U[0]) - 1;
    cout << "n_knot " << n_knot << endl;
    int p = 3, d = 3;


    for (double u = 0; u <=1; u+=0.01)

    {

        double s[p];

        BSplinePoint(u, U, n_knot, p, P, d, s);


        //for (uint i = 0; i < p; i++)
             cout << s[2] <<  ", " ;



        //int i = WhichSpan(u, U, n_knot, p);

       // double B[p];

       // BasisFuns(i, u, p, U, B);

       // for (uint i = 0; i < p; i++)
         //   cout << B[2] <<  ", " ;

        cout << endl;

    }



    return 0;
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
- value of the B-spline at u
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




