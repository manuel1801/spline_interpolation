#ifndef SPLINEINTERPOLATION_HPP
#define SPLINEINTERPOLATION_HPP

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

class SplineInterpolation
{

public:
    SplineInterpolation();
    ~SplineInterpolation();

    bool generate_cubic_b_spline(std::vector<Eigen::Vector3f> q);
    bool get_waypoint_at(double u, std::vector<Eigen::Vector3f> &);

    //private:
    void BasisFuns(int i, double u, int p, double U[], double B[]);
    int WhichSpan(double u, double U[], int n_knot, int p);
    void solveMatrix(double *a, double *b, double *c, double *d, int n);
    void BSplinePoint(double u, double U[], int n_knot, int p, std::vector<Eigen::Vector3f> P, int d, double s[]);

    int p, dim, n, n_knot;          // p:polynomal degree, dim:cart space dimension
    Eigen::Vector3f t_0, t_n;       // can be removed
    std::vector<Eigen::Vector3f> P; // controll points
    // extern double u_knot[n_knot];
    std::vector<double> u_knot;
};

#endif // SPLINEINTERPOLATION_HPP
