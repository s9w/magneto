#ifndef ICING_HELPERS_H
#define ICING_HELPERS_H

#include <boost/math/special_functions/erf.hpp>
#include <boost/math/distributions/normal.hpp>

std::vector<double> linspace(double x_min, double x_max, unsigned int N){
    std::vector<double> T_vec;
    double dT = (x_max - x_min)/(N -1);
    if(N==1)
        dT=0;
    for(int i=0; i< N; ++i)
        T_vec.push_back(x_min + i*dT);
    return T_vec;
}

std::vector<double> normalSpace(double x_min, double x_max, unsigned int N, double mu, double sigma){
    if(N==1)
        return std::vector<double> {x_min};
    boost::math::normal normal_dist(mu, sigma);
    double area_min = boost::math::cdf(normal_dist, x_min);
    double area_max = boost::math::cdf(normal_dist, x_max);
    double A = area_max - area_min;
    double A0 = A/(N-1.0);
    double area;
    std::vector<double> x;
    for(int i=0; i<N; ++i) {
        area = i*A0 + area_min;
        x.push_back(sqrt(2.0)*sigma* boost::math::erf_inv(area*2.0 - 1.0) + mu);
    }
    return x;
}

#include <iostream>
#include <vector>

#endif //ICING_HELPERS_H
