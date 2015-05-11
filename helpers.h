#ifndef ICING_HELPERS_H
#define ICING_HELPERS_H

std::vector<double> linspace(double x_min, double x_max, unsigned int N);
std::vector<double> normalSpace(double x_min, double x_max, unsigned int N, double mu, double sigma);
std::string to_string(double const & value);

#endif
