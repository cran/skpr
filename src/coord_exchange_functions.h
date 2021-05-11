
#include <function>

using namespace std;

//Return xmin: fa > fb < fc
double parabolicInterpolation(double a, double b, double c, std::function<double(double)> f) {
  double fa = f(a);
  double fb = f(b);
  double fc = f(c);
  double numerator = (fa-fb)*(c-b)*(c-b) - (fc-fb)*(b-a)*(b-a);
  double denominator = (fa-fb)*(c-b) + (fc-fb)*(b-a);
  double xmin = b + 1.0 / 2.0* numerator / denominator;
}
