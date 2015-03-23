#include <armadillo>
#include <iostream>
using namespace std;
using namespace arma;


struct MomentumVector
{
  int id;
  vector<double> p;
};

int main() {
	mat test = randu<mat>(3,3);
	cout << test << endl;

	mat A = randu<mat>(5,5);
	vec b = randu<vec>(5);
	mat B = randu<mat>(5,5);

	vec x = solve(A, b);
	cout << x << endl;

	MomentumVector p1;
	p1.p = {1,2,3,4};

	cout << p1.p[0] << p1.p[1] << p1.p[2] << p1.p[3] << endl;

	return 0;
}