// C++ minimization.

#include <iostream>
#include <vector>
#include <armadillo>
#include <math.h>
#include <string>
#include <sstream>
#include <algorithm>
#include <iterator>
using namespace std;
using namespace arma;



// Implementation of Nelder-Mead Simplex method:
double * alloc_vector(int cols)
{
	return (double *) malloc(sizeof(double) * cols);
}
void free_vector(double * vector , int cols)
{
	free(vector);
}
double ** alloc_matrix(int rows, int cols)
{
	int	i;
	double ** matrix = (double **) malloc(sizeof(double *) * rows);
	for (i = 0; i < rows; i++)
		matrix[i] = alloc_vector(cols);
	return matrix;
}
void free_matrix(double ** matrix, int rows, int cols)
{
	int	i;
	for (i =0; i < rows; i++)
		free_vector(matrix[i], cols);
	free(matrix);
}
double ** make_simplex(double * point, int dim)
{
	int i, j;
	double ** simplex = alloc_matrix(dim + 1, dim);
	for (i = 0; i < dim + 1; i++)
		for (j = 0; j < dim; j++)
			simplex[i][j] = point[j];
	for (i = 0; i < dim; i++)
		simplex[i][i] += 1.0;
	return simplex;
}
void evaluate_simplex(double ** simplex, int dim,double * fx,  double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &),
	int Nevents, int j, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists)
{
	int i;
	for (i = 0; i < dim + 1; i++)
	fx[i] = (*func)(simplex[i], Nevents, j, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists);
}
void simplex_extremes(double *fx, int dim, int & ihi, int & ilo,int & inhi)
{
	int i;
	if (fx[0] > fx[1])
	{ ihi = 0; ilo = inhi = 1; }
	else
	{ ihi = 1; ilo = inhi = 0; }
	for (i = 2; i < dim + 1; i++)
		if (fx[i] <= fx[ilo])
			ilo = i;
		else if (fx[i] > fx[ihi])
			{ inhi = ihi; ihi = i; }
		else if (fx[i] > fx[inhi])
			inhi = i;
}
void simplex_bearings(double ** simplex, int dim,double * midpoint, double * line, int ihi)
{
	int i, j;
	for (j = 0; j < dim; j++)
		midpoint[j] = 0.0;
	for (i = 0; i < dim + 1; i++)
		if (i != ihi)
			for (j = 0; j < dim; j++)
				midpoint[j] += simplex[i][j];
	
	for (j = 0; j < dim; j++)
	{
		midpoint[j] /= dim;
		line[j] = simplex[ihi][j] - midpoint[j];
	}
}
int update_simplex(double * point, int dim, double & fmax,double * midpoint, double * line, double scale, double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &),
	int Nevents, int j, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists)
{
	int i, update =	0; 
	double * next = alloc_vector(dim), fx;
	for (i = 0; i < dim; i++)
		next[i] = midpoint[i] + scale * line[i];
	fx = (*func)(next, Nevents, j, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists);
	if (fx < fmax)
	{
		for (i = 0; i < dim; i++)	
			point[i] = next[i];
		fmax = fx;
		update = 1;
	}
	free_vector(next, dim);
	return update;
}

void contract_simplex(double ** simplex, int dim, double * fx, int ilo, double (*func)(double *, int))
{
	int i, j;
	for (i = 0; i < dim + 1; i++)
		if (i != ilo)
		{
			for (j = 0; j < dim; j++)
				simplex[i][j] = (simplex[ilo][j]+simplex[i][j])*0.5;
			fx[i] = (*func)(simplex[i], dim);
		}
}


#define ZEPS 1e-10
int check_tol(double fmax, double fmin, double ftol)
{
double delta = fabs(fmax - fmin);
double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
return (delta < (accuracy + ZEPS));
}


double amoeba(double *point, double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &), 
	double tol, 
	int Nevents, int j, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists)
{
	// Usage: Point is an allocated dim-dimensional array of doubles
	// to be filled with coordinates of the best-fit point,
	// func is the function to minimize. func takes the coordinates 
	// and the dim as input.
	int dim = 4;
	int ihi, ilo, inhi, j;
	double fmin;
	double * fx = alloc_vector(dim + 1);
	double * midpoint = alloc_vector(dim);
	double * line = alloc_vector(dim);
	double ** simplex = make_simplex(point, dim);
	evaluate_simplex(simplex, dim, fx, func);

	while (true)
	{
		simplex_extremes(fx, dim, ihi, ilo, inhi);
		simplex_bearings(simplex, dim, midpoint, line, ihi);
		if (check_tol(fx[ihi], fx[ilo], tol)) break;
		update_simplex(simplex[ihi], dim, fx[ihi],
		midpoint, line, -1.0, func);
		if (fx[ihi] < fx[ilo])
			update_simplex(simplex[ihi], dim, fx[ihi], midpoint, line, -2.0, func);
		else if (fx[ihi] >= fx[inhi])
			if (!update_simplex(simplex[ihi], dim, fx[ihi], midpoint, line, 0.5, func))
				contract_simplex(simplex, dim, fx, ilo, func);
	}

	for (j = 0; j < dim; j++)
		point[j] = simplex[ilo][j];
	fmin = fx[ilo];
	free_vector(fx, dim);
	free_vector(midpoint, dim);
	free_vector(line, dim);
	free_matrix(simplex, dim + 1, dim);
	return fmin;
}







double testfun(double *x, int dim)
{
	return x[1]*x[1] + x[0]*x[0];
}


struct MomentumVector
{
  int id;
  arma::vec p;
};

double minkowskidot(vec a, vec b)
{
	return a[3]*b[3]-a[0]*b[0]-a[1]*b[1]-a[2]*b[2];
}

double xisquared(double *Masses, int Nevents, int j, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists)
{
	vec M;
	M << Masses[0]*Masses[0] << Masses[1]*Masses[1] << Masses[2]*Masses[2] << Masses[3]*Masses[3] 
		<< Masses[0]*Masses[0] << Masses[1]*Masses[1] << Masses[2]*Masses[2] << Masses[3]*Masses[3];
	M = M/pow(Mnorm, 2);

	double xisquared = 0;

	if combinatorics
	{
		for (int iEvent = j, iEvent < j+Nevents, iEvent++)
		{
			int Ncombinations;
			if all_leptons_equal_list[i]
				Ncombinations = 8;
			else 
				Ncombinations = 16;
	
			for (iCombinations = 0, iCombinations < Ncombinations, iCombinations++)
			{
				vector<double> xisquared_current_list;
				double xisquared_current;
				vec P;
	
				P = D_lists[iCombinations][iEvent]*M + E_lists[iCombinations][iEvent];
	
				xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - M[3], 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - M[3], 2);
				xisquared_current_list.push_back(xisquared_current);
			}
	
			xisquared = xisquared + *min_element(xisquared_current_list, xisquared_current_list+Ncombinations);
		}
	}
	else // NOT combinatorics
	{
		for (int iEvent = j, iEvent < j+Nevents, iEvent++)
		{
			double xisquared_current;
			xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - M[3], 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - M[3], 2);
			xisquared = xisquared + xisquared_current_list;	
		}
		
	}

	return xisquared;

}

void best_fit(int Nbins, int Nevents, double *Masses, bool combinatorics, )
{
	int N = Nbins*Nevents;

	// Define permutation matrices
	mat permute23;
	mat permute67;
	mat permute23and67;	
	mat B;

	permute23 		 <<	1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 1 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 1 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 1 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 1 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 1 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << endr;
	permute67 	 << 	1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 1 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 1 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 1 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 1 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 1 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << endr;
	permute23and67 	  <<1 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 1 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 1 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 1 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 1 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 1 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 0 << 1 << endr;

	B 		 <<		   -1 << 1 << 0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << -1 << 1 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << -1 << 1 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 <<  0 << 0 << 0 << 0 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << -1 << 1 << 0 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << -1 << 1 << 0 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << -1 << 1 << endr << 
						0 << 0 << 0 << 0 << 0 << 0 << 0 << 0 << endr;				


	// Declare vectors of matrices to be stored

	vector<mat> D11_list, D12_list, D13_list, D14_list, D21_list, D22_list, D23_list, D24_list, D31_list, D32_list, D33_list, D34_list, D41_list, D42_list, D43_list, D44_list;
	vector<vec> E11_list, E12_list, E13_list, E14_list, E21_list, E22_list, E23_list, E24_list, E31_list, E32_list, E33_list, E34_list, E41_list, E42_list, E43_list, E44_list; // D11 is case 1 unpermuted, D12 is case 1 with permutation of 2 and 3, etc.

	vector<vector<mat>> D_lists;
	vector<vector<vec>> E_lists;

	vector<bool> all_leptons_equal_list;


	string line;
	ifstream events ("../python/Pythia_cascade_events_20150120.dat");
	if (events.is_open())
	{
		for (int iEvent = 0; iEvent < N; iEvent++) 
		{
			vector<string> particle;
			MomentumVector p1, p2, p3, p4, p5, p6, p7, p8;

			for (int iParticle = 0; iParticle < 9; iParticle++)
			{
				getline(events,line);
				// cout << iParticle << " " << line << endl;
				istringstream iss(line);
				particle = {istream_iterator<string>{iss}, istream_iterator<string>{}};

				if (iParticle == 1)
				{	
					p1.id = stoi(particle[0]);
					p1.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 2)
				{
					p2.id = stoi(particle[0]);
					// cout << abs(p2.id) << endl;
					p2.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 3)
				{
					p3.id = stoi(particle[0]);
					// cout << abs(p3.id) << endl;
					p3.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 4)
				{
					p4.id = stoi(particle[0]);
					p4.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 5)
				{
					p5.id = stoi(particle[0]);
					p5.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 6)
				{
					p6.id = stoi(particle[0]);
					// cout << abs(p6.id) << endl;
					p6.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 7)
				{
					p7.id = stoi(particle[0]);
					// cout << abs(p7.id) << endl;
					p7.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
				if (iParticle == 8)
				{
					p8.id = stoi(particle[0]);
					p8.p << stod(particle[1]) << stod(particle[2]) << stod(particle[3]) << stod(particle[4]);
				}
			}

			all_leptons_equal_list.push_back((abs(p2.id)==abs(p6.id)));



			double m1squared = minkowskidot(p1.p, p1.p);
			double m2squared = minkowskidot(p2.p, p2.p);
			double m3squared = minkowskidot(p3.p, p3.p);
			double m4squared = minkowskidot(p4.p, p4.p);
			double m5squared = minkowskidot(p5.p, p5.p);
			double m6squared = minkowskidot(p6.p, p6.p);
			double m7squared = minkowskidot(p7.p, p7.p);
			double m8squared = minkowskidot(p8.p, p8.p);
			double pxmiss 	 = p4.p[0] + p8.p[0];
			double pymiss 	 = p4.p[1] + p8.p[1];

			mat A1, A2, A3, A4;
			vec C1, C2, C3, C4;

			A1 <<	p1.p[1] << p1.p[2] << p1.p[3] << -p1.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
					p2.p[1] << p2.p[2] << p2.p[3] << -p2.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
					p3.p[1] << p3.p[2] << p3.p[3] << -p3.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
					0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
					0 	<< 0 	<< 0 	<< 0 	<< p5.p[1] << p5.p[2] << p5.p[3] << -p5.p[0] << endr <<
					0 	<< 0 	<< 0 	<< 0 	<< p6.p[1] << p6.p[2] << p6.p[3] << -p6.p[0] << endr <<
					0 	<< 0 	<< 0 	<< 0 	<< p7.p[1] << p7.p[2] << p7.p[3] << -p7.p[0] << endr <<
					0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
			A1 = A1*2/pow(Mnorm, 2)
			C1 	<< 	2*minkowskidot(p1.p, p2.p) + 2*minkowskidot(p1.p, p3.p) + m1squared
				<<	2*minkowskidot(p2.p, p3.p) + m2squared
				<<	m3squared
				<<	pow(pxmiss, 2)
				<<	2*minkowskidot(p5.p, p6.p) + 2*minkowskidot(p5.p, p7.p) + m5squared
				<<	2*minkowskidot(p6.p, p7.p) + m6squared
				<< 	m7squared
				<< 	pow(pymiss, 2);

			D11_list.push_back(inv(A1)*B);
			D12_list.push_back(inv(permute23*A1)*B);
			D12_list.push_back(inv(permute67*A1)*B);
			D12_list.push_back(inv(permute23and67*A1)*B);
			E11_list.push_back(inv(A1)*C1);
			E12_list.push_back(inv(permute23*A1)*C1);
			E12_list.push_back(inv(permute67*A1)*C1);
			E12_list.push_back(inv(permute23and67*A1)*C1);



			if combinatorics
			{

				A2 <<	p5.p[1] << p5.p[2] << p5.p[3] << -p5.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
						p2.p[1] << p2.p[2] << p2.p[3] << -p2.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
						p3.p[1] << p3.p[2] << p3.p[3] << -p3.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
						0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
						0 	<< 0 	<< 0 	<< 0 	<< p1.p[1] << p1.p[2] << p1.p[3] << -p1.p[0] << endr <<
						0 	<< 0 	<< 0 	<< 0 	<< p6.p[1] << p6.p[2] << p6.p[3] << -p6.p[0] << endr <<
						0 	<< 0 	<< 0 	<< 0 	<< p7.p[1] << p7.p[2] << p7.p[3] << -p7.p[0] << endr <<
						0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
				A2 = A2*2/pow(Mnorm, 2)
				C2 	<< 	2*minkowskidot(p5.p, p2.p) + 2*minkowskidot(p5.p, p3.p) + m5squared
					<<	2*minkowskidot(p2.p, p3.p) + m2squared
					<<	m3squared
					<<	pow(pxmiss, 2)
					<<	2*minkowskidot(p1.p, p6.p) + 2*minkowskidot(p1.p, p7.p) + m1squared
					<<	2*minkowskidot(p6.p, p7.p) + m6squared
					<< 	m7squared
					<< 	pow(pymiss, 2);
	
				D21_list.push_back(inv(A2)*B);
				D22_list.push_back(inv(permute23*A2)*B);
				D22_list.push_back(inv(permute67*A2)*B);
				D22_list.push_back(inv(permute23and67*A2)*B);
				E21_list.push_back(inv(A2)*C2);
				E22_list.push_back(inv(permute23*A2)*C2);
				E22_list.push_back(inv(permute67*A2)*C2);
				E22_list.push_back(inv(permute23and67*A2)*C2);
	
				if all_leptons_equal_list[iEvent]
				{
					if (p3.id*p7.id < 0)
					{
						p67_was_flipped = true;
						MomentumVector p6copy = p6;
						p6 = p7;
						p7 = p6copy;
					}
	
					A3 <<	p1.p[1] << p1.p[2] << p1.p[3] << -p1.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
							p2.p[1] << p2.p[2] << p2.p[3] << -p2.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
							p7.p[1] << p7.p[2] << p7.p[3] << -p7.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
							0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p5.p[1] << p5.p[2] << p5.p[3] << -p5.p[0] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p6.p[1] << p6.p[2] << p6.p[3] << -p6.p[0] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p3.p[1] << p3.p[2] << p3.p[3] << -p3.p[0] << endr <<
							0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
					A3 = A3*2/pow(Mnorm, 2)
					C3 	<< 	2*minkowskidot(p1.p, p2.p) + 2*minkowskidot(p1.p, p7.p) + m1squared
						<<	2*minkowskidot(p2.p, p7.p) + m2squared
						<<	m7squared
						<<	pow(pxmiss, 2)
						<<	2*minkowskidot(p5.p, p6.p) + 2*minkowskidot(p5.p, p3.p) + m5squared
						<<	2*minkowskidot(p6.p, p3.p) + m6squared
						<< 	m3squared
						<< 	pow(pymiss, 2);
	
					D31_list.push_back(inv(A3)*B);
					D32_list.push_back(inv(permute23*A3)*B);
					D32_list.push_back(inv(permute67*A3)*B);
					D32_list.push_back(inv(permute23and67*A3)*B);
					E31_list.push_back(inv(A3)*C3);
					E32_list.push_back(inv(permute23*A3)*C3);
					E32_list.push_back(inv(permute67*A3)*C3);
					E32_list.push_back(inv(permute23and67*A3)*C3);	
	
					A4 <<	p5.p[1] << p5.p[2] << p5.p[3] << -p5.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
							p2.p[1] << p2.p[2] << p2.p[3] << -p2.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
							p7.p[1] << p7.p[2] << p7.p[3] << -p7.p[0] << 	0 <<	0 <<	0 <<	0 << endr <<
							0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p1.p[1] << p1.p[2] << p1.p[3] << -p1.p[0] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p6.p[1] << p6.p[2] << p6.p[3] << -p6.p[0] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p3.p[1] << p3.p[2] << p3.p[3] << -p3.p[0] << endr <<
							0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
					A4 = A4*2/pow(Mnorm, 2)
					C4 	<< 	2*minkowskidot(p5.p, p2.p) + 2*minkowskidot(p5.p, p7.p) + m5squared
						<<	2*minkowskidot(p2.p, p7.p) + m2squared
						<<	m7squared
						<<	pow(pxmiss, 2)
						<<	2*minkowskidot(p1.p, p6.p) + 2*minkowskidot(p1.p, p3.p) + m1squared
						<<	2*minkowskidot(p6.p, p3.p) + m6squared
						<< 	m3squared
						<< 	pow(pymiss, 2);
	
					D41_list.push_back(inv(A4)*B);
					D42_list.push_back(inv(permute23*A4)*B);
					D42_list.push_back(inv(permute67*A4)*B);
					D42_list.push_back(inv(permute23and67*A4)*B);
					E41_list.push_back(inv(A4)*C4);
					E42_list.push_back(inv(permute23*A4)*C4);
					E42_list.push_back(inv(permute67*A4)*C4);
					E42_list.push_back(inv(permute23and67*A4)*C4);
				}
				else // NOT same-flavour lepton chains, no need for 3 and 4.
				{
					D31_list.push_back(0);
					D32_list.push_back(0);
					D33_list.push_back(0);
					D34_list.push_back(0);
					D41_list.push_back(0);
					D42_list.push_back(0);
					D43_list.push_back(0);
					D44_list.push_back(0);
					E31_list.push_back(0);
					E32_list.push_back(0);
					E33_list.push_back(0);
					E34_list.push_back(0);
					E41_list.push_back(0);
					E42_list.push_back(0);
					E43_list.push_back(0);
					E44_list.push_back(0);
				}

			} // END IF combinatorics
			else
			{
				// Don't need to do anythinge else
			}



		// END FOR loop over events
		}
	    events.close();
	}
	else cout << "Unable to open file"; 

	// Store matrices and vectors compactly
	D_lists = {D11_list, D12_list, D13_list, D14_list, D21_list, D22_list, D23_list, D24_list, D31_list, D32_list, D33_list, D34_list, D41_list, D42_list, D43_list, D44_list}
	E_lists = {E11_list, E12_list, E13_list, E14_list, E21_list, E22_list, E23_list, E24_list, E31_list, E32_list, E33_list, E34_list, E41_list, E42_list, E43_list, E44_list}


	// Finished with making the D and E matrices. Now to minimize xisquared!


	vector<double> best_fit_value;
	vector<double[4]> best_fit_point; // figure out how to store this
	for (int iBin=0, iBin<Nbins, iBin++)
	{
		double best_fit_current;
		double Masses_current[4] = Masses_initial;
		best_fit_current = amoeba(Masses_current, xisquared, combinatorics=true);

		best_fit_value.push_back(best_fit_current);
		best_fit_point.push_back(Masses_current);
	}




}




int main()
{


	// double xvector [2] = {-30,2};
	// cout << amoeba(xvector, 2, testfun, 0.00001) << endl;
	// cout << xvector[0] << "\t" << xvector[1] << endl;

	int Nbins = 1;
	int Nevents = 1;
	bool combinatorics = false;
	double Masses_initial[4] = {568, 180, 144, 97};

	best_fit(Nbins, Nevents, *Masses_initial, combinatorics);



	return 1;
}