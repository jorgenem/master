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



// Method for couting vectors of doubles 
std::ostream& operator<<(std::ostream &stream, vector<double> vec) {
	stream << "[";
	for (auto val:vec)
	{
		stream << val << ", " ;
	}
	return stream << "]";
}
// Method for couting vectors of ints 
std::ostream& operator<<(std::ostream &stream, vector<int> vec) {
	stream << "[";
	for (auto val:vec)
	{
		stream << val << ", " ;
	}
	return stream << "]";
}
// Method for couting vectors of bool 
std::ostream& operator<<(std::ostream &stream, vector<bool> vec) {
	stream << "[";
	for (auto val:vec)
	{
		if (val)
			stream << 1 << ", " ;
		else
			stream << 0 << ", " ;
	}
	return stream << "]";
}
// Count the number of true in a vector<bool>
int number_of_true(vector<bool> a)
{
	int number = 0;
	for (auto ai:a)
		if (ai)
		{
			number += 1;
		}
	return number;

}


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
		// simplex[i][i] += 1.0;
		simplex[i][i] *= 1.1;
	return simplex;
}
void evaluate_simplex(double ** simplex, int dim,double * fx,  double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &, vector<bool> &),
	int Nevents, int jBin, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists, vector<bool> &correct_combinatorics)
{
	int i;
	for (i = 0; i < dim + 1; i++)
	{
		correct_combinatorics.clear();
		fx[i] = (*func)(simplex[i], Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);
	}
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
int update_simplex(double * point, int dim, double & fmax,double * midpoint, double * line, double scale, double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &, vector<bool> &),
	int Nevents, int jBin, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists, vector<bool> &correct_combinatorics)
{
	int i, update =	0; 
	double * next = alloc_vector(dim), fx;
	for (i = 0; i < dim; i++)
		next[i] = midpoint[i] + scale * line[i];
	correct_combinatorics.clear();
	fx = (*func)(next, Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);
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

void contract_simplex(double ** simplex, int dim, double * fx, int ilo, double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &, vector<bool> &),	int Nevents, int jBin, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists, vector<bool> &correct_combinatorics)
{
	int i, j;
	for (i = 0; i < dim + 1; i++)
		if (i != ilo)
		{
			for (j = 0; j < dim; j++)
				simplex[i][j] = (simplex[ilo][j]+simplex[i][j])*0.5;
			correct_combinatorics.clear();
			fx[i] = (*func)(simplex[i], Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);
		}
}


#define ZEPS 1e-30
int check_tol(double fmax, double fmin, double ftol)
{
double delta = fabs(fmax - fmin);
double accuracy = (fabs(fmax) + fabs(fmin)) * ftol;
// cout << delta << ", " << accuracy << ", " << ftol << endl;
return (delta < (accuracy + ZEPS));
}

bool amoeba(double *point, double &fmin, double (*func)(double *, int, int, double, bool, vector<bool> &, vector<vector<mat>> &, vector<vector<vec>> &, vector<bool> &), 
	double tol, int maxiter,
	int Nevents, int jBin, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists, vector<bool> &correct_combinatorics)
{
	// Usage: Point is an allocated dim-dimensional array of doubles
	// to be filled with coordinates of the best-fit point,
	// func is the function to minimize. 
	int dim = 3; // MODIFIED TO FIT MD
	int ihi, ilo, inhi, j;
	// double fmin;
	double * fx = alloc_vector(dim + 1);
	double * midpoint = alloc_vector(dim);
	double * line = alloc_vector(dim);
	double ** simplex = make_simplex(point, dim);
	evaluate_simplex(simplex, dim, fx, func, 
		Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);

	int iter = 0;
	while (iter < maxiter)
	{
		simplex_extremes(fx, dim, ihi, ilo, inhi);
		simplex_bearings(simplex, dim, midpoint, line, ihi);
		if (check_tol(fx[ihi], fx[ilo], tol)) { /*cout << "below tol = " << tol << endl;*/ break; }
		update_simplex(simplex[ihi], dim, fx[ihi],
		midpoint, line, -1.0, func, 
		Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);
		if (fx[ihi] < fx[ilo])
			update_simplex(simplex[ihi], dim, fx[ihi], midpoint, line, -2.0, func,
				Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);
		else if (fx[ihi] >= fx[inhi])
			if (!update_simplex(simplex[ihi], dim, fx[ihi], midpoint, line, 0.5, func, Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics))
				contract_simplex(simplex, dim, fx, ilo, func, Nevents, jBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics);
		iter += 1;
	}

	for (j = 0; j < dim; j++)
		point[j] = simplex[ilo][j];
	fmin = fx[ilo];
	free_vector(fx, dim);
	free_vector(midpoint, dim);
	free_vector(line, dim);
	free_matrix(simplex, dim + 1, dim);

	if (iter < maxiter)
	{
		return true;
	}
	else
		return false;
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

double xisquared(double *Masses, int Nevents, int j, double Mnorm, bool combinatorics, vector<bool> &all_leptons_equal_list, vector<vector<mat>> &D_lists, vector<vector<vec>> &E_lists, vector<bool> &correct_combinatorics)
{

	mat B;
	B 		 <<		    -1 <<  0 <<  0  << endr << 
						0  << -1 <<  0  << endr << 
						0  <<  0 << -1  << endr << 
						0  <<  0 <<  0  << endr << 
						-1 <<  0 <<  0  << endr << 
						0  << -1 <<  0  << endr << 
						0  <<  0 << -1  << endr << 
						0  <<  0 <<  0  << endr;	

	vec M;
	M << Masses[0]
	  << Masses[1]
	  << Masses[2];
	M = M/pow(Mnorm, 2);

	double xisquared = 0;

	// Avoid regions of unphysical mass combinations by adding a huge contribution in a continuous way
	double hugefactor = 100000.0;
	if (Masses[0] < 0) xisquared = xisquared + hugefactor*M[0]*M[0];
	if (Masses[1] < 0) xisquared = xisquared + hugefactor*M[1]*M[1];
	if (Masses[2] < 0) xisquared = xisquared + hugefactor*M[2]*M[2];

	// Calculate current estimate for LSP mass from dilepton invariant mass edge
	double mllinv = 80; // Calculated from true masses using formula
	double MLSPsq = M[2]*(M[1]/(mllinv*mllinv) - 1.0);

	if (combinatorics)
	{
		for (int iEvent = j*Nevents; iEvent < (j+1)*Nevents; iEvent++)
		{
			int Ncombinations;
			if (all_leptons_equal_list[iEvent])
				Ncombinations = 4;
			else 
				Ncombinations = 2;
	
			vector<double> xisquared_current_list;
			// xisquared_current_list.clear();
			int iSmallest = 0;
			for (int iCombinations = 0; iCombinations < Ncombinations; iCombinations++)
			{

				double xisquared_current;
				vec P;
	
				// cout << "D has n_cols = " << D_lists[iCombinations][iEvent].n_cols << endl;
				P = D_lists[iCombinations][iEvent]*M + E_lists[iCombinations][iEvent];
				xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);

				// Add contributions from the three closest wrong combos, where the leptons are flipped inside chains
				P = D_lists[1+4*iCombinations][iEvent]*M + E_lists[1+4*iCombinations][iEvent];
				xisquared_current = xisquared_current + pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
				P = D_lists[2+4*iCombinations][iEvent]*M + E_lists[2+4*iCombinations][iEvent];
				xisquared_current = xisquared_current + pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
				P = D_lists[3+4*iCombinations][iEvent]*M + E_lists[3+4*iCombinations][iEvent];
				xisquared_current = xisquared_current + pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);

				xisquared_current_list.push_back(xisquared_current);
				if (xisquared_current_list[iCombinations] < xisquared_current_list[iSmallest])
					iSmallest = iCombinations;

			}
			// cout << "hei" << endl;
			if (iSmallest == 0)
				correct_combinatorics.push_back(1);
			else
				correct_combinatorics.push_back(0);

			
			// cout << "iEvent = " << iEvent << ", xisquared_current_list = " << xisquared_current_list << ", iSmallest = " << iSmallest << endl;

			xisquared = xisquared + *min_element(xisquared_current_list.begin(), xisquared_current_list.end());
		}
	}
	else // NOT combinatorics
	{
		for (int iEvent = j*Nevents; iEvent < (j+1)*Nevents; iEvent++)
		{
			double xisquared_current;
			vec P;

			// Check detA cut condition (HACK: Using correct_combinatorics to store passcut true/false)
			if (!correct_combinatorics[iEvent])
			{
				continue;
			}

			P = D_lists[0][iEvent]*M + E_lists[0][iEvent];
			// cout << "j = " << j << ", D(3,2) = " << D_lists[0][iEvent](3,2) << endl;

			// cout << "M[3] = " << M[3] << endl;

			xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);

			// Add contributions from the three first wrong combos, where the leptons are flipped inside chains
			P = D_lists[1][iEvent]*M + E_lists[1][iEvent];
			xisquared_current = xisquared_current + pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
			P = D_lists[2][iEvent]*M + E_lists[2][iEvent];
			xisquared_current = xisquared_current + pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
			P = D_lists[3][iEvent]*M + E_lists[3][iEvent];
			xisquared_current = xisquared_current + pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);

			xisquared = xisquared + xisquared_current;	

			// REMOVE COMMENT IF YOU WANT TO SUM ALL EIGHT COMBINATIONS
			// P = D_lists[4][iEvent]*M + E_lists[4][iEvent];
			// xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
			// xisquared = xisquared + xisquared_current;	
			// P = D_lists[5][iEvent]*M + E_lists[5][iEvent];
			// xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
			// xisquared = xisquared + xisquared_current;	
			// P = D_lists[6][iEvent]*M + E_lists[6][iEvent];
			// xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
			// xisquared = xisquared + xisquared_current;	
			// P = D_lists[7][iEvent]*M + E_lists[7][iEvent];
			// xisquared_current = pow(P[3]*P[3] - P[0]*P[0] - P[1]*P[1] - P[2]*P[2] - MLSPsq, 2) + pow(P[7]*P[7] - P[4]*P[4] - P[5]*P[5] - P[6]*P[6] - MLSPsq, 2);
			// xisquared = xisquared + xisquared_current;	

		}
		
	}


	// cout << "xisquared evaluated to = " << xisquared << endl;
	// xisquared = pow(M[0]-569*569/(Mnorm*Mnorm),4)+pow(M[1],4)+pow(M[2],4)+pow(M[3],4);

	// A test to try to make the xisquared steeper for increased resolution:
	// xisquared = pow(xisquared, 10);

	return xisquared;

}

void best_fit(int Nbins, int Nevents, string eventfile, vector<double> masses_initial, double tol, int maxiter, bool combinatorics, double Mnorm, vector<double> &best_fit_value, vector<vector<double> > &best_fit_point, vector<double> &correct_combinatorics_fraction, vector<int> &bin_number)
{
	int N = Nbins*Nevents;
	// cout << "N = " << endl;

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

	B 		 <<		    -1 <<  0 <<  0  << endr << 
						0  << -1 <<  0  << endr << 
						0  <<  0 << -1  << endr << 
						0  <<  0 <<  0  << endr << 
						-1 <<  0 <<  0  << endr << 
						0  << -1 <<  0  << endr << 
						0  <<  0 << -1  << endr << 
						0  <<  0 <<  0  << endr;				
 

	// Declare vectors of matrices to be stored

	vector<mat> D11_list, D12_list, D13_list, D14_list, D21_list, D22_list, D23_list, D24_list, D31_list, D32_list, D33_list, D34_list, D41_list, D42_list, D43_list, D44_list;
	vector<vec> E11_list, E12_list, E13_list, E14_list, E21_list, E22_list, E23_list, E24_list, E31_list, E32_list, E33_list, E34_list, E41_list, E42_list, E43_list, E44_list; // D11 is case 1 unpermuted, D12 is case 1 with permutation of 2 and 3, etc.

	vector<vector<mat>> D_lists;
	vector<vector<vec>> E_lists;

	vector<bool> all_leptons_equal_list;

	// Declare vector of detA values and correct_combinatorics (which as a hack is used to send in passcut true/false)
	vector<vector<double> > subdetAlist;
	vector<bool> correct_combinatorics;
	double subdetAcut = 1.0; // Value of detA cut (cut events with abs(detA) LOWER than this!)
						  // (note: Cut scale varies depending on Mnorm)


	string line;
	ifstream events (eventfile);

	if (events.is_open())
	{
		int shift = 0;
		for (int iEvent = 0; iEvent < shift; iEvent++)
		{
			for (int iParticle = 0; iParticle < 9; iParticle++)
			{
				getline(events,line);
			}
		}
		for (int iEvent = 0; iEvent < N; iEvent++) // if shift!=0, iEvent will not match the actual chain events from pythia/herwig
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
			// cout << (abs(p2.id)==abs(p6.id)) << endl;
			// cout << all_leptons_equal_list[iEvent] << endl;
			// cout << p1.id << ", " << p2.id << ", " << p3.id << ", " << p4.id << ", " << p5.id << ", " << p6.id << ", " << p7.id << ", " << p8.id << endl; 

			// if (iEvent == 0)
				// cout << p1.p << endl;

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


			A1 <<	p1.p[0] << p1.p[1] << p1.p[2] << -p1.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
					p2.p[0] << p2.p[1] << p2.p[2] << -p2.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
					p3.p[0] << p3.p[1] << p3.p[2] << -p3.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
					0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
					0 	<< 0 	<< 0 	<< 0 	<< p5.p[0] << p5.p[1] << p5.p[2] << -p5.p[3] << endr <<
					0 	<< 0 	<< 0 	<< 0 	<< p6.p[0] << p6.p[1] << p6.p[2] << -p6.p[3] << endr <<
					0 	<< 0 	<< 0 	<< 0 	<< p7.p[0] << p7.p[1] << p7.p[2] << -p7.p[3] << endr <<
					0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
			A1 = A1*2/Mnorm;
			C1 	<< 	2*minkowskidot(p1.p, p2.p) + 2*minkowskidot(p1.p, p3.p) + m1squared
				<<	2*minkowskidot(p2.p, p3.p) + m2squared
				<<	m3squared
				<<	pow(pxmiss, 2)
				<<	2*minkowskidot(p5.p, p6.p) + 2*minkowskidot(p5.p, p7.p) + m5squared
				<<	2*minkowskidot(p6.p, p7.p) + m6squared
				<< 	m7squared
				<< 	pow(pymiss, 2);
			C1 = C1/pow(Mnorm, 2);

			// Store detA value
			subdetAlist.push_back(vector<double> {det(A1.submat(0,0,2,2)), det(A1.submat(4,4,6,6))});

			// Check subdetA cut condition
			if (subdetAcut > abs(subdetAlist[iEvent][0]) || subdetAcut > abs(subdetAlist[iEvent][1]))
			{
				// HACK: Use correct_combinatorics vector to store detAcutcounter true/false value.
				// TRUE = PASSED CUT
				correct_combinatorics.push_back(0);
			}
			else
			{
				correct_combinatorics.push_back(1);
			}


			// mat test = inv(A1)*A1;
			// cout << test(0,0) << test(0,1) << endl;

			// Check if event is uninvertible
			// if (abs(det(A1))<1e-5)
			// {
			// 	cout << "det(A1) = " << det(A1) << endl;
			// 	cout << "Event number " << iEvent << endl;
			// }

			D11_list.push_back(inv(A1)*B);
			D12_list.push_back(inv(permute23*A1)*B);
			D13_list.push_back(inv(permute67*A1)*B);
			D14_list.push_back(inv(permute23and67*A1)*B);
			E11_list.push_back(inv(A1)*C1);
			E12_list.push_back(inv(permute23*A1)*C1);
			E13_list.push_back(inv(permute67*A1)*C1);
			E14_list.push_back(inv(permute23and67*A1)*C1);



			// if (combinatorics)
			if (true) // Evaluate all combinations irrespective -- doesn't take long
			{

				A2 <<	p5.p[0] << p5.p[1] << p5.p[2] << -p5.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
						p2.p[0] << p2.p[1] << p2.p[2] << -p2.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
						p3.p[0] << p3.p[1] << p3.p[2] << -p3.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
						0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
						0 	<< 0 	<< 0 	<< 0 	<< p1.p[0] << p1.p[1] << p1.p[2] << -p1.p[3] << endr <<
						0 	<< 0 	<< 0 	<< 0 	<< p6.p[0] << p6.p[1] << p6.p[2] << -p6.p[3] << endr <<
						0 	<< 0 	<< 0 	<< 0 	<< p7.p[0] << p7.p[1] << p7.p[2] << -p7.p[3] << endr <<
						0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
				A2 = A2*2/Mnorm;
				C2 	<< 	2*minkowskidot(p5.p, p2.p) + 2*minkowskidot(p5.p, p3.p) + m5squared
					<<	2*minkowskidot(p2.p, p3.p) + m2squared
					<<	m3squared
					<<	pow(pxmiss, 2)
					<<	2*minkowskidot(p1.p, p6.p) + 2*minkowskidot(p1.p, p7.p) + m1squared
					<<	2*minkowskidot(p6.p, p7.p) + m6squared
					<< 	m7squared
					<< 	pow(pymiss, 2);
				C2 = C2/pow(Mnorm, 2);
	
				D21_list.push_back(inv(A2)*B);
				D22_list.push_back(inv(permute23*A2)*B);
				D23_list.push_back(inv(permute67*A2)*B);
				D24_list.push_back(inv(permute23and67*A2)*B);
				E21_list.push_back(inv(A2)*C2);
				E22_list.push_back(inv(permute23*A2)*C2);
				E23_list.push_back(inv(permute67*A2)*C2);
				E24_list.push_back(inv(permute23and67*A2)*C2);
	
				if (all_leptons_equal_list[iEvent])
				{
					bool p67_was_flipped = false;
					if (p3.id*p7.id < 0)
					{
						p67_was_flipped = true;
						MomentumVector p6copy = p6;
						p6 = p7;
						p7 = p6copy;
					}
	
					A3 <<	p1.p[0] << p1.p[1] << p1.p[2] << -p1.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
							p2.p[0] << p2.p[1] << p2.p[2] << -p2.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
							p7.p[0] << p7.p[1] << p7.p[2] << -p7.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
							0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p5.p[0] << p5.p[1] << p5.p[2] << -p5.p[3] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p6.p[0] << p6.p[1] << p6.p[2] << -p6.p[3] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p3.p[0] << p3.p[1] << p3.p[2] << -p3.p[3] << endr <<
							0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
					A3 = A3*2/Mnorm;
					C3 	<< 	2*minkowskidot(p1.p, p2.p) + 2*minkowskidot(p1.p, p7.p) + m1squared
						<<	2*minkowskidot(p2.p, p7.p) + m2squared
						<<	m7squared
						<<	pow(pxmiss, 2)
						<<	2*minkowskidot(p5.p, p6.p) + 2*minkowskidot(p5.p, p3.p) + m5squared
						<<	2*minkowskidot(p6.p, p3.p) + m6squared
						<< 	m3squared
						<< 	pow(pymiss, 2);
					C3 = C3/pow(Mnorm, 2);
	
					D31_list.push_back(inv(A3)*B);
					D32_list.push_back(inv(permute23*A3)*B);
					D33_list.push_back(inv(permute67*A3)*B);
					D34_list.push_back(inv(permute23and67*A3)*B);
					E31_list.push_back(inv(A3)*C3);
					E32_list.push_back(inv(permute23*A3)*C3);
					E33_list.push_back(inv(permute67*A3)*C3);
					E34_list.push_back(inv(permute23and67*A3)*C3);	
	
					A4 <<	p5.p[0] << p5.p[1] << p5.p[2] << -p5.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
							p2.p[0] << p2.p[1] << p2.p[2] << -p2.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
							p7.p[0] << p7.p[1] << p7.p[2] << -p7.p[3] << 	0 <<	0 <<	0 <<	0 << endr <<
							0.5*pxmiss		<< 0 	<< 0 	<< 0 	<< 0.5*pxmiss 	<< 0 	<< 0 	<< 0 << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p1.p[0] << p1.p[1] << p1.p[2] << -p1.p[3] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p6.p[0] << p6.p[1] << p6.p[2] << -p6.p[3] << endr <<
							0 	<< 0 	<< 0 	<< 0 	<< p3.p[0] << p3.p[1] << p3.p[2] << -p3.p[3] << endr <<
							0 		<< 0.5*pymiss	<< 0 	<< 0 	<< 0 	<< 0.5*pymiss << 0 	<< 0 	<< endr;
					A4 = A4*2/Mnorm;
					C4 	<< 	2*minkowskidot(p5.p, p2.p) + 2*minkowskidot(p5.p, p7.p) + m5squared
						<<	2*minkowskidot(p2.p, p7.p) + m2squared
						<<	m7squared
						<<	pow(pxmiss, 2)
						<<	2*minkowskidot(p1.p, p6.p) + 2*minkowskidot(p1.p, p3.p) + m1squared
						<<	2*minkowskidot(p6.p, p3.p) + m6squared
						<< 	m3squared
						<< 	pow(pymiss, 2);
					C4 = C4/pow(Mnorm, 2);
	
					D41_list.push_back(inv(A4)*B);
					D42_list.push_back(inv(permute23*A4)*B);
					D43_list.push_back(inv(permute67*A4)*B);
					D44_list.push_back(inv(permute23and67*A4)*B);
					E41_list.push_back(inv(A4)*C4);
					E42_list.push_back(inv(permute23*A4)*C4);
					E43_list.push_back(inv(permute67*A4)*C4);
					E44_list.push_back(inv(permute23and67*A4)*C4);
				}
				else // same-flavour lepton chains, no need for 3 and 4.
				{
					mat dummy;
					dummy << 0 << endr;
					D31_list.push_back(dummy);
					D32_list.push_back(dummy);
					D33_list.push_back(dummy);
					D34_list.push_back(dummy);
					D41_list.push_back(dummy);
					D42_list.push_back(dummy);
					D43_list.push_back(dummy);
					D44_list.push_back(dummy);
					E31_list.push_back(dummy);
					E32_list.push_back(dummy);
					E33_list.push_back(dummy);
					E34_list.push_back(dummy);
					E41_list.push_back(dummy);
					E42_list.push_back(dummy);
					E43_list.push_back(dummy);
					E44_list.push_back(dummy);
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
	else cout << "Unable to open file" << endl; 

	// Store matrices and vectors compactly
	D_lists = {D11_list, D12_list, D13_list, D14_list, D21_list, D22_list, D23_list, D24_list, D31_list, D32_list, D33_list, D34_list, D41_list, D42_list, D43_list, D44_list};
	E_lists = {E11_list, E12_list, E13_list, E14_list, E21_list, E22_list, E23_list, E24_list, E31_list, E32_list, E33_list, E34_list, E41_list, E42_list, E43_list, E44_list};




	// Finished with making the D and E matrices. Now to minimize xisquared!


	// DEBUG evaluate xisquared
	// vector<double> masses_exact = {(565.312+570.734)/2, 180.337, 144.06, 97.0071979};
	// double debug_xisquared = xisquared(&masses_exact[0], 1, 0, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists);
	// cout << "DEBUG xisquared-true = " << debug_xisquared << endl;
	// DEBUG check detA bool (hack, stored in correct_combo)
	// cout << detAlist << endl;
	cout << correct_combinatorics << endl;

	// double tol = 0.1;
	// double maxiter = 500;
	for (int iBin=0; iBin<Nbins; iBin++)
	{
		cout << "Minimizing bin number " << iBin+1 << endl;
		double best_fit_current;
		vector<double> masses_current = masses_initial;
		// cout << xisquared(&masses_initial[0], Nevents, iBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists) << endl;
		// double *Masses_current = Masses_initial;




		double fmin;
		if (amoeba(&masses_current[0], fmin, xisquared, tol, maxiter, Nevents, iBin, Mnorm, combinatorics, all_leptons_equal_list, D_lists, E_lists, correct_combinatorics)) 
		{


			// Calculate masses from MLSP + squared-diffs
			double mllinv = 80; // Calculated from true masses using formula
			double MLSPsq = masses_current[2]*(masses_current[1]/(mllinv*mllinv) - 1.0);
			if (MLSPsq < 0)
				continue;
			vector<double> Masses = {0,0,0,0};
			Masses[3] = sqrt(MLSPsq);
			Masses[2] = sqrt(MLSPsq + masses_current[2]);
			Masses[1] = sqrt(MLSPsq + masses_current[2] + masses_current[1]);
			Masses[0] = sqrt(MLSPsq + masses_current[2] + masses_current[1] + masses_current[0]);
			if (std::isnan(Masses[3]) || std::isnan(Masses[2]) || std::isnan(Masses[1]) || std::isnan(Masses[0])) // Drop bin if any masses are NaN
				continue;

			best_fit_point.push_back(Masses);
			bin_number.push_back(iBin);

			best_fit_value.push_back(fmin);


			// HACK: Count detA passcut fraction, calling it correct_combinatorics
			vector<bool> correct_combinatorics_current;
			for (int iEvent = Nevents*iBin; iEvent < Nevents*(iBin+1); iEvent++)
			{
				correct_combinatorics_current.push_back(correct_combinatorics[iEvent]);
			}
			double fraction_of_correct_combinatorics = number_of_true(correct_combinatorics_current)/(double)Nevents;
			cout << fraction_of_correct_combinatorics << endl;
			correct_combinatorics_fraction.push_back(fraction_of_correct_combinatorics);
		}

		// cout << "correct_combinatorics = " << correct_combinatorics << endl;
		// double fraction_of_correct_combinatorics = number_of_true(correct_combinatorics)/(double)Nevents;
		// cout << fraction_of_correct_combinatorics << endl;
		// correct_combinatorics_fraction.push_back(fraction_of_correct_combinatorics);

	}



}



int main()
{


	int Nbins = 100;
	int Nevents = 25;
	bool combinatorics = false;
	vector<double> masses_initial = {568, 180, 144, 97};
	// vector<double> masses_initial = {400, 300, 200, 100};
	// vector<double> masses_initial = {800, 500, 300, 50};
	// vector<double> masses_initial = {1000, 100, 80, 30};
	double Mnorm = 100;
	double tol = 1e-12;
	double maxiter = 2000;

	// Calculate mass-squared diff vector from masses_initial
	vector<double> massdiff = {	masses_initial[0]*masses_initial[0] - masses_initial[1]*masses_initial[1], 
								masses_initial[1]*masses_initial[1] - masses_initial[2]*masses_initial[2],
								masses_initial[2]*masses_initial[2] - masses_initial[3]*masses_initial[3] };

	vector<double> best_fit_value;
	vector<vector<double> > best_fit_point; 
	vector<double> correct_combinatorics_fraction;
	vector<int> bin_number;


	string eventfile;
	// eventfile = "../events/simple_2500_events_no_mass_smearing.dat";
	// eventfile = "../events/simple_2500_events_gauss_and_exp_mass_smearing.dat";
	// eventfile = "../events/Pythia_cascade_events_no_ISR_or_FSR_20150120_only_opposite_flavour_leptons.dat";
	// eventfile = "../events/Pythia_cascade_10000_events_everything_turned_on_20150210_only_opposite_flavour_leptons.dat";
	// eventfile = "../events/herwigpp_only_OFL_20150305.dat";
	// eventfile = "../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL.dat";
	// eventfile = "../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL-10percent_momentum_smearing.dat";	
	eventfile = "../events/herwigpp-9563-events-complete-momcons-20150314_only_OFL-5percent_WEBBERmomentum_smearing.dat";
	// eventfile = "../events/HERWIG-events-10pmomsmear.dat";
	// eventfile = "../events/HERWIG-events.dat";



	best_fit(Nbins, Nevents, eventfile, masses_initial, tol, maxiter, combinatorics, Mnorm, best_fit_value, best_fit_point, correct_combinatorics_fraction, bin_number);

	cout << "correct_combinatorics_fraction = " << endl << correct_combinatorics_fraction << endl;

	int Naccepted = best_fit_value.size();

	for (int iBin = 0; iBin<Naccepted; iBin++)
	{
		cout << iBin+1 << "\t " << best_fit_value[iBin] << "\t ";
		cout << best_fit_point[iBin][0] << "\t " << best_fit_point[iBin][1] << "\t " << best_fit_point[iBin][2] << "\t " << best_fit_point[iBin][3] << "\t " << correct_combinatorics_fraction[iBin] << endl;
	}


	/** Make and open text output file */
	ofstream textOutput;
	textOutput.open("../best_fit_results/TEMP-subdetAcut.dat", ios::out);
	// textOutput.open("../best_fit_results/best_fit_100_bins_simple_combinatorics-OFF_massinit-571-181-145-98.dat", ios::out);

	textOutput << "# MASS DIFF FIT. combinatorics = " << combinatorics << " (true/false = 1/0)" << endl;
	textOutput << "# Event file name = " << eventfile << ", SIMPLEX tolerance = " << tol << endl;
	for (int iBin = 0; iBin < Naccepted; iBin++)
	{
		textOutput << iBin+1 << "\t" << best_fit_point[iBin][0] << "\t" << best_fit_point[iBin][1] << "\t" << best_fit_point[iBin][2] << "\t" << best_fit_point[iBin][3] << "\t " << 0 << "\t " << best_fit_value[iBin] << "\t" << correct_combinatorics_fraction[iBin] << "\t " << bin_number[iBin] << endl;
	}
	textOutput.close();

	cout << "Mean correct-combo fraction = " << accumulate(correct_combinatorics_fraction.begin(), correct_combinatorics_fraction.end(), 0.0)/correct_combinatorics_fraction.size() << endl;

	return 1;
}