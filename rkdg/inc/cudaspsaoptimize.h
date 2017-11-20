#pragma once
#include "cudarkdgsolver.h"
#include "kernelfunctions.cuh"
#include <random>
#include<math.h>
using namespace std;

class CSpsaOptimize {
private:
	int Niter;
	double SPSA_alpha;
	double SPSA_gamma;
	double SPSA_A;
	double SPSA_a;
	double SPSA_c;

	int IMAX;
	int DesignVariableNum;
	double radial;
	double opt_tol;

	CCUDARkdgSolver *solver;
	CUnstructuredGrid *grid;
	CTriangleInfo *triangle_infos;

	double *DesignVariable, *DesignVariable_perturb;
	double *initVertice;

	double distance2nodes2D(double *vert, int i, int j);
	double diffmax(double x);
public:
	CSpsaOptimize(CCUDARkdgSolver *solver);

	void RunOptimize();

	void MeshDeformationRBFI2D(double *verinit);

	inline double RBF(double x);

	void GaussElimination(int N, double **A, double *B);

	void CostFunction(bool ifpertuib, double cl0, double cd0, double area0, 
		double c1, double cd, double area, double &punfac1, double &punfac2, double &oj_val);

	void AccuracyTest_Airfoil(double &Cl, double &Cd);

	double hh(double x, int n, int nb);

	void hh_shape2para2D();

	void hh_para2shape2D(double *newvert, double *designVariable);

	void obj_val_output(int iter, const double *obj_val, double gradnorm, double cd, double cl, double area, double punfac1, double punfac2);

	void foil_output(int iter);

	double Uh(int e, double* ul, double x, double y);
};