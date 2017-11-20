#include "../inc/spsaoptimize.h"

CSpsaOptimize::CSpsaOptimize(CCUDARkdgSolver *solver):
Niter(50),
SPSA_alpha(1.0),
SPSA_gamma(0.167),
SPSA_A(300.0),
SPSA_a(0.167),
SPSA_c(1.0e-4),
IMAX(1),
DesignVariableNum(40),
radial(20.0),
opt_tol(1.0e-6)
{
	this->solver = solver;
	grid = &solver->grid;
	triangle_infos = &grid->triangle_infos;
	RunOptimize();
}

void CSpsaOptimize::RunOptimize() {

	ofstream of("output/air");
	/*for (int i = 0; i < grid->upper_node.size(); i++) {
		of<<grid->vertice[grid->upper_node.at(i)].getX()<<"  "<<grid->vertice[grid->upper_node.at(i)].getY()<<endl;
	}
	for (int i = grid->lower_node.size() - 1; i >= 0; i--) {
		of<<grid->vertice[grid->lower_node.at(i)].getX()<<"  "<<grid->vertice[grid->lower_node.at(i)].getY()<<endl;
	}*/
	/*for (int i = 0; i < grid->farfiled.size(); ++i) {
		of<<grid->vertice[grid->farfiled[i]].getX()<<"  "<<grid->vertice[grid->farfiled[i]].getY()<<endl;
	}*/
	int elem_l;
	int tnum = grid->getTriangleNumber();
	int num = grid->getCellNumber();
	for(int f = 0; f < grid->wall_elem.size(); ++f)
	{
		elem_l = grid->wall_elem[f];
		//loop for each GPs
		if (grid->tri_neighbour[elem_l] >= tnum) {
			of<<grid->vertice[grid->edge[grid->tri_edge[elem_l*TRIANGLE_EDGES]].getStart()].getX()<<"  "
				<<grid->vertice[grid->edge[grid->tri_edge[elem_l*TRIANGLE_EDGES]].getStart()].getY()<<endl;
		} else if (grid->tri_neighbour[elem_l + num] >= tnum) {
			of<<grid->vertice[grid->edge[grid->tri_edge[elem_l*TRIANGLE_EDGES+1]].getStart()].getX()<<"  "
				<<grid->vertice[grid->edge[grid->tri_edge[elem_l*TRIANGLE_EDGES+1]].getStart()].getY()<<endl;
		} else {
			of<<grid->vertice[grid->edge[grid->tri_edge[elem_l*TRIANGLE_EDGES+2]].getStart()].getX()<<"  "
				<<grid->vertice[grid->edge[grid->tri_edge[elem_l*TRIANGLE_EDGES+2]].getStart()].getY()<<endl;
		}
	
	}
	of.close();
	cout<<"finish"<<endl;
	getchar();
	random_device rd;
	mt19937 gen(rd());
	bernoulli_distribution distribution(0.5);

	int iter = 0;
	double ak, ck;
	double *delta_k, *delta_k_minus;
	vector<CVertice2D> newvert;
	bool ifperturb;
	double *obj_val = new double[IMAX];
	double cl0(0), cd0(0), area0(0), cl(0), cd(0), area(0),gradnorm(0.0), obj_val_plus(0.0), obj_val_minus(0.0);
	double punfac1(5.0), punfac2(5.0);
	double *grad;
	delta_k = new double[DesignVariableNum];
	delta_k_minus = new double[DesignVariableNum];
	grad = new double[DesignVariableNum];
	DesignVariable = new double[DesignVariableNum];
	DesignVariable_perturb = new double[DesignVariableNum];
	newvert.resize(grid->getVerticeNumber());
	initVertice.resize(grid->getVerticeNumber());

	for (int i = 0; i < grid->getVerticeNumber(); ++i) {
		initVertice[i].setVertice(grid->vertice[i].getX(), grid->vertice[i].getY());
	}

	for (int i = 0; i < DesignVariableNum; ++i) {
		DesignVariable[i] = 0.0;
	}
	AccuracyTest_Airfoil(cl0, cd0);
	cout<<"cl: "<<cl0<<"    cd: "<<cd0;
	getchar();
	area0 = grid->calArea();
	obj_val[0] = 1.0;
	foil_output(iter);
	obj_val_output(iter, obj_val, gradnorm, cd0, cl0, area0, punfac1, punfac2);
	for (iter = 1; iter < IMAX; iter++) {

		for (int i = 0; i < grid->getVerticeNumber(); ++i) {
			newvert[i].setVertice(grid->vertice[i].getX(), grid->vertice[i].getY());
		}

		ck  = SPSA_c / std::pow((double)iter, SPSA_gamma);
		ak  = SPSA_a / std::pow((double)(iter+SPSA_A), SPSA_alpha);
		for(int i = 0; i < DesignVariableNum; ++i)
		{
			delta_k[i] = (double)(2 * distribution(gen) - 1);
			delta_k_minus[i] = 1.0 / delta_k[i];
		}
		ifperturb = true;
		for(int i = 0; i < DesignVariableNum; ++i) {
			DesignVariable_perturb[i] = DesignVariable[i] + ck * delta_k[i];
		}
		hh_para2shape2D(initVertice, DesignVariable_perturb);
		MeshDeformationRBFI2D(newvert);
		solver->runNext();
		AccuracyTest_Airfoil(cl, cd);
		area = grid->calArea();
		CostFunction(ifperturb, cl0, cd0, area0, cl, cd, area, punfac1, punfac2, obj_val_plus);
		for (int i = 0; i < DesignVariableNum; ++i) {
			DesignVariable_perturb[i] = DesignVariable[i] - ck * delta_k[i];
		}
		hh_para2shape2D(newvert, DesignVariable_perturb);
		MeshDeformationRBFI2D(newvert);
		solver->runNext();
		AccuracyTest_Airfoil(cl, cd);
		area = grid->calArea();
		CostFunction(true, cl0, cd0, area0, cl, cd, area, punfac1, punfac2, obj_val_minus);

		ifperturb = false;
		double coef = (obj_val_plus - obj_val_minus) / (2.0 * ck);
		for (int i = 0; i < DesignVariableNum; ++i) {
			grad[i] = coef * delta_k_minus[i];
		}
		gradnorm = 0.0;
		for (int i = 0; i < DesignVariableNum; ++i) {
			gradnorm = gradnorm + grad[i] * grad[i];
		}
		gradnorm = std::sqrt(gradnorm);
		for (int i = 0; i < DesignVariableNum; ++i) {
			DesignVariable[i] = DesignVariable[i] - ak * grad[i] / gradnorm;
		}

		hh_para2shape2D(newvert, DesignVariable);
		MeshDeformationRBFI2D(newvert);
		solver->runNext();
		AccuracyTest_Airfoil(cl, cd);
		area = grid->calArea();
		CostFunction(ifperturb, cl0, cd0, area0, cl, cd, area, punfac1, punfac2, obj_val[iter]);
		foil_output(iter);
		obj_val_output(iter, obj_val, gradnorm, cd, cl, area, punfac1,punfac2);
		
		if (iter >= Niter) {
			int j;
			for (j = 1; j <= Niter; ++j) {
				if( std::fabs(obj_val[iter] - obj_val[iter-j]) >= opt_tol)
					break;
			}
			if(j == Niter+1)
			{
				std::cout << "convergence at step " << iter << std::endl;
				return ;
			}
		}
	}
	solver->outputSolution("output/solution2.plt");
}

void CSpsaOptimize::AccuracyTest_Airfoil(double &Cl, double &Cd) {
	int elem_l(0),face_l(0);
	int node0(0),node1(0);
	double gpx(0), gpy(0), gpw(0);
	double nx(0),ny(0);
	double forcePX(0),forcePY(0),forceX(0),forceY(0);
	double rho(0),u(0),v(0),pre(0),rhom(0),rhou(0),rhov(0), rhoE(0);
	
	int cell_num = grid->getCellNumber();
	int tnum = grid->getTriangleNumber();
	this->grid = &solver->grid;
	double  ut = solver->gamma*solver->pref/solver->rhoref*solver->mach*solver->mach;

	Cl = 0.0;
	Cd = 0.0;

	double q_inf = 0.5 * solver->rhoref * ut;

	forcePX = 0;
	forcePY = 0;
	int gauss_index;
	for(int f = 0; f < grid->wall_elem.size(); ++f)
	{
		elem_l = grid->wall_elem[f];
		//loop for each GPs
		for(int p = 0; p < EDGE_GPOINTS; ++p)
		{
			if (grid->tri_neighbour[elem_l] >= tnum) {
				gauss_index = p;
				nx = triangle_infos->outer_normal_vector[elem_l];
				ny = triangle_infos->outer_normal_vector[elem_l+1*cell_num];
			} else if (grid->tri_neighbour[elem_l+cell_num] >= tnum) {
				gauss_index = p + EDGE_GPOINTS;
				nx = triangle_infos->outer_normal_vector[elem_l+2*cell_num];
				ny = triangle_infos->outer_normal_vector[elem_l+3*cell_num];
			} else if (grid->tri_neighbour[elem_l+2*cell_num] >= tnum) {
				gauss_index = p + 2 * EDGE_GPOINTS;
				nx = triangle_infos->outer_normal_vector[elem_l+4*cell_num];
				ny = triangle_infos->outer_normal_vector[elem_l+5*cell_num];
			}
			
			gpw = triangle_infos->edge_gauss_weight[elem_l+gauss_index*cell_num];

			//calculate physical variables
			rho  = Uh(elem_l, solver->_freedom_rho, gauss_index);
			rhou = Uh(elem_l, solver->_freedom_rhou, gauss_index);
			rhov = Uh(elem_l, solver->_freedom_rhov, gauss_index);
			rhoE = Uh(elem_l, solver->_freedom_rhoE, gauss_index);

			rhom = 1.0/rho;
			u    = rhou * rhom;
			v    = rhov * rhom;

			pre  = (solver->gamma - 1.0)*(rhoE - 0.5 * rho * (u*u + v*v));

			//Force by pressure
			forcePX += gpw * pre * nx;
			forcePY += gpw * pre * ny;
		}
	}
	forceX = forcePX;
	forceY = forcePY;
	Cl = forceY * cos(solver->alpha) - forceX * sin(solver->alpha);
	Cd = forceY * sin(solver->alpha) + forceX * cos(solver->alpha);

	Cl = Cl / q_inf;
	Cd = Cd / q_inf;
}


//void CSpsaOptimize::AccuracyTest_Airfoil(double &Cl, double &Cd) {
//	int elem_l(0),face_l(0);
//	int node0(0),node1(0);
//	double gpx(0), gpy(0), gpw(0);
//	double nx(0),ny(0);
//	double forcePX(0),forcePY(0),forceX(0),forceY(0);
//	double rho(0),u(0),v(0),pre(0),rhom(0),rhou(0),rhov(0), rhoE(0);
//
//	int cell_num = grid->getCellNumber();
//	int tnum = grid->getTriangleNumber();
//	this->grid = &solver->grid;
//	double  ut = sqrt(solver->gamma*solver->pref/solver->rhoref)*solver->mach;
//
//	Cl = 0.0;
//	Cd = 0.0;
//
//	double q_inf = 0.5 * solver->rhoref * ut * ut;
//
//	forcePX = 0;
//	forcePY = 0;
//
//	for(int f = 0; f < grid->wall_elem.size(); ++f)
//	{
//		elem_l = grid->wall_elem[f];
//		//loop for each GPs
//		for(int p = 0; p < EDGE_GPOINTS; ++p)
//		{
//			if (grid->tri_neighbour[elem_l] >= tnum) {
//				gpx = triangle_infos->edge_gauss_vertice[elem_l*EDGE_GPOINTS*TRIANGLE_EDGES+p].getX();
//				gpy = triangle_infos->edge_gauss_vertice[elem_l*EDGE_GPOINTS*TRIANGLE_EDGES+p].getY();
//				gpw = triangle_infos->edge_gauss_weight[elem_l+p*cell_num];
//				nx = triangle_infos->outer_normal_vector[elem_l];
//				ny = triangle_infos->outer_normal_vector[elem_l+1*cell_num];
//			} else if (grid->tri_neighbour[elem_l+cell_num] >= tnum) {
//				gpx = triangle_infos->edge_gauss_vertice[elem_l*EDGE_GPOINTS*TRIANGLE_EDGES+ EDGE_GPOINTS + p].getX();
//				gpy = triangle_infos->edge_gauss_vertice[elem_l*EDGE_GPOINTS*TRIANGLE_EDGES+ EDGE_GPOINTS + p].getY();
//				gpw = triangle_infos->edge_gauss_weight[elem_l+(EDGE_GPOINTS + p)*cell_num];
//				nx = triangle_infos->outer_normal_vector[elem_l+2*cell_num];
//				ny = triangle_infos->outer_normal_vector[elem_l+3*cell_num];
//			} else if (grid->tri_neighbour[elem_l+2*cell_num] >= tnum) {
//				gpx = triangle_infos->edge_gauss_vertice[elem_l*EDGE_GPOINTS*TRIANGLE_EDGES+ 2*EDGE_GPOINTS + p].getX();
//				gpy = triangle_infos->edge_gauss_vertice[elem_l*EDGE_GPOINTS*TRIANGLE_EDGES+ 2*EDGE_GPOINTS + p].getY();
//				gpw = triangle_infos->edge_gauss_weight[elem_l+(2*EDGE_GPOINTS + p)*cell_num];
//				nx = triangle_infos->outer_normal_vector[elem_l+4*cell_num];
//				ny = triangle_infos->outer_normal_vector[elem_l+5*cell_num];
//			}
//
//			//calculate physical variables
//			rho  = Uh(elem_l, solver->_freedom_rho, gpx, gpy);
//			rhou = Uh(elem_l, solver->_freedom_rhou, gpx, gpy);
//			rhov = Uh(elem_l, solver->_freedom_rhov, gpx, gpy);
//			rhoE = Uh(elem_l, solver->_freedom_rhoE, gpx, gpy);
//
//			rhom = 1.0/rho;
//			u    = rhou * rhom;
//			v    = rhov * rhom;
//
//			pre  = (solver->gamma - 1.0)*(rhoE - 0.5 * rho * (u*u + v*v));
//
//			//Force by pressure
//			forcePX += gpw * pre * nx;
//			forcePY += gpw * pre * ny;
//		}
//	}
//	forceX = forcePX;
//	forceY = forcePY;
//	Cl = forceY * cos(solver->alpha) - forceX * sin(solver->alpha);
//	Cd = forceY * sin(solver->alpha) + forceX * cos(solver->alpha);
//
//	Cl = Cl / q_inf;
//	Cd = Cd / q_inf;
//}


void CSpsaOptimize::MeshDeformationRBFI2D(vector<CVertice2D> verinit)
{
	int nwi, nwj, nii;
	double disp, dispx, dispy;
	double *disx, *disy;
	double *coefx, *coefy;
	int wallNum = grid->wall.size();
	int farfiledNum = grid->farfiled.size();
	int N = wallNum + farfiledNum + 3;
	double **matrix;
	matrix = (double **) malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		matrix[i] = (double *)malloc(N * sizeof(double));
	}
	disx = (double *) malloc(N * sizeof(double));
	disy = (double *) malloc(N * sizeof(double));
	coefx = (double *) malloc(N * sizeof(double));
	coefy = (double *) malloc(N * sizeof(double));
	for (int i = 0; i < N; ++i) {
		disx[i] = 0.0;
		disy[i] = 0.0;
	}
	for (int i = 0; i < N; i++) {
		for (int j = 0; j < N; j++) {
			matrix[i][j] = 0;
		}
	}

	for (int i = 0; i < wallNum + farfiledNum; ++i) {
		if (i < wallNum) {
			nwi = grid->wall[i];
		} else {
			nwi = grid->farfiled[i-wallNum]; 
		}
		for (int j = i; j < wallNum + farfiledNum; ++j) {
			if(j < wallNum)
				nwj  = grid->wall[j]; 
			else
				nwj  = grid->farfiled[j-wallNum];
			disp = distance2nodes2D(verinit, nwi, nwj);
			if(j == i) {
				matrix[i][j] = RBF(disp);
			} else{
				matrix[i][j] = RBF(disp);
				matrix[j][i] = RBF(disp);
			}
		}
		matrix[i][N-3]   = 1.0;   //P_b and P_b'
		matrix[N-3][i]   = 1.0;
		
		matrix[i][N-2] = verinit[nwi].getX();
		matrix[N-2][i] = verinit[nwi].getX();

		matrix[i][N-1] = verinit[nwi].getY();
		matrix[N-1][i] = verinit[nwi].getY();
	}
	
	// Step 2: Calculate the displacement, new - old (only wall boundary moves)
	for(int i = 0; i < wallNum; ++i){
		nwi = grid->wall[i];
		disx[i] = grid->vertice[nwi].getX() - verinit[nwi].getX();
		disy[i] = grid->vertice[nwi].getY() - verinit[nwi].getY();
	}
	// Step 3: Solve the equations for the coef
	GaussElimination(N, matrix, disx);
	GaussElimination(N, matrix, disy); //disx / disy is the solution returned
	coefx = disx;
	coefy = disy;
	// Step 4: Generate new inner nodes coordinate
	for (int i = 0; i < grid->getVerticeNumber(); ++i) {
		if (grid->ver_flag[i] == VERTICE_INTERIOR) {
			dispx = 0.0;
			dispy = 0.0;
			for(int j = 0; j < wallNum+farfiledNum; ++j){
				if(j < wallNum)
					nwj = grid->wall[j];
				else
					nwj = grid->farfiled[j-wallNum];
				disp  = distance2nodes2D(verinit, i, nwj);
				dispx = dispx + coefx[j] * RBF(disp);
				dispy = dispy + coefy[j] * RBF(disp);
			}
			dispx = dispx + 1 * coefx[N-3] + verinit[i].getX() * coefx[N-2] + verinit[i].getY() * coefx[N-1];
			dispy = dispy + 1 * coefy[N-3] + verinit[i].getX() * coefy[N-2] + verinit[i].getY() * coefy[N-1];
			
			grid->vertice[i].setX(verinit[i].getX() + dispx);
			grid->vertice[i].setY(verinit[i].getY() + dispy);
		}
	}
	// Step 5: Far field nodes does not move
}

double CSpsaOptimize::distance2nodes2D(vector<CVertice2D> vert, int i, int j) {
	double dis(0.0);
	dis = (vert[i].getX() - vert[j].getX())
			* (vert[i].getX() -vert[j].getX()) +
		  (vert[i].getY() - vert[j].getY())
		    * (vert[i].getY() - vert[j].getY());
	dis = sqrt(dis);
	return dis;
}

double CSpsaOptimize::RBF(double x) {
	if(x >= 0.0 && x <= radial)
		return (1.0 - x/radial) * (1.0 - x/radial);
	else 
		return 0.0;
}

void CSpsaOptimize::GaussElimination(int N, double **A, double *B) {
	int i(0), j(0), k(0), n(0);
	double **a;
	double *b, *x;
	a = (double **) malloc(N * sizeof(double*));
	for (int i = 0; i < N; i++) {
		a[i] = (double *)malloc(N * sizeof(double));
	}
	b = (double*) malloc(N * sizeof(double));
	x = (double*) malloc(N * sizeof(double));
	for (i = 0; i < N; ++i) {
		for (int j = 0; j < N; ++j) {
			a[i][j] = A[i][j];
		}
	}

	for (j = 0; j < N; ++j) {
		b[j] = B[j];
	}

	int num(0);
	double maj(0.0);
	double temp(0.0), sum(0.0);
	double *m;
	m = (double *) malloc(N * sizeof(double));
	for (k = 0; k < N - 1; ++k) {
		num = k;
		maj = a[k][k];
		for (i = k + 1; i < N; ++i) {
			if (fabs(a[i][k]) > fabs(maj)) {
				num = i;
				maj = a[i][k];
			}
		}

		for (j = k; j < N; ++j) {
			temp = a[k][j];
			a[k][j] = a[num][j];
			a[num][j] = temp;
		}
		temp = b[k];
		b[k] = b[num];
		b[num] = temp;

		for (i = k + 1; i < N; ++i) {
			m[i] = a[i][k] / a[k][k];

			for (j = k + 1; j < N; ++j) {
				a[i][j] = a[i][j] - m[i] * a[k][j];
			}
			b[i] = b[i] - m[i] * b[k];
		}
	}

	x[N - 1] = b[N - 1] / a[N - 1][N - 1];
	for (n = N - 2; n > -1; --n) {
		sum = 0.0;
		for (j = n + 1; j < N; j++) {
			sum = sum + a[n][j] * x[j];
		}
		
		x[n] = (b[n] - sum) / a[n][n];
	}

	for (i = 0; i < N; ++i) {

		
		B[i] = x[i];
	}
}

void CSpsaOptimize::CostFunction(bool ifpertuib, double cl0, double cd0, double area0, 
		double cl, double cd, double area, double &punfac1, double &punfac2, double &obj_val) {
	double dcl, darea;
	switch(1) {
	case 1:
		dcl = 1.0 - cl / cl0;
		darea = 1.0 - area / area0;
		obj_val = cd / cd0 + punfac1 * diffmax(dcl) + punfac2 * diffmax(darea);
		break;
	case 2:
		dcl     = 1.0 - cl/cl0;
		darea   = 1.0 - area/area0;
		obj_val = 0.0 - (cl/cd)/(cl0/cd0) + punfac1*diffmax(dcl) + punfac2*diffmax(darea);
		break;

	default:
		cout << "OptCase is not supported!" << std::endl;
		exit(1);
	}
	return;
}

double CSpsaOptimize::diffmax(double x) {
	double eps = 1.0e-6;
	if(x > eps)
		return x;
	else if(x < 0)
		return 0;
	else
		return (x*x)/(2*eps);	
}

double CSpsaOptimize::hh(double x, int n, int nb) {
	double pi = 4.0 * std::atan(1.0);
	double fx;
	double t, pn, bn;

	//pn = 0.5 * (1.0 - std::cos((double)n*pi / (double)(nb+1))); //两头密，中间稀疏
	pn = 1.0 / (double)nb * ((double)n - 0.5);  //均匀分布
	t  = 3.0;                 // t should be an integer
	bn = log(0.5) / log(pn);
	fx = pow(sin(pi * pow(x,bn)), t);

	if(n == nb)
		fx = 8.0 * x * (1.0-x) * exp(-10.0 * (1.0-x)); // reset at trailing edge(?)

	return fx;
}

void CSpsaOptimize::hh_para2shape2D(vector<CVertice2D> vertemp, double *designVariable) {
	int upperbasisNum = DesignVariableNum / 2;
	int lowerbasisNum = DesignVariableNum - upperbasisNum;
	int ni, upperNum, lowerNum;

	upperNum = grid->upper_node.size();
	for(int i = 0; i < upperNum; ++i){
		ni = grid->upper_node.at(i);
		grid->vertice[ni].setY(vertemp[ni].getY());
		for(int j = 1; j <= upperbasisNum; ++j) {
			grid->vertice[ni].setY(grid->vertice[ni].getY() + designVariable[j - 1] 
				* hh(vertemp[ni].getX(), j, upperbasisNum));
		}
	}
	
	lowerNum = grid->lower_node.size();
	for(int i = 0; i < lowerNum; ++i){
		ni = grid->lower_node.at(i);
		grid->vertice[ni].setY(vertemp[ni].getY());
		for(int j = upperbasisNum+1; j <= DesignVariableNum; ++j) {
			grid->vertice[ni].setY(grid->vertice[ni].getY() + designVariable[j - 1]
				* hh(vertemp[ni].getX(),j - upperbasisNum, lowerbasisNum));
		}
	}
}

void CSpsaOptimize::obj_val_output(int iter, const double *obj_val, double gradnorm, double cd, double cl, double area, double punfac1, double punfac2) {
	ofstream outfile;
	if (iter == 0) {
		outfile.open("output/obj_val.plt", ios::out);
		outfile<< "Title = Objective Value History"<<endl;
		outfile<< "Variables=iter, obj_val,  gradnorm,   Cd,   Cl,  area, punfac1, punfac2"<<std::endl;
	}
	else {
		outfile.open("output/obj_val.plt", ios::app);
	}

	outfile << std::setw(6)  << iter;
	outfile << std::setiosflags(ios::fixed);
	outfile << std::setw(10) << obj_val[iter]; 
	outfile << std::setw(15) << gradnorm;
	outfile << std::setw(10) << cd;
	outfile << std::setw(10) << cl;
	outfile << std::setw(10) << area;
	outfile << std::resetiosflags(ios::fixed);
	outfile << std::setw(10) << punfac1;
	outfile << std::setw(10) << punfac2 << std::endl;
	outfile.close();
}

void CSpsaOptimize::foil_output(int iter) {
	int upperNum = grid->upper_node.size();
	int lowerNum = grid->lower_node.size();
	ofstream outfile;
	if(iter == 0)
	{
		outfile.open("output/foil.plt",ios::out);
		outfile << "Title = airfoil" << endl;
		outfile << "Variables=X,Y" << endl;
	}
	else
	{
		outfile.open("output/foil.plt",ios::app);
	}

	outfile << "ZONE T=T" << iter << "  F=POINT" << std::endl;
	for (int i = 0; i < grid->upper_node.size(); i++) {
		outfile<<grid->vertice[grid->upper_node.at(i)].getX()<<"  "<<grid->vertice[grid->upper_node.at(i)].getY()<<endl;
	}
	for (int i = grid->lower_node.size() - 1; i >= 0; i--) {
		outfile<<grid->vertice[grid->lower_node.at(i)].getX()<<"  "<<grid->vertice[grid->lower_node.at(i)].getY()<<endl;
	}
	outfile.close();
}

double CSpsaOptimize::Uh(int e,double* ul, int gauss_index) {
	int cell_num = grid->getCellNumber();
	double s(0);
	s = ul[e]
		+ ul[e + cell_num] * triangle_infos->edge_bf_value[e+gauss_index*cell_num+1*TRIANGLE_EDGES*EDGE_GPOINTS*cell_num]
		+ ul[e + 2*cell_num] * triangle_infos->edge_bf_value[e+gauss_index*cell_num+2*TRIANGLE_EDGES*EDGE_GPOINTS*cell_num]
		+ ul[e + 3*cell_num] * triangle_infos->edge_bf_value[e+gauss_index*cell_num+3*TRIANGLE_EDGES*EDGE_GPOINTS*cell_num]
		+ ul[e + 4*cell_num] * triangle_infos->edge_bf_value[e+gauss_index*cell_num+4*TRIANGLE_EDGES*EDGE_GPOINTS*cell_num]
		+ ul[e + 5*cell_num] * triangle_infos->edge_bf_value[e+gauss_index*cell_num+5*TRIANGLE_EDGES*EDGE_GPOINTS*cell_num];
	return s;
}


//double CSpsaOptimize::Uh(int e,double* ul, double x, double y) {
//	int cell_num = grid->getCellNumber();
//	double s(0);
//	s = ul[e]
//	+ ul[e + cell_num] * triangle_infos->basisFunction(e, 1, x, y)
//		+ ul[e + 2*cell_num] * triangle_infos->basisFunction(e, 2, x, y)
//		+ ul[e + 3*cell_num] * triangle_infos->basisFunction(e, 3, x, y)
//		+ ul[e + 4*cell_num] * triangle_infos->basisFunction(e, 4, x, y)
//		+ ul[e + 5*cell_num] * triangle_infos->basisFunction(e, 5, x, y);
//	return s;
//}