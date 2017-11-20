#include "../inc/cudarkdgsolver.h"

CCUDARkdgSolver::CCUDARkdgSolver():
title("Unknown case"),
alpha(5.0),
gamma(1.4),
mach(0.4),
cfl(0.18),
rhoref(1.0),
pref(1.0),

_terminal_time(15),
log_history('Y'),
print_interval(1000),

gridconf("input/mesh.conf"),
solution_file("output/solution.dat"),
log_file("output/log.dat"),
residual_file("output/residual.dat"),

threads_per_block(512),
reduction_threads(512),

_freedom_rho(NULL),
_freedom_rhou(NULL),
_freedom_rhov(NULL),
_freedom_rhoE(NULL),

_dt(NULL),
_residual(NULL)
{}


void CCUDARkdgSolver::detectCUDADevice( void )
{
	int count(0);

	cudaGetDeviceCount( &count );

	if ( 0==count )
		throw CMyException("No device surpports CUDA found!");

	cudaDeviceProp prop;

	bool double_support(false);

	for ( int i=0; i<count; ++i )
	{
		cudaGetDeviceProperties( &prop, i );
		if ( prop.major>1 )
		{
			double_support = true;
			break;
		}
	}

	if ( !double_support )
		throw CMyException("No device has capability of 2.0 or higher is found!");

	memset( &prop, 0, sizeof(cudaDeviceProp) );
	prop.major = 2;
	prop.minor = 0;

	int devid;

	cudaChooseDevice(&devid, &prop);

	cout<<"\nThere are "<<count<<" device surpports CUDA, and the "<<devid+1<<"th device will be used."<<endl;
}

void CCUDARkdgSolver::initConfig(void)
{
	string conf_items[] = {
		"title",
		"gamma", "alpha", "mach", "cfl", "rhoref","pref",
		"time",  
		"gridconf", "logfile",  "solutionfile", "residualfile", 
		"threadsperblock", "reductionthreads",
		"loghistory", "printinterval"
	};

	CConfig program_conf(config_file, conf_items, 16);

	program_conf.parseConfigFile();

	// ת�����ò���
	if ( program_conf.config_items["title"]!="" )
		title = program_conf.config_items["title"];

	if ( program_conf.config_items["gamma"]!="" )
		gamma = atof(program_conf.config_items["gamma"].c_str());
		
	if ( program_conf.config_items["alpha"]!="" )
		alpha = atof(program_conf.config_items["alpha"].c_str())*atan(1.0)*4 / 180;
		
	if ( program_conf.config_items["mach"]!="" )
		mach  = atof(program_conf.config_items["mach"].c_str());
	
	if ( program_conf.config_items["cfl"]!="" )
		cfl   = atof(program_conf.config_items["cfl"].c_str());
		
	if ( program_conf.config_items["rhoref"]!="" )
		rhoref = atof(program_conf.config_items["rhoref"].c_str());
	
	if ( program_conf.config_items["pref"]!="" )
		pref = atof(program_conf.config_items["pref"].c_str());

	if ( program_conf.config_items["time"]!="" )
		_terminal_time = atof(program_conf.config_items["time"].c_str());

	if ( program_conf.config_items["gridconf"]!="" )
		gridconf = program_conf.config_items["gridconf"];

	if ( program_conf.config_items["solutionfile"]!="" )
		solution_file = program_conf.config_items["solutionfile"];
		
	if ( program_conf.config_items["logfile"]!="" )
		log_file = program_conf.config_items["logfile"];

	if ( program_conf.config_items["residualfile"]!="" )
		residual_file = program_conf.config_items["residualfile"];
		
	if ( program_conf.config_items["loghistory"]!="" )
		log_history = toupper(program_conf.config_items["loghistory"].at(0));
		
	if ( program_conf.config_items["threadsperblock"]!="" )	
		threads_per_block = atoi(program_conf.config_items["threadsperblock"].c_str());

	if ( program_conf.config_items["printinterval"]!="" )	
		print_interval = abs(atoi(program_conf.config_items["printinterval"].c_str()));
	
	if ( program_conf.config_items["reductionthreads"]!="" )
		reduction_threads = atoi(program_conf.config_items["reductionthreads"].c_str());
}

void CCUDARkdgSolver::run(void)
{
	ofstream fout(log_file.c_str());
	if ( !fout )
		throw CMyException("Failed to open log file: "+log_file);
	
	CMyTime mt;
	
	
	fout<<mt.getCurrentTime()<<": programs starts"<<endl;

	// ���CUDA�豸
	detectCUDADevice();
	fout<<mt.getCurrentTime()<<": Device with capability of 2.0 is found."<<endl;

	
	// ��ʼ���������ò������������
	initConfig();
	fout<<mt.getCurrentTime()<<": initialize configure from file."<<endl<<endl;
	fout<<"Title: "<<title<<endl<<endl;
	printConfig(cout);
	printConfig(fout);
	
	
	fout<<mt.getCurrentTime()<<": reading grid information."<<endl;

	// ��ʼ������
	grid.config_file = gridconf;
	grid.initializeGrid();
	grid.outputGrid();
	grid.outputGridWithGhostCells("output/ghostmesh.plt");
	
	fout<<mt.getCurrentTime()<<": complete grid initialization."<<endl;

	// ���������ζ����Ƿ���ʱ������
	grid.testTrianglesAntiwise();
	
	// ��ʼ�������ϻ���������Ϣ
	grid.triangle_infos.allocateMemory(grid.getCellNumber());
	grid.initializeTriangleInfos();
	
	fout<<mt.getCurrentTime()<<": complete grid information initialization."<<endl;

	// �������Ԫ
	grid.markBoundaryTriangles();
	
	// ����GPU�ڴ�
	_cuarrays.allocateMemory(grid.getCellNumber(), grid.wall_elem.size());

	// �����ǵ�Ԫ��Ϣ���͵�GPU
	copyTriangleInfosToGPU();

	// ��ʼ��RKDG���ɶȣ�������ʼ�����ݴ���GPU
	initRKDG();

	fout<<mt.getCurrentTime()<<": program initialization complete."<<endl;

	fout<<mt.getCurrentTime()<<": begin to solve flow."<<endl<<endl;
	/** ʱ���ƽ�*/
	mt.beginTimer();
	//rkdgAdvance();
	cudaDeviceSynchronize();
	mt.endTimer();
	
	fout<<"RKDG performance:"<<endl;
	fout<<"CPU time:  "<<mt.getCPUElapsedTime()<<" s"<<endl;
	fout<<"wall time: "<<mt.getWallElapsedTime()<<" s"<<endl<<endl;
	
	fout<<mt.getCurrentTime()<<": complete solving flow."<<endl;

	// �������ɶȵ�����
	copyFreedomToHost();
	
	// �����
	outputSolution(solution_file);

	fout<<mt.getCurrentTime()<<": complete solution output."<<endl;
	
	fout.close();
}

void CCUDARkdgSolver::runNext() {
	grid.initGhostGrid();
	grid.initializeTriangleInfos();
	copyTriangleInfosToGPU();
	//initRKDG();
	rkdgAdvance();
	copyFreedomToHost();
}

void CCUDARkdgSolver::copyFreedomToHost()
{
	size_t size = sizeof(double)*grid.getCellNumber();
	size_t pitch = _cuarrays.getDoublePitch();

	/*cudaMemcpy2D(_freedom_rho,  size, _cuarrays.freedom_rho,  pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
	cudaMemcpy2D(_freedom_rhou, size, _cuarrays.freedom_rhou, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			  
	cudaMemcpy2D(_freedom_rhov, size, _cuarrays.freedom_rhov, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);			 
	cudaMemcpy2D(_freedom_rhoE, size, _cuarrays.freedom_rhoE, pitch, size, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

	cudaDeviceSynchronize();*/

	ifstream ifs("input/rho.dat");
	double temp;
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j  < grid.getCellNumber(); j++) {
			ifs>>_freedom_rho[j+i*grid.getCellNumber()];
		}
	}
	ifs.close();

	ifs.open("input/rhou.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j  < grid.getCellNumber(); j++) {
			ifs>>_freedom_rhou[j+i*grid.getCellNumber()];
		}
	}
	ifs.close();

	ifs.open("input/rhov.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j  < grid.getCellNumber(); j++) {
			ifs>>_freedom_rhov[j+i*grid.getCellNumber()];
		}
	}
	ifs.close();

	ifs.open("input/rhoE.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j  < grid.getCellNumber(); j++) {
			ifs>>_freedom_rhoE[j+i*grid.getCellNumber()];
		}
	}
	ifs.close();

	cudaMemcpy2D(_cuarrays.freedom_rho,  pitch, _freedom_rho,  size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2D(_cuarrays.freedom_rhou, pitch, _freedom_rhou, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			  
	cudaMemcpy2D(_cuarrays.freedom_rhov, pitch, _freedom_rhov, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);			 
	cudaMemcpy2D(_cuarrays.freedom_rhoE, pitch, _freedom_rhoE, size, size, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);

}

void CCUDARkdgSolver::copyTriangleInfosToGPU(void)
{
	int num = grid.getCellNumber();
	size_t int_pitch    = _cuarrays.getIntPitch();
	size_t double_pitch = _cuarrays.getDoublePitch();
	
	cudaMemcpy2DAsync(_cuarrays.neighbour, int_pitch, grid.tri_neighbour, sizeof(int)*num, sizeof(int)*num, TRIANGLE_EDGES, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.sharedEdge, int_pitch, grid.tri_sharedEdge, sizeof(int)*num, sizeof(int)*num, TRIANGLE_EDGES, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.triangle_flag, int_pitch, grid.tri_flag, sizeof(int)*num, sizeof(int)*num, 1, cudaMemcpyHostToDevice);

	size_t gsize = sizeof(double)*num;

	cudaMemcpy2DAsync(_cuarrays.area, double_pitch, grid.triangle_infos.area, gsize, gsize, 1, cudaMemcpyHostToDevice);
	
	cudaMemcpy2DAsync(_cuarrays.perimeter, double_pitch, grid.triangle_infos.perimeter, gsize, gsize, 1, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.outer_normal_vector, double_pitch, grid.triangle_infos.outer_normal_vector, gsize, gsize, TRIANGLE_EDGES*2, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.mass_coeff, double_pitch, grid.triangle_infos.mass_coeff, gsize, gsize, BASIS_FUNCTIONS, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.vol_bf_value, double_pitch, grid.triangle_infos.vol_bf_value, gsize, gsize, VOLUME_GPOINTS*BASIS_FUNCTIONS, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.vol_bdf_value, double_pitch, grid.triangle_infos.vol_bdf_value, gsize, gsize, VOLUME_GPOINTS*BASIS_FUNCTIONS*2, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.edge_bf_value, double_pitch, grid.triangle_infos.edge_bf_value, gsize, gsize, TRIANGLE_EDGES*EDGE_GPOINTS*BASIS_FUNCTIONS, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.vol_gauss_weight, double_pitch, grid.triangle_infos.vol_gauss_weight, gsize, gsize,  VOLUME_GPOINTS, cudaMemcpyHostToDevice);

	cudaMemcpy2DAsync(_cuarrays.edge_gauss_weight, double_pitch, grid.triangle_infos.edge_gauss_weight, gsize, gsize, EDGE_GPOINTS*TRIANGLE_EDGES, cudaMemcpyHostToDevice);

	cudaMemcpyAsync(_cuarrays.airfoil, grid.airfoil, grid.wall_elem.size()*sizeof(int), cudaMemcpyHostToDevice);

	if ( cudaPeekAtLastError()!=cudaSuccess )
	{
		throw CMyException(cudaGetErrorString(cudaPeekAtLastError()));
	}
}

void CCUDARkdgSolver::initRKDG()
{
	int num = grid.getCellNumber();
	double  ut = sqrt(gamma*pref/rhoref)*mach;

	double u = ut * cos(alpha);
	double v = ut * sin(alpha);

	// �����ڴ�
	_freedom_rho  = new double[num*BASIS_FUNCTIONS];
	_freedom_rhou = new double[num*BASIS_FUNCTIONS];
	_freedom_rhov = new double[num*BASIS_FUNCTIONS];
	_freedom_rhoE = new double[num*BASIS_FUNCTIONS];

	cudaHostAlloc((void**)&_dt, sizeof(double), cudaHostAllocDefault);
	cudaHostAlloc((void**)&_residual, sizeof(double)*RESIDUAL_VARS, cudaHostAllocDefault);

	if ( cudaPeekAtLastError()!=cudaSuccess )
		throw CMyException(cudaGetErrorString(cudaPeekAtLastError()));

	// ��ʼ�����ɶȵ�ֵ
	for ( int i=0; i<num; ++i )
	{
		_freedom_rho[i] = rhoref;
		_freedom_rhou[i] = rhoref *u;
		_freedom_rhov[i] = rhoref*v;
		_freedom_rhoE[i] = rhoref*(ut*ut)/2 + pref/(rhoref*(gamma-1));
	}

	int dev_pitch = _cuarrays.getDoublePitch();
	int host_pitch = sizeof(double)*num;
	cudaMemsetAsync(_cuarrays.freedom_rho,  0, dev_pitch*BASIS_FUNCTIONS);
	cudaMemsetAsync(_cuarrays.freedom_rhou, 0, dev_pitch*BASIS_FUNCTIONS);
	cudaMemsetAsync(_cuarrays.freedom_rhov, 0, dev_pitch*BASIS_FUNCTIONS);
	cudaMemsetAsync(_cuarrays.freedom_rhoE, 0, dev_pitch*BASIS_FUNCTIONS);

	cudaMemcpy2DAsync(_cuarrays.freedom_rho,  dev_pitch, _freedom_rho,  host_pitch, host_pitch, 1, cudaMemcpyHostToDevice);														   
	cudaMemcpy2DAsync(_cuarrays.freedom_rhou, dev_pitch, _freedom_rhou, host_pitch, host_pitch, 1, cudaMemcpyHostToDevice);														   
	cudaMemcpy2DAsync(_cuarrays.freedom_rhov, dev_pitch, _freedom_rhov, host_pitch, host_pitch, 1, cudaMemcpyHostToDevice);														   
	cudaMemcpy2DAsync(_cuarrays.freedom_rhoE, dev_pitch, _freedom_rhoE, host_pitch, host_pitch, 1, cudaMemcpyHostToDevice);
}

void CCUDARkdgSolver::getTimeStep(int tnum)
{
	kernel_getTimeStep<<<1,reduction_threads, sizeof(double)*reduction_threads>>>(
		tnum, gamma, cfl, _cuarrays.ddt,

		_cuarrays.freedom_rho,  _cuarrays.freedom_rhou,
		_cuarrays.freedom_rhov, _cuarrays.freedom_rhoE,
		
		_cuarrays.perimeter, _cuarrays.area
		);
}

void CCUDARkdgSolver::calculateConVars(int tnum, int double_pitch, int blocks)
{
	size_t size = sizeof(double)*threads_per_block*CONSERVATIVE_VARS;

	kernel_calculateConVars<<<blocks,threads_per_block, size>>>(
		tnum, double_pitch,
		_cuarrays.freedom_rho,   _cuarrays.freedom_rhou,
		_cuarrays.freedom_rhov,  _cuarrays.freedom_rhoE,

		_cuarrays.convar_rho_vol,   _cuarrays.convar_rhou_vol,
		_cuarrays.convar_rhov_vol,  _cuarrays.convar_rhoE_vol,

		_cuarrays.convar_rho_edge,  _cuarrays.convar_rhou_edge, 
		_cuarrays.convar_rhov_edge, _cuarrays.convar_rhoE_edge, 

		_cuarrays.vol_bf_value,  _cuarrays.edge_bf_value
		);
}

void CCUDARkdgSolver::boundaryCondition(int tnum, int num, int double_pitch, double rho, double rhou, double rhov, double rhoE)
{
	// �߽絥Ԫ����Ŀ���٣�Ϊ��������ܣ���ÿ��block����ɼ����߳������������
	int threads = 64;

	int blocks = ((num-tnum)%threads) ? (num-tnum)/threads+1 : (num-tnum)/threads;

	kernel_boundaryCondition<<<blocks,threads>>>(
		tnum, num, double_pitch,
		rho, rhou, rhov, rhoE,

		_cuarrays.convar_rho_edge,   _cuarrays.convar_rhou_edge,
		_cuarrays.convar_rhov_edge,  _cuarrays.convar_rhoE_edge,

		_cuarrays.freedom_rho,   _cuarrays.freedom_rhou,
		_cuarrays.freedom_rhov,  _cuarrays.freedom_rhoE,

		_cuarrays.neighbour,     _cuarrays.sharedEdge,
		  
		_cuarrays.triangle_flag, _cuarrays.outer_normal_vector
		);
}

void CCUDARkdgSolver::calculateVolumeRHS(int tnum, int double_pitch, int blocks)
{
	size_t size = sizeof(double)*threads_per_block*VOLUME_GPOINTS;
	
	kernel_calculateVolumeRHS<<<blocks, threads_per_block, size>>>(
		tnum, double_pitch, gamma, 

		_cuarrays.convar_rho_vol, _cuarrays.convar_rhou_vol, 
		_cuarrays.convar_rhov_vol, _cuarrays.convar_rhoE_vol,

		_cuarrays.rhs_volume_rho, _cuarrays.rhs_volume_rhou, 
		_cuarrays.rhs_volume_rhov, _cuarrays.rhs_volume_rhoE,

		_cuarrays.vol_gauss_weight, _cuarrays.vol_bdf_value
		);
}

void CCUDARkdgSolver::calculateLFCoeff(int tnum, int ipitch_num, int dpitch_num, int blocks)
{

	kernel_calculateLFCoeff<<<blocks, threads_per_block>>>(
		tnum, ipitch_num, dpitch_num, gamma, 
		_cuarrays.outer_normal_vector,  _cuarrays.neighbour, 

		_cuarrays.freedom_rho, _cuarrays.freedom_rhou, 
		_cuarrays.freedom_rhov, _cuarrays.freedom_rhoE, 

		_cuarrays.lfflux_coeff
		);
}

void CCUDARkdgSolver::calculateEdgeFG(int tnum, int num, int double_pitch, int blocks)
{
	// �˴���Ҫ����ĵ�Ԫ������������һ�����Ӷ��߳̿���Ҫ���¶���
	blocks = (num%threads_per_block) ? num/threads_per_block+1 : num/threads_per_block;
	
	kernel_calculateEdgeFG<<<blocks, threads_per_block>>>(
		tnum, num, double_pitch, gamma, 

		_cuarrays.convar_rho_edge,  _cuarrays.convar_rhou_edge,
		_cuarrays.convar_rhov_edge, _cuarrays.convar_rhoE_edge,

		_cuarrays.fedge_rho,  _cuarrays.fedge_rhou,
		_cuarrays.fedge_rhov, _cuarrays.fedge_rhoE,

		_cuarrays.gedge_rho,  _cuarrays.gedge_rhou,
		_cuarrays.gedge_rhov, _cuarrays.gedge_rhoE
		);
}

void CCUDARkdgSolver::calculateFlux(int tnum, int int_pitch, int double_pitch, int blocks)
{
/*
	kernel_calculateFlux<<<blocks, threads_per_block>>>(
		tnum, int_pitch, double_pitch,
		_cuarrays.neighbour, _cuarrays.sharedEdge,

		_cuarrays.convar_rho_edge,  _cuarrays.convar_rhou_edge,
		_cuarrays.convar_rhov_edge, _cuarrays.convar_rhoE_edge,

		_cuarrays.fedge_rho,  _cuarrays.fedge_rhou,
		_cuarrays.fedge_rhov, _cuarrays.fedge_rhoE,

		_cuarrays.gedge_rho,  _cuarrays.gedge_rhou,
		_cuarrays.gedge_rhov, _cuarrays.gedge_rhoE,

		_cuarrays.outer_normal_vector, _cuarrays.lfflux_coeff,

		_cuarrays.lfflux_rho,  _cuarrays.lfflux_rhou,
		_cuarrays.lfflux_rhov, _cuarrays.lfflux_rhoE
		);
*/
	kernel_calculateFlux<<<blocks, threads_per_block>>>(
		tnum, int_pitch, double_pitch,
		_cuarrays.neighbour, _cuarrays.sharedEdge,

		_cuarrays.convar_rho_edge,  _cuarrays.convar_rhou_edge,
//		_cuarrays.convar_rhov_edge, _cuarrays.convar_rhoE_edge,

		_cuarrays.fedge_rho,  _cuarrays.fedge_rhou,
//		_cuarrays.fedge_rhov, _cuarrays.fedge_rhoE,

		_cuarrays.gedge_rho,  _cuarrays.gedge_rhou,
//		_cuarrays.gedge_rhov, _cuarrays.gedge_rhoE,

		_cuarrays.outer_normal_vector, _cuarrays.lfflux_coeff,

		_cuarrays.lfflux_rho,  _cuarrays.lfflux_rhou
//		_cuarrays.lfflux_rhov, _cuarrays.lfflux_rhoE
		);
		
	kernel_calculateFlux<<<blocks, threads_per_block>>>(
		tnum, int_pitch, double_pitch,
		_cuarrays.neighbour, _cuarrays.sharedEdge,

//		_cuarrays.convar_rho_edge,  _cuarrays.convar_rhou_edge,
		_cuarrays.convar_rhov_edge, _cuarrays.convar_rhoE_edge,

//		_cuarrays.fedge_rho,  _cuarrays.fedge_rhou,
		_cuarrays.fedge_rhov, _cuarrays.fedge_rhoE,

//		_cuarrays.gedge_rho,  _cuarrays.gedge_rhou,
		_cuarrays.gedge_rhov, _cuarrays.gedge_rhoE,

		_cuarrays.outer_normal_vector, _cuarrays.lfflux_coeff,

//		_cuarrays.lfflux_rho,  _cuarrays.lfflux_rhou,
		_cuarrays.lfflux_rhov, _cuarrays.lfflux_rhoE
		);
}

void CCUDARkdgSolver::calculateEdgeRHS(int tnum, int double_pitch, int blocks)
{
	size_t size = sizeof(double)*threads_per_block*TRIANGLE_EDGES*EDGE_GPOINTS;

	kernel_calculateEdgeRHS<<<blocks, threads_per_block, size>>>(
		tnum, double_pitch,
		_cuarrays.edge_gauss_weight, _cuarrays.edge_bf_value,

		_cuarrays.lfflux_rho,  _cuarrays.lfflux_rhou, 
		_cuarrays.lfflux_rhov, _cuarrays.lfflux_rhoE, 

		_cuarrays.rhs_edge_rho,  _cuarrays.rhs_edge_rhou, 
		_cuarrays.rhs_edge_rhov, _cuarrays.rhs_edge_rhoE,

		_cuarrays.rhs_volume_rho,  _cuarrays.rhs_volume_rhou,
		_cuarrays.rhs_volume_rhov, _cuarrays.rhs_volume_rhoE
		);

}

void CCUDARkdgSolver::rkdgStepOne(double dt, int tnum, int double_pitch, int blocks)
{

	kernel_rkdgStepOne<<<blocks, threads_per_block>>>(
		tnum, double_pitch, dt, _cuarrays.mass_coeff,

		_cuarrays.freedom_rho,  _cuarrays.freedom_rhou, 
		_cuarrays.freedom_rhov, _cuarrays.freedom_rhoE, 

		_cuarrays.rhs_edge_rho,  _cuarrays.rhs_edge_rhou,
		_cuarrays.rhs_edge_rhov, _cuarrays.rhs_edge_rhoE
		);
}

void CCUDARkdgSolver::rkdgStepTwo(double dt, int tnum, int double_pitch, int blocks)
{

	kernel_rkdgStepTwo<<<blocks, threads_per_block>>>(
		tnum, double_pitch, dt, _cuarrays.mass_coeff,
		_cuarrays.freedom_rho,  _cuarrays.freedom_rhou, 
		_cuarrays.freedom_rhov, _cuarrays.freedom_rhoE, 

		_cuarrays.rhs_edge_rho,  _cuarrays.rhs_edge_rhou,
		_cuarrays.rhs_edge_rhov, _cuarrays.rhs_edge_rhoE,

		_cuarrays.freedom_rho_old,  _cuarrays.freedom_rhou_old, 
		_cuarrays.freedom_rhov_old, _cuarrays.freedom_rhoE_old
		);
}

void CCUDARkdgSolver::rkdgStepThree(double dt, int tnum, int double_pitch, int blocks)
{

	kernel_rkdgStepThree<<<blocks, threads_per_block>>>(
		tnum, double_pitch, dt, _cuarrays.mass_coeff,

		_cuarrays.freedom_rho,     _cuarrays.freedom_rhou, 
		_cuarrays.freedom_rhov,    _cuarrays.freedom_rhoE, 

		_cuarrays.rhs_edge_rho,    _cuarrays.rhs_edge_rhou,
		_cuarrays.rhs_edge_rhov,   _cuarrays.rhs_edge_rhoE,

		_cuarrays.freedom_rho_old,     _cuarrays.freedom_rhou_old, 
		_cuarrays.freedom_rhov_old,    _cuarrays.freedom_rhoE_old
		);
}

void CCUDARkdgSolver::calculateResidual(int tnum)
{
	kernel_calculateResidual<<<1,reduction_threads, sizeof(double)*reduction_threads*RESIDUAL_VARS>>>(
		tnum, 
		_cuarrays.freedom_rho, _cuarrays.freedom_rhoE,

		_cuarrays.freedom_rho_old, _cuarrays.freedom_rhoE_old,

		_cuarrays.residual
		);
}

void CCUDARkdgSolver::rkdgAdvance(void)
{
	ofstream fout;
	if ( log_history=='Y' )
	{
		fout.open(residual_file.c_str());
		if ( !fout )
			throw CMyException("Failed to open residual log file: "+residual_file);

		fout<<"N, rho"<<endl;
	}
	
	double resi(0);
	double nt(0);
	int count(0);

	int tnum = grid.getTriangleNumber();
	int num  = grid.getCellNumber();

	int blocks = (tnum%threads_per_block) ? tnum/threads_per_block+1 : tnum/threads_per_block;

	double ut   = sqrt(gamma*pref/rhoref)*mach;
	double rhou = rhoref*ut*cos(alpha);
	double rhov = rhoref*ut*sin(alpha);
	double rhoE = 0.5*rhoref*(ut*ut) + pref/(gamma-1);

	bool copy(false);

	cudaError_t error;
	size_t pitch = _cuarrays.getDoublePitch();
	int pitch_num = pitch / sizeof(double);
	int ipitch_num = _cuarrays.getIntPitch() / sizeof(int);

	cudaEvent_t time_start, time_stop;

	cudaEventCreateWithFlags(&time_start, cudaEventDisableTiming|cudaEventBlockingSync);
	cudaEventCreateWithFlags(&time_stop,  cudaEventDisableTiming|cudaEventBlockingSync);

	if ( log_history=='Y' )
		copy = true;

	// ȷ��֮ǰCUDA�ĳ�ʼ���������Ѿ����
	cudaDeviceSynchronize();

	do 
	{
		++ count;

		cudaEventRecord(time_start);
		// ���㵱ǰʱ�䲽��
		getTimeStep(tnum);

		cudaEventRecord(time_stop);

		// ��������ɶ�
		cudaMemcpy2DAsync(_cuarrays.freedom_rho_old,  pitch, _cuarrays.freedom_rho,  pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(_cuarrays.freedom_rhou_old, pitch, _cuarrays.freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(_cuarrays.freedom_rhov_old, pitch, _cuarrays.freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							 
		cudaMemcpy2DAsync(_cuarrays.freedom_rhoE_old, pitch, _cuarrays.freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);
			
		for ( int i=0; i<RUNGE_KUTTA_STEPS; ++i )
		{
			// �����غ�����ֵ
			calculateConVars(tnum, pitch_num, blocks);

			// ����߽�����
			boundaryCondition(tnum, num, pitch_num, rhoref, rhou, rhov, rhoE);

			// ��������ֲв�
			calculateVolumeRHS(tnum, pitch_num, blocks);

			// ����LFͨ��ϵ��
			calculateLFCoeff(tnum, ipitch_num, pitch_num, blocks);

			// ����f, g�ڱ��ϵ�ֵ
			calculateEdgeFG(tnum, num, pitch_num, blocks);

			calculateFlux(tnum, ipitch_num, pitch_num, blocks);

			// �����߻��ֲв�
			calculateEdgeRHS(tnum, pitch_num, blocks);

			// ʱ���ƽ�
			switch (i)
			{
			case 0:
				cudaEventSynchronize(time_stop);

				// ��ʱ�䲽�����͵�����
				cudaMemcpy(_dt, _cuarrays.ddt, sizeof(double), cudaMemcpyDeviceToHost);

				if ( 0==(count-1)%2000 )
					cout<<"Step: "<<count<<", time step: "<<_dt[0]<<endl;

				if ( (_terminal_time-nt)<_dt[0] )
				{
					_dt[0] = _terminal_time -  nt;
				}

				// ʱ�䲽�ƽ�
				rkdgStepOne(_dt[0], tnum, pitch_num, blocks);

				break;

			case 1:
				rkdgStepTwo(_dt[0], tnum, pitch_num, blocks);
				break;

			case 2:
				rkdgStepThree(_dt[0], tnum, pitch_num, blocks);
				break;

			default:
				throw CMyException("impossible case!");
				break;
			}
		}

		
		if ( copy && (count-1) )
		{
			// ���Ʋв�����
			cudaMemcpy(_residual, _cuarrays.residual,
						sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);

			if ( 0==(count-1)%2000 )
				cout<<"Current time: "<<nt<<"   rhomax: "<<_residual[0]/rhoref<<"   E: "<<_residual[1]/rhoE<<endl;

			fout<<count<<"   "<<log(_residual[0]/rhoref)/log(10.0)<<endl;
			
			resi = log(_residual[0]/rhoref)/log(10.0);
		}

		// ����в�
		calculateResidual(tnum);

		
		// ��ʱ�ƽ�
		nt += _dt[0];

		error = cudaPeekAtLastError();
		if ( error!=cudaSuccess )
			throw CMyException(cudaGetErrorString(error));

	} while ( resi > -6.5 );
	cudaDeviceSynchronize();
	
	if ( copy )
	{
		// ���Ʋв�����
		cudaMemcpy(_residual, _cuarrays.residual,
			sizeof(double)*RESIDUAL_VARS, cudaMemcpyDeviceToHost);

		if ( 0==(count-1)%print_interval )
			cout<<"��ǰʱ�䣺 "<<nt-_dt[0]<<"   rhomax: "<<_residual[0]/rhoref<<"   E: "<<_residual[1]/rhoE<<endl;

		fout<<count<<"   "<<log(_residual[0]/rhoref)/log(10.0)<<endl;
	}


	cudaEventDestroy(time_start);
	cudaEventDestroy(time_stop);

	if ( log_history=='Y' )
		fout.close();

}

void CCUDARkdgSolver::outputSolution(string solfile)
{
	ofstream fos("rho.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < grid.getCellNumber(); j++) {
			fos<<_freedom_rho[j+i*grid.getCellNumber()]<<endl;
		}
	}
	fos.close();

	fos.open("rhou.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < grid.getCellNumber(); j++) {
			fos<<_freedom_rhou[j+i*grid.getCellNumber()]<<endl;
		}
	}
	fos.close();

	fos.open("rhov.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < grid.getCellNumber(); j++) {
			fos<<_freedom_rhov[j+i*grid.getCellNumber()]<<endl;
		}
	}
	fos.close();

	fos.open("rhoE.dat");
	for (int i = 0; i < 6; i++) {
		for (int j = 0; j < grid.getCellNumber(); j++) {
			fos<<_freedom_rhoE[j+i*grid.getCellNumber()]<<endl;
		}
	}
	fos.close();


	ofstream fout(solfile.c_str());
	if ( !fout )
	{
		cout<<"Failed to open solution file: "<<solution_file<<" and output will be omitted."<<endl;
		return;
	}

	int i;
	int vnum, tnum;
	double rho, u, v, rhoE, p, a, ma;
	vnum = grid.getVerticeNumber();
	tnum = grid.getTriangleNumber();

	fout<<"TITLE=RKDG"<<endl;
	fout<<"VARIABLES=X , Y , rho , u , v , p, Ma , FLAG"<<endl;
	fout<<"ZONE T= T1 N= "<<vnum<<" , E= "<<tnum<<" , ZONETYPE=FETRIANGLE"<<endl;
	fout<<"DATAPACKING=BLOCK"<<endl;
	fout<<"VARLOCATION=([1-2]=NODAL,[3-8]=CELLCENTERED)"<<endl;
	fout<<"DT=(SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE SINGLE)"<<endl;

	for ( i=0; i<vnum; ++i )
	{
		fout<<grid.vertice[i].getX()<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}

	for ( i=0; i<vnum; ++i )
	{
		fout<<grid.vertice[i].getY()<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	for ( i=0; i<tnum; ++i )
	{
		fout<<_freedom_rho[i]<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	for ( i=0; i<tnum; ++i )
	{
		fout<<_freedom_rhou[i]/_freedom_rho[i]<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	for ( i=0; i<tnum; ++i )
	{
		fout<<_freedom_rhov[i]/_freedom_rho[i]<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	for ( i=0; i<tnum; ++i )
	{
		rho = _freedom_rho[i];
		u	= _freedom_rhou[i]/rho;
		v	= _freedom_rhov[i]/rho;
		rhoE = _freedom_rhoE[i];

		p = (gamma-1)*(rhoE-0.5*rho*(u*u+v*v));

		fout<<p<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	for ( i=0; i<tnum; ++i )
	{
		rho = _freedom_rho[i];
		u	= _freedom_rhou[i]/rho;
		v	= _freedom_rhov[i]/rho;
		rhoE = _freedom_rhoE[i];

		p = (gamma-1)*(rhoE-0.5*rho*(u*u+v*v));
		a = sqrt(gamma*p/rho);
		ma = sqrt(u*u+v*v)/a;

		fout<<ma<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	// ���������
	for ( i=0; i<tnum; ++i )
	{
		fout<<"1"<<"  ";
		if ( i%6==0 )
		{
			fout<<endl;
		}
	}
	fout<<endl;

	for ( i=0; i<tnum; ++i )
	{
		fout<<grid.tri_vertice[3*i]+1<<"    "<<grid.tri_vertice[3*i+1]+1<<"    "<<grid.tri_vertice[3*i+2]+1<<endl;
	}
	fout.close();
}



void CCUDARkdgSolver::printConfig( ostream& out )
{
	if ( !out )
	{
		cerr<<"Invalid output stream and output will be omitted."<<endl;
		return;
	}

	// �����������
	out<<"===================="<<endl;
	out<<"Program configures: "<<endl;
	out<<"===================="<<endl;
	out<<"gamma:              "<<gamma<<endl;
	out<<"alpha:              "<<alpha*180/(4*atan(1.0))<<endl;
	out<<"mach:               "<<mach<<endl;
	out<<"cfl:                "<<cfl<<endl;
	out<<"rhoref:             "<<rhoref<<endl;
	out<<"pref:               "<<pref<<endl;
	out<<"time: 			  "<<_terminal_time<<endl;
	out<<"===================="<<endl;
	out<<"gridconf:           "<<gridconf<<endl;
	out<<"solution:           "<<solution_file<<endl;
	out<<"residualfile:       "<<residual_file<<endl;
	out<<"printinterval:      "<<print_interval<<endl;
	out<<"loghistory:         "<<log_history<<endl;
	out<<"===================="<<endl;
	out<<"threads_per_block:  "<<threads_per_block<<endl;
	out<<"reduction_threads:  "<<reduction_threads<<endl;
	out<<"===================="<<endl<<endl;
}


CCUDARkdgSolver::~CCUDARkdgSolver()
{
	delete []_freedom_rho;
	delete []_freedom_rhou;
	delete []_freedom_rhov;
	delete []_freedom_rhoE;

	cudaFreeHost(_residual);
	cudaFreeHost(_dt);
}