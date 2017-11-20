#include "cudarkdgsolvertester.h"

CCUDARkdgSolverTester::CCUDARkdgSolverTester()
{

}

void CCUDARkdgSolverTester::initEnvironment()
{
	// 初始化程序配置并输出程序配置
	initConfig();

	// 初始化网格
	grid.config_file = gridconf;
	grid.initializeGrid();

	// 初始化网格上基函数等信息
	grid.triangle_infos.allocateMemory(grid.getCellNumber());
	grid.initializeTriangleInfos();

	// 标记网格单元
	grid.markBoundaryTriangles();

	// 分配GPU内存
	_cuarrays.allocateMemory(grid.getCellNumber());

	// 初始化GPU上数组的值
	//_cuarrays.memsetArrays(grid.getCellNumber());

	// 将三角单元信息传送到GPU
	copyTriangleInfosToGPU();

	// 初始化RKDG自由度，并将初始化数据传到GPU
	initRKDG();

}

void CCUDARkdgSolverTester::runTest()
{
	double nt(0);
	int count(0);

	int tnum = grid.getTriangleNumber();
	int num  = grid.getCellNumber();

	int blocks = (tnum%threads_per_block) ? tnum/threads_per_block+1 : tnum/threads_per_block;

	double ut = sqrt(gamma*pref/rhoref)*mach;
	double rhou = rhoref*ut*cos(alpha);
	double rhov = rhoref*ut*sin(alpha);
	double rhoE  = 0.5*rhoref*(ut*ut) + pref/(gamma-1);

	cudaError_t error;
	size_t pitch = _cuarrays.getDoublePitch();
	size_t host_pitch = sizeof(double)*num;
	int pitch_num = pitch / sizeof(double);
	int ipitch_num = _cuarrays.getIntPitch() / sizeof(int);

	int i;

	// 确保之前CUDA的初始化工作都已经完成
	cudaDeviceSynchronize();

	cout.width(15);
	cout.precision(12);

	do 
	{
		++ count;

		// 计算当前时间步长
		getTimeStep(tnum);

		cudaDeviceSynchronize();
		error = cudaPeekAtLastError();
		if ( error!=cudaSuccess )
			throw CMyException(cudaGetErrorString(error));

		// 将时间步长传送到本地
		cudaMemcpy(_dt, _cuarrays.ddt, sizeof(double), cudaMemcpyDeviceToHost);
		cout<<"Time step: "<<_dt<<endl;

		// 保存旧自由度
		cudaMemcpy2DAsync(_cuarrays.freedom_rho_old,  pitch, _cuarrays.freedom_rho,  pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(_cuarrays.freedom_rhou_old, pitch, _cuarrays.freedom_rhou, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							  
		cudaMemcpy2DAsync(_cuarrays.freedom_rhov_old, pitch, _cuarrays.freedom_rhov, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);							 
		cudaMemcpy2DAsync(_cuarrays.freedom_rhoE_old, pitch, _cuarrays.freedom_rhoE, pitch, pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToDevice);

		for ( i=0; i<RUNGE_KUTTA_STEPS; ++i )
		{
			// 计算守恒量的值
			calculateConVars(tnum, pitch_num, blocks);

			cudaMemcpy2DAsync(_freedom_rho, host_pitch, _cuarrays.convar_rho_vol,   pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.convar_rhou_vol, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.convar_rhov_vol, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.convar_rhoE_vol, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/convarVolRho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/convarVolRhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/convarVolRhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/convarVolRhoE.dat");

			cudaMemcpy2DAsync(_freedom_rho, host_pitch, _cuarrays.convar_rho_edge,   pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.convar_rhou_edge, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.convar_rhov_edge, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.convar_rhoE_edge, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/convarEdgerho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/convarEdgerhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/convarEdgerhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/convarEdgerhoE.dat");



			// 处理边界条件
			boundaryCondition(tnum, num, pitch_num, rhoref, rhou, rhov, rhoE);

			cudaMemcpy2DAsync(_freedom_rho, host_pitch, _cuarrays.freedom_rho,   pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.freedom_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.freedom_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.freedom_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/bndFreedomRho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/bndFreedomRhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/bndFreedomRhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/bndFreedomRhoE.dat");

			cudaMemcpy2DAsync(_freedom_rho, host_pitch, _cuarrays.convar_rho_edge,   pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.convar_rhou_edge, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.convar_rhov_edge, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.convar_rhoE_edge, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/convarEdgerho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/convarEdgerhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/convarEdgerhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/convarEdgerhoE.dat");



			// 计算体积分残差
			calculateVolumeRHS(tnum, pitch_num, blocks);

			cudaMemcpy2DAsync(_freedom_rho, host_pitch, _cuarrays.rhs_volume_rho,   pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.rhs_volume_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.rhs_volume_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.rhs_volume_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/RHSVolRho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/RHSVolRhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/RHSVolRhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/RHSVolRhoE.dat");


			// 计算LF通量系数
			calculateLFCoeff(tnum, ipitch_num, pitch_num, blocks);
			cudaMemcpy2DAsync(_freedom_rho, host_pitch, _cuarrays.lfflux_coeff,   pitch, host_pitch, TRIANGLE_EDGES, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, TRIANGLE_EDGES, _freedom_rho,  "output/lffluxCoeff.dat");




			// 计算f, g在边上的值
			calculateEdgeFG(tnum, num, pitch_num, blocks);
			cudaMemcpy2DAsync(_freedom_rho,  host_pitch, _cuarrays.fedge_rho,  pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.fedge_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.fedge_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.fedge_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/edgeFrho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/edgeFrhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/edgeFrhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/edgeFrhoE.dat");

			cudaMemcpy2DAsync(_freedom_rho,  host_pitch, _cuarrays.gedge_rho,  pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.gedge_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.gedge_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.gedge_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/edgeGrho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/edgeGrhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/edgeGrhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/edgeGrhoE.dat");




			calculateFlux(tnum, ipitch_num, pitch_num, blocks);

			cudaMemcpy2DAsync(_freedom_rho,  host_pitch, _cuarrays.lfflux_rho,  pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.lfflux_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.lfflux_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
			cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.lfflux_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

			cudaDeviceSynchronize();
			error = cudaPeekAtLastError();
			if ( error!=cudaSuccess )
				throw CMyException(cudaGetErrorString(error));

			output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/lfFluxrho.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/lfFluxrhou.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/lfFluxrhov.dat");
			output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/lfFluxrhoE.dat");


			// 计算线积分残差
			calculateEdgeRHS(tnum, pitch_num, blocks);

			// 时间推进
			switch (i)
			{
			case 0:
				// 时间步推进
				rkdgStepOne(_dt[0], tnum, pitch_num, blocks);

				cudaMemcpy2DAsync(_freedom_rho,  host_pitch, _cuarrays.freedom_rho,  pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.freedom_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.freedom_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.freedom_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

				cudaDeviceSynchronize();
				error = cudaPeekAtLastError();
				if ( error!=cudaSuccess )
					throw CMyException(cudaGetErrorString(error));

				output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/rkdgT1rho.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/rkdgT1rhou.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/rkdgT1rhov.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/rkdgT1rhoE.dat");

				break;

			case 1:
				rkdgStepTwo(_dt[0], tnum, pitch_num, blocks);

				cudaMemcpy2DAsync(_freedom_rho,  host_pitch, _cuarrays.freedom_rho,  pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.freedom_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.freedom_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.freedom_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

				cudaDeviceSynchronize();
				error = cudaPeekAtLastError();
				if ( error!=cudaSuccess )
					throw CMyException(cudaGetErrorString(error));

				output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/rkdgT2rho.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/rkdgT2rhou.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/rkdgT2rhov.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/rkdgT2rhoE.dat");

				break;

			case 2:
				rkdgStepThree(_dt[0], tnum, pitch_num, blocks);

				cudaMemcpy2DAsync(_freedom_rho,  host_pitch, _cuarrays.freedom_rho,  pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhou, host_pitch, _cuarrays.freedom_rhou, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhov, host_pitch, _cuarrays.freedom_rhov, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);
				cudaMemcpy2DAsync(_freedom_rhoE, host_pitch, _cuarrays.freedom_rhoE, pitch, host_pitch, BASIS_FUNCTIONS, cudaMemcpyDeviceToHost);

				cudaDeviceSynchronize();
				error = cudaPeekAtLastError();
				if ( error!=cudaSuccess )
					throw CMyException(cudaGetErrorString(error));

				output2DVar(num, BASIS_FUNCTIONS, _freedom_rho,  "output/rkdgT3rho.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhou, "output/rkdgT3rhou.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhov, "output/rkdgT3rhov.dat");
				output2DVar(num, BASIS_FUNCTIONS, _freedom_rhoE, "output/rkdgT3rhoE.dat");


				break;

			default:
				throw CMyException("impossible case!");
				break;
			}
		}


		// 计算残差
		calculateResidual(tnum);


		// 计时推进
		nt += _dt[0];

		error = cudaPeekAtLastError();
		if ( error!=cudaSuccess )
			throw CMyException(cudaGetErrorString(error));

	} while ( nt<_terminal_time );
}

void CCUDARkdgSolverTester::output2DVar(int rnum, int cnum, const double *var, const string &output)
{

	ofstream fout(output.c_str());
	if (!fout)
	{
		cout<<"文件： "<<output<<" 打开失败，文件输出过程将被忽略"<<endl;
		return;
	}

	fout.width(15);
	fout.precision(12);
	for ( int i=0; i<rnum; ++i )
	{
		for ( int j=0; j<cnum; ++j )
		{
			fout<<(fabs(var[i+j*rnum])>1e-12? var[i+j*rnum]:0)<<"   ";
		}
		fout<<endl;
	}

	fout.close();
}


CCUDARkdgSolverTester::~CCUDARkdgSolverTester()
{

}



