/**
 * 测试函数的入口文件
 * 该文件根据需要进行单元函数测试
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
 */
 
#include "../inc/cppstdheaders.h"
#include "cudarkdgsolvertester.h"

using namespace std;

int testmain(void)
{
	CCUDARkdgSolverTester solver;
	
	try {
		solver.config_file = "input/main.conf";

		solver.initEnvironment();
		
//		solver.grid.testArea();
//		solver.grid.testEdgeBfValue();
//		solver.grid.testEdgeGaussWeight();
//		solver.grid.testMassCoeff();
//		solver.grid.testNeighbour();
//		solver.grid.testOuterNormalVector();
//		solver.grid.testPerimeter();
//		solver.grid.testSharedEdge();
//		solver.grid.testTriangleFlag();
//		solver.grid.testVolumeBdfValue();
//		solver.grid.testVolumeBfValue();
//		solver.grid.testVolumeGaussWeight();
		// 网格信息一致 tlanyan 2013-5-10 17:14
		


//		solver.testConVars();

//		solver.testVolumeRHS();

// 		solver.testLFCoeff();

//		solver.testEdgeFG();

//		solver.testFlux();

//		solver.testRKDGStepOne();

//		solver.testRKDGStepTwo();

		solver.testRKDGStepThree();

	}
	catch ( const exception& e )
   	{
		cout<<endl<<"程序运行中发生错误，错误描述: "<<endl;
		cout<<e.what()<<endl;
		cout<<endl<<"请根据错误修改程序，按回车键结束程序"<<endl;
		getchar();
		
		exit(-1);
	}
	
	cout<<"测试结束，按回车键退出程序"<<endl;
	getchar();
	return 0;
}