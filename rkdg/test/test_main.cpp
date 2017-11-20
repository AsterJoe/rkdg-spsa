/**
 * ���Ժ���������ļ�
 * ���ļ�������Ҫ���е�Ԫ��������
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
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
		// ������Ϣһ�� tlanyan 2013-5-10 17:14
		


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
		cout<<endl<<"���������з������󣬴�������: "<<endl;
		cout<<e.what()<<endl;
		cout<<endl<<"����ݴ����޸ĳ��򣬰��س�����������"<<endl;
		getchar();
		
		exit(-1);
	}
	
	cout<<"���Խ��������س����˳�����"<<endl;
	getchar();
	return 0;
}