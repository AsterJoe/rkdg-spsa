/**
 * ��GPU��ʹ��RKDG���������ά�ǽṹ������NACA0012�����������ص����£�
 * 1). ��ÿ�����ǵ�Ԫ��ʹ�ö�����������ʽ��Ϊ������
 * 2). �����������������ʱ���ƽ�
 *
 * ��л�� cjbuaa �ṩ�㷨�Ͳο�Դ����
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */

#include<mpi.h>
#include "../inc/cppstdheaders.h"
#include "../inc/cudarkdgsolver.h"
#include "../inc/spsaoptimize.h"
using namespace std;

int main(int argc, char *argv[])
{
	int myid;
	int nprocs;

	/*MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	if (myid == 0) {*/
		CCUDARkdgSolver solver;
		try {
			solver.config_file = "input/main.conf";

			solver.run();
		}
		catch ( const exception& e )
		{
			cout<<endl<<"Error occured, description: "<<endl;
			cout<<e.what()<<endl;
			cout<<endl<<"Please modify the program according to the hint."<<endl;
			getchar();

			exit(-1);
		}
	//}
	CSpsaOptimize optimize(&solver);

	//MPI_Finalize();
	cout<<"complete solving flow, press ENTER to exit."<<endl;
	getchar();
	return 0;
}