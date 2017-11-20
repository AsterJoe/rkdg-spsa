/**
 * 在GPU上使用RKDG方法计算二维非结构网格上NACA0012机翼扰流，特点如下：
 * 1). 在每个三角单元上使用二次正交多项式作为基函数
 * 2). 三阶龙格库塔法进行时间推进
 *
 * 感谢： cjbuaa 提供算法和参考源程序
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention 你可以无限制的使用本文件（开源的或商业的目的），唯一的要求是保留作者信息、版权说明以及该使用注意。
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