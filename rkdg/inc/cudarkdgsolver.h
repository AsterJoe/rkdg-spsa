/**
 * CUDAʵ�ֵĶ�ά�ǽṹ������RKDG�㷨�ⷨ��
 * �ⷨ�����������������̣���ȡ�����ļ�����������ӿڶ����������񣬽���RKDG����ʱ���ƽ����Լ�������
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 
#pragma once

#include "cppstdheaders.h"
#include "config.h"
#include "unstructuredgrid.h"
#include "cudaarrays.h"
#include "kernelfunctions.cuh"
#include "defines.h"
#include "myexception.h"
#include "mytime.h"
using namespace std;

/**
 * ��ά�ǽṹ�����ϻ���CUDA��RKDG�ⷨ��
 * �ⷨ�����������������̣���ȡ�����ļ�����������ӿڶ����������񣬽���RKDG����ʱ���ƽ���������
 * ���г�Ա�����п����׳���׼�쳣�����ڵ��õ�ʱ��׽��
 * @example �����÷�����
 *  CCUDARkdgSolver solver;
 *  solver.config_file = config_filename
 *  solver.run();
 *  ����config_filename�ǳ���������ļ������ڸ������ļ����밴�ո�ʽ��д��������
 */
class CCUDARkdgSolver {
	// ����
	public:
		string title;				/**< �������� */

		double alpha;				/**< ���͹��� */
		double gamma;				/**< ������ȱ� */
		double mach;				/**< ��������� */
		double cfl;					/**< cfl�� */
		double rhoref;				/**< ��ʼ�ο��ܶ� */
		double pref;				/**< ��ʼ�ο�ѹǿ */
		
		char log_history;			/**< �Ƿ��¼������ʷ */
		int print_interval;			/**< ��Ļ��Ϣ������ */
	
		string config_file;			/**< ���������ļ� */
		string gridconf;			/**< �������������ļ� */
		string solution_file;		/**< ����Ľ�����ļ� */
		string log_file;			/**< ������־�ļ� */
		string residual_file;		/**< �в��¼�ļ� */
		
		CUnstructuredGrid grid;		/**< ����ӿ� */
		
		int threads_per_block;		/**< ÿ��CUDA�����̵߳���Ŀ */
		int reduction_threads;		/**< �ھۺϺ����ͨ��ֻ��һ���߳̿飩���߳���Ŀ */

//	protected:
		double _terminal_time;		/**< �����ʱ�� */

		double *_freedom_rho;		/**< rho���ɶ� */
		double *_freedom_rhou;		/**< rhou���ɶ� */
		double *_freedom_rhov;		/**< rhov���ɶ� */
		double *_freedom_rhoE;		/**< rhoE���ɶ� */

		double *_dt;				/**< ʱ�䲽�� */
		double *_residual;			/**< �в� */

		CCUDAArrays _cuarrays;		/**< GPU�����ʹ�õ����� */
		
	private:

	// ����
	public:
		/** Ĭ�Ϲ��캯�� */
		CCUDARkdgSolver();
		
		/** ��ʼ���������� */
		void initConfig(void);

		/** ���г���������� */
		void run(void);

		void runNext(void);
		
		/** ���CUDA�豸 */
		void detectCUDADevice(void);

		/** ������������Ϣ��GPU */
		void copyTriangleInfosToGPU(void);

		/** ��ʼ��RKDG���ɶȲ������ݴ��͵�GPU�� */
		void initRKDG(void);

		/** rkdgʱ���ƽ� */
		void rkdgAdvance(void);

		/**
		 * ���ʱ�䲽��
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 */
		void getTimeStep(int tnum);

		/**
		 *����غ���,
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] double_pitch integer �����˫��������Ŀ��
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void calculateConVars(int tnum, int double_pitch, int blocks);

		/**
		 * ����߽�����
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] num integer ����Ԫ����
		 * @param[in] double_pitch integer ������˫����������
		 * @param[in] rho double Զ��ѹǿ
		 * @param[in] rhou double Զ������
		 * @param[in] rhov double Զ������
		 * @param[in] rhoE double Զ������
		 */
		void boundaryCondition(int tnum, int num, int double_pitch, double rho, double rhou, double rhov, double rhoE);

		/**
		 * ���������
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] double_pitch integer �����˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void calculateVolumeRHS(int tnum, int double_pitch, int blocks);

		/**
		 * ����LFͨ��ϵ��
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] int_pitch integer ����������������
		 * @param[in] double_pitch integer ������˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void calculateLFCoeff(int tnum, int int_pitch, int double_pitch, int blocks);

		/**
		 * ������ϵ�FG
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] num integer ����Ԫ����
		 * @param[in] double_pitch integer ������˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void calculateEdgeFG(int tnum, int num, int double_pitch, int blocks);

		/**
		 * ����ͨ��
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] int_pitch integer ���������������
		 * @param[in] double_pitch integer �����˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void calculateFlux(int tnum, int int_pitch, int double_pitch, int blocks);

		/**
		 * ������ϵĲв�
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] double_pitch integer �����˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void calculateEdgeRHS(int tnum, int double_pitch, int blocks);

		/**
		 * rkdg ʱ�䲽��һ���ƽ�
		 * @param[in] dt double ʱ�䲽��
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] double_pitch integer �����˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void rkdgStepOne(double dt, int tnum, int double_pitch, int blocks);

		/**
		 * rkdg ʱ�䲽��һ���ƽ�
		 * @param[in] dt double ʱ�䲽��
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] double_pitch integer �����˫����������
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void rkdgStepTwo(double dt, int tnum, int double_pitch, int blocks);

		/**
		 * rkdg ʱ�䲽��һ���ƽ�
		 * @param[in] dt double ʱ�䲽��
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 * @param[in] double_pitch integer �����˫��������Ŀ��
		 * @param[in] blocks integer �߳̿����Ŀ
		 */
		void rkdgStepThree(double dt, int tnum, int double_pitch, int blocks);

		/** 
		 * �����غ����в�
		 * @param[in] tnum integer ���ǵ�Ԫ��Ŀ
		 */
		void calculateResidual(int tnum);

		/** �������ɶȵ����� */
		void copyFreedomToHost(void);

		/** ����� */
		//void outputSolution(void);
		void outputSolution(string solfile);

		/**
		 * �������������Ϣ
		 * @param[in] out ostream �����
		 */
		void printConfig(ostream& out);
	
		/** �������� */
		~CCUDARkdgSolver();
	protected:
	
	private:
	
};