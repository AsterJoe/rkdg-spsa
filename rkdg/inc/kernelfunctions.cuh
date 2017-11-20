/**
 * ����CUDA�ĺ˺���
 * ����CUDA�ĺ˺�������Ϊȫ�ֺ������Ӷ���Ҫ������ⲿ��������
 *
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 
#pragma once 
#include "cppstdheaders.h"
#include "cuda_runtime.h"
#include "defines.h"
using namespace std;

extern "C"{

	/**
	 * ��GPU�����ʱ�䲽��
	 * @param[in] tnum integer ʵ����������Ŀ
	 * @param[in] gamma double ������ȱ�
	 * @param[out] ddt double ʱ�䲽��
	 * @param[in] freedom_rho double[] rho���ɶ�
	 * @param[in] freedom_rhou double[] rhou���ɶ�
	 * @param[in] freedom_rhov double[] rhov���ɶ�
	 * @param[in] freedom_rhoE double[] rhoE���ɶ�
	 * @param[in] perimeter double[] �������ܳ�
	 * @param[in] area double[] ���������
	 */
	__global__ void kernel_getTimeStep(int tnum, double gamma, double cfl, double * ddt,

				const double *freedom_rho, const double *freedom_rhou, 
				const double *freedom_rhov, const double *freedom_rhoE, 

				const double *perimeter, const double *area
				);

	/**
	 * �����غ���
	 * @param[in] tnum integer ʵ����������Ŀ
	 * @param[in] double_pitch integer �����ĸ�����������
	 * @param[in] freedom_rho double[] rho���ɶ�
	 * @param[in] freedom_rhou double[] rhou���ɶ�
	 * @param[in] freedom_rhov double[] rhov���ɶ�
	 * @param[in] freedom_rhoE double[] rhoE���ɶ�
	 * @param[in] convar_rho_vol double[] ��Ԫ�����غ���rho��ֵ
	 * @param[in] convar_rhou_vol double[] ��Ԫ�����غ���rhou��ֵ
	 * @param[in] convar_rhov_vol double[] ��Ԫ�����غ���rhov��ֵ
	 * @param[in] convar_rhoE_vol double[] ��Ԫ�����غ���rhoE��ֵ
	 * @param[in] convar_rho_edge double[] ��Ԫ�����غ���rho��ֵ
	 * @param[in] convar_rhou_edge double[] ��Ԫ�����غ���rhou��ֵ
	 * @param[in] convar_rhov_edge double[] ��Ԫ�����غ���rhov��ֵ
	 * @param[in] convar_rhoE_edge double[] ��Ԫ�����غ���rhoE��ֵ
	 * @param[in] vol_bf_value double[] �������ڵ�Ԫ���ڵ�ֵ
	 * @param[in] edge_bf_value double[] �������ڵ�Ԫ���ϵ�ֵ
	 */
	__global__ void kernel_calculateConVars(int tnum, int double_pitch,
		const double *freedom_rho,  const double *freedom_rhou,
		const double *freedom_rhov, const double *freedom_rhoE,

		double *convar_rho_vol,   double *convar_rhou_vol,
		double *convar_rhov_vol,  double *convar_rhoE_vol,

		double *convar_rho_edge,  double *convar_rhou_edge, 
		double *convar_rhov_edge, double *convar_rhoE_edge, 

		const double *vol_bf_value, const double *edge_bf_value
		);

	/** 
	 * ����߽�����
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] num integer  �����е�Ԫ����
	 * @param[in] double_pitch integer  ������˫����������
	 * @param[in] rho double Զ���ܶ�
	 * @param[in] rhou double Զ���غ���rhou
	 * @param[in] rhov double Զ���غ���rhov
	 * @param[in] rhoE double Զ���غ���rhoE
	 * @param[in|out] freedom_rho double[] rho���ɶ�
	 * @param[in|out] freedom_rhou double[] rhou���ɶ�
	 * @param[in|out] freedom_rhov double[] rhov���ɶ�
	 * @param[in|out] freedom_rhoE double[] rhoE���ɶ�
	 * @param[in|out] convar_rho_edge double[] �߽����غ���rho��ֵ
	 * @param[in|out] convar_rhou_edge double[] �߽����غ���rhou��ֵ
	 * @param[in|out] convar_rhov_edge double[] �߽����غ���rhov��ֵ
	 * @param[in|out] convar_rhoE_edge double[] �߽����غ���rhoE��ֵ
	 * @param[in] neighbour int[] ��Ԫ�ھ���Ϣ
	 * @param[in] sharedEdge int[] ��Ԫ�������Ϣ
	 * @param[in] triangle_flag int[] ��Ԫ���
	 * @param[in] outer_normal double[] ��Ԫ���ⷨ����
	 */
	__global__ void kernel_boundaryCondition(int tnum, int num, int double_pitch, 
		double rho, double rhou, double rhov, double rhoE,

		double *convar_rho_edge,  double *convar_rhou_edge, 
		double *convar_rhov_edge, double *convar_rhoE_edge,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE,

		const int *neighbour, const int *sharedEdge, 

		const int *triangle_flag, const double *outer_normal
		);

	/** 
	 * ����������غ�в�
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] double_pitch integer  ������˫���ȸ�����������
	 * @param[in] gamma double ������ȱ�
	 * @param[in] convar_rho_vol double[] ��Ԫ���غ���rho��ֵ
	 * @param[in] convar_rhou_vol double[] ��Ԫ���غ���rhou��ֵ
	 * @param[in] convar_rhov_vol double[] ��Ԫ���غ���rhov��ֵ
	 * @param[in] convar_rhoE_vol double[] ��Ԫ���غ���rhoE��ֵ
	 * @param[out] rhs_volume_rho double[] �غ���rho�ڵ�Ԫ����ֲв�
	 * @param[out] rhs_volume_rhou double[] �غ���rhou�ڵ�Ԫ����ֲв�
	 * @param[out] rhs_volume_rhov double[] �غ���rhov�ڵ�Ԫ����ֲв�
	 * @param[out] rhs_volume_rhoE double[] �غ���rhoE�ڵ�Ԫ����ֲв�
	 * @param[in] vol_gauss_weight double[] ����ָ�˹Ȩ��
	 * @param[in] vol_bf_value double[] �����������ڵ�Ԫ���ϵ�ֵ
	 */
	__global__ void kernel_calculateVolumeRHS(int tnum, int double_pitch, double gamma,

		const double *convar_rho_vol,  const double *convar_rhou_vol, 
		const double *convar_rhov_vol, const double *convar_rhoE_vol,

		double *rhs_volume_rho,  double *rhs_volume_rhou, 
		double *rhs_volume_rhov, double *rhs_volume_rhoE,

		const double *vol_gauss_weight, const double *vol_bdf_value
		);

	/** 
	 * ����LFͨ��ϵ��
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] int_pitch integer  ����������������
	 * @param[in] double_pitch integer  ������˫����������
	 * @param[in] outer_normal_vector double [] ��Ԫ�ߵ��ⷨ����
	 * @param[in] neighbour int[] ��Ԫ�������Ϣ
	 * @param[in] freedom_rho double[] �غ���rho���ɶ�
	 * @param[in] freedom_rhou double[] �غ���rhou���ɶ�
	 * @param[in] freedom_rhov double[] �غ���rhov���ɶ�
	 * @param[in] freedom_rhoE double[] �غ���rhoE���ɶ�
	 * @param[out] lfflux_coeff double[] LFͨ��ϵ��
	 */
	__global__ void kernel_calculateLFCoeff(
		int tnum, int int_pitch, int double_pitch, double gamma, 

		const double *outer_normal_vector,  const int *neighbour, 

		const double* freedom_rho, const double *freedom_rhou, 
		const double *freedom_rhov, const double *freedom_rhoE, 

		double *lfflux_coeff
		);

	/** 
	 * �������FG��ֵ
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] num integer  ��������
	 * @param[in] double_pitch integer  ������������
	 * @param[in] gamma double ������ȱ�
	 * @param[in] convar_rho_edge double[] ��Ԫ�����غ���rho��ֵ
	 * @param[in] convar_rhou_edge double[] ��Ԫ�����غ���rhou��ֵ
	 * @param[in] convar_rhov_edge double[] ��Ԫ�����غ���rhov��ֵ
	 * @param[in] convar_rhoE_edge double[] ��Ԫ�����غ���rhoE��ֵ
	 * @param[out] fedge_rho double[] ��Ԫ����fͨ��
	 * @param[out] fedge_rhou double[] ��Ԫ����fͨ��
	 * @param[out] fedge_rhov double[] ��Ԫ����fͨ��
	 * @param[out] fedge_rhoE double[] ��Ԫ����fͨ��
	 * @param[out] gedge_rho double[] ��Ԫ����fͨ��
	 * @param[out] gedge_rhou double[] ��Ԫ����fͨ��
	 * @param[out] gedge_rhov double[] ��Ԫ����fͨ��
	 * @param[out] gedge_rhoE double[] ��Ԫ����fͨ��
	 */
	__global__ void kernel_calculateEdgeFG(
		 int tnum, int num, int double_pitch, double gamma, 

		 const double *convar_rho_edge, const double *convar_rhou_edge,
		 const double *convar_rhov_edge, const double *convar_rhoE_edge,

		 double *fedge_rho, double *fedge_rhou,
		 double *fedge_rhov, double *fedge_rhoE,

		 double *gedge_rho, double *gedge_rhou,
		 double *gedge_rhov, double *gedge_rhoE
		);


	/** 
	 * �������ͨ����ֵ
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] int_pitch integer  ���������������
	 * @param[in] double_pitch integer �����˫����������
	 * @param[in] neighbour int[] ��Ԫ�ھ���Ϣ
	 * @param[in] sharedEdge int[] ��Ԫ�������Ϣ
	 * @param[in] convar_rho_edge double[] ��Ԫ�����غ���rho��ֵ
	 * @param[in] convar_rhou_edge double[] ��Ԫ�����غ���rhou��ֵ
	 * @param[in] convar_rhov_edge double[] ��Ԫ�����غ���rhov��ֵ
	 * @param[in] convar_rhoE_edge double[] ��Ԫ�����غ���rhoE��ֵ
	 * @param[in] fedge_rho double[] ��Ԫ����fͨ��
	 * @param[in] fedge_rhou double[] ��Ԫ����fͨ��
	 * @param[in] fedge_rhov double[] ��Ԫ����fͨ��
	 * @param[in] fedge_rhoE double[] ��Ԫ����fͨ��
	 * @param[in] gedge_rho double[] ��Ԫ����fͨ��
	 * @param[in] gedge_rhou double[] ��Ԫ����fͨ��
	 * @param[in] gedge_rhov double[] ��Ԫ����fͨ��
	 * @param[in] gedge_rhoE double[] ��Ԫ����fͨ��
	 * @param[in] outer_normal_vector double[] ��Ԫ���ⷨ����
	 * @param[in] lfflux_coeff double[] ��Ԫ���ϵ�LFͨ��ϵ��
	 * @param[out] lfflux_rho double[] ��Ԫ���ϸ����㴦��ͨ��
	 * @param[out] lfflux_rhou double[] ��Ԫ���ϸ����㴦��ͨ��
	 * @param[out] lfflux_rhov double[] ��Ԫ���ϸ����㴦��ͨ��
	 * @param[out] lfflux_rhoE double[] ��Ԫ���ϸ����㴦��ͨ��
	 */
/*	__global__ void kernel_calculateFlux(
		int tnum, int int_pitch, int double_pitch, 

		const int *neighbour, const int *sharedEdge,

		const double *convar_rho_edge, const double *convar_rhou_edge,
		const double *convar_rhov_edge, const double *convar_rhoE_edge,
		
		const double *fedge_rho, const double *fedge_rhou,
		const double *fedge_rhov, const double *fedge_RhoE,
		
		const double *gedge_rho, const double *gedge_rhou,
		const double *gedge_rhov, const double *gedge_rhoE,

		const double *outer_normal_vector, const double *lfflux_coeff,

		double *lfflux_rho, double *lfflux_rhou,
		double *lfflux_rhov, double *lfflux_rhoE
		);
*/
	__global__ void kernel_calculateFlux(
		int tnum, int int_pitch, int double_pitch, 

		const int *neighbour, const int *sharedEdge,

		const double *convar_rho_edge, const double *convar_rhou_edge,
//		const double *convar_rhov_edge, const double *convar_rhoE_edge,
		
		const double *fedge_rho, const double *fedge_rhou,
//		const double *fedge_rhov, const double *fedge_RhoE,
		
		const double *gedge_rho, const double *gedge_rhou,
//		const double *gedge_rhov, const double *gedge_rhoE,

		const double *outer_normal_vector, const double *lfflux_coeff,

		double *lfflux_rho, double *lfflux_rhou
//		double *lfflux_rhov, double *lfflux_rhoE
		);
	
	/** 
	 * ���㵥Ԫ�������ϵĻ��ֲв�
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] double_pitch integer  ������˫����������
	 * @param[in] edge_gauss_weight double[] ��Ԫ�ߵĻ���Ȩ��
	 * @param[in] edge_bf_value double[] �������ڱ��ϵ�ֵ
	 * @param[out] rhs_edge_rho double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[out] rhs_edge_rhou double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[out] rhs_edge_rhov double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[out] rhs_edge_rhoE double[] ��Ԫ���ϵĻ��ֲв�
	 */
	__global__ void kernel_calculateEdgeRHS(
		int tnum, int double_pitch, 

		const double *edge_gauss_weight, const double *edge_bf_value,

		const double *lfflux_rho,  const double *lfflux_rhou, 
		const double *lfflux_rhov, const double *lfflux_rhoE, 

		double *rhs_edge_rho,  double *rhs_edge_rhou, 
		double *rhs_edge_rhov, double *rhs_edge_rhoE,

		const double *rhs_volume_rho,  const double *rhs_volume_rhou, 
		const double *rhs_volume_rhov, const double *rhs_volume_rhoE
		);

	/** 
	 * rkdgʱ���ƽ���һ��
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] double_pitch integer  ������˫����������
	 * @param[in] dt  double ʱ�䲽��
	 * @param[in] mass_coeff double[] ����������ϵ��
	 * @param[in|out] freedom_rho double[] �غ���rho���ɶ�
	 * @param[in|out] freedom_rhou double[] �غ���rhou���ɶ�
	 * @param[in|out] freedom_rhov double[] �غ���rhov���ɶ�
	 * @param[in|out] freedom_rhoE double[] �غ���rhoE���ɶ�
	 * @param[in] rhs_edge_rho double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhou double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhov double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhoE double[] ��Ԫ���ϵĻ��ֲв�
	 */
	__global__ void kernel_rkdgStepOne(
		int tnum, int double_pitch, double dt, const double *mass_coeff,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE, 

		const double *rhs_edge_rho,  const double *rhs_edge_rhou,
		const double *rhs_edge_rhov, const double *rhs_edge_rhoE
		);

	/** 
	 * rkdgʱ���ƽ��ڶ���
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] double_pitch integer  ������˫����������
	 * @param[in] dt  double ʱ�䲽��
	 * @param[in] mass_coeff double[] ����������ϵ��
	 * @param[out] freedom_rho double[] �غ���rho���ɶ�
	 * @param[out] freedom_rhou double[] �غ���rhou���ɶ�
	 * @param[out] freedom_rhov double[] �غ���rhov���ɶ�
	 * @param[out] freedom_rhoE double[] �غ���rhoE���ɶ�
	 * @param[in] rhs_edge_rho double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhou double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhov double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhoE double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] freedom_rho_old  double[] rho�����ɶ�
	 * @param[in] freedom_rhou_old double[] rhou�����ɶ�
	 * @param[in] freedom_rhov_old double[] rhov�����ɶ�
	 * @param[in] freedom_rhoE_old double[] rhoE�����ɶ�
	 */
	__global__ void kernel_rkdgStepTwo(
		int tnum, int double_pitch, double dt, const double *mass_coeff,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE, 

		const double *rhs_edge_rho,  const double *rhs_edge_rhou,
		const double *rhs_edge_rhov, const double *rhs_edge_rhoE,

		const double *freedom_rho_old,  const double *freedom_rhou_old, 
		const double *freedom_rhov_old, const double *freedom_rhoE_old
		);

	/** 
	 * rkdgʱ���ƽ�������
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] double_pitch integer  ������˫����������
	 * @param[in] dt  double ʱ�䲽��
	 * @param[in] mass_coeff double[] ����������ϵ��
	 * @param[in|out] freedom_rho double[] �غ���rho���ɶ�
	 * @param[in|out] freedom_rhou double[] �غ���rhou���ɶ�
	 * @param[in|out] freedom_rhov double[] �غ���rhov���ɶ�
	 * @param[in|out] freedom_rhoE double[] �غ���rhoE���ɶ�
	 * @param[in] rhs_edge_rho double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhou double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhov double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] rhs_edge_rhoE double[] ��Ԫ���ϵĻ��ֲв�
	 * @param[in] freedom_rho_old double[] rho�����ɶ�
	 * @param[in] freedom_rhou_old double[] rhou�����ɶ�
	 * @param[in] freedom_rhov_old double[] rhov�����ɶ�
	 * @param[in] freedom_rhoE_old double[] rhoE�����ɶ�
	 */
	__global__ void kernel_rkdgStepThree(
		int tnum, int double_pitch, double dt, const double *mass_coeff,

		double *freedom_rho,  double *freedom_rhou, 
		double *freedom_rhov, double *freedom_rhoE, 

		const double *rhs_edge_rho,    const double *rhs_edge_rhou,
		const double *rhs_edge_rhov,   const double *rhs_edge_rhoE,

		const double *freedom_rho_old,     const double *freedom_rhou_old, 
		const double *freedom_rhov_old,    const double *freedom_rhoE_old
		);

	/** 
	 * ����в�
	 * @param[in] tnum integer ������������ʵ��������Ŀ
	 * @param[in] freedom_rho double[] �غ���rho�����ɶ�
	 * @param[in] freedom_rhoE double[] �غ���rhoE�����ɶ�
	 * @param[in] freedom_rho_old  double[] �����ɶ�
	 * @param[in] freedom_rhoE_old double[] �����ɶ�
	 * @param[out] residual double[] �в�
	 */
	__global__ void kernel_calculateResidual(
		int tnum, 

		const double *freedom_rho,     const double *freedom_rhoE,
		const double *freedom_rho_old, const double *freedom_rhoE_old,

		double *residual
		);


	__global__ void kernel_airfoil(
		double gamma, double alpha, double *cl, double *cd,
		int tnum, int double_pitch, int int_pitch, double q_inf,
		const int *airfoil, const double *outer_normal_vector,
		const int *neighbour, const double *edge_gauss_weight,
		const double *freedom_rho,     const double *freedom_rhou,
		const double *freedom_rhov, const double *freedom_rhoE,
		const double *edge_bf_value
		);

	__device__ double uh(
		int e, int double_pitch, const double* ul, 
		int gauss_index, const double *edge_bf_value);

};
