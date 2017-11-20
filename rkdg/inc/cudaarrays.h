/**
 * ������Ҫ��GPU��ʹ�õ�����
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://www.tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 
#pragma once

#include "cppstdheaders.h"
#include "cuda_runtime.h"
#include "defines.h"
#include "myexception.h"
using namespace std;

/**
 * ������Ҫ��GPU��ʹ�õ����飬��Щ���齫��GPU�����ʱ��ʹ��
 */
class CCUDAArrays {
	// ����
	public:
		// ��Ԫ��Ϣ
		int *neighbour;							/**< ��Ԫ�ھ���Ϣ */
		int *sharedEdge;						/**< ��Ԫ�������Ϣ */
		int *triangle_flag;						/**< ��Ԫ��� */

		double *area;							/**< ����Ԫ��� */
		double *perimeter;						/**< �ܳ� */
		double *outer_normal_vector;			/**< ��Ԫ�ߵ��ⷨ���� */
		double *mass_coeff;						/**< ������������ϵ�� */
		double *vol_bf_value;					/**< �����������˹���ֽڵ��ϵ�ֵ */
		double *vol_bdf_value;					/**< �������ĵ��������˹���ֽڵ��ֵ */
		double *edge_bf_value; 					/**< �������ڱ��ϸ�˹���ֽڵ��ֵ */
		double *vol_gauss_weight;				/**< ��˹����ֵ�Ȩ�� */
		double *edge_gauss_weight;				/**< ��˹�߻��ֵ�Ȩ�� */

		// ��������
		double *freedom_rho;			/**< rho���ɶ� */
		double *freedom_rhou;			/**< rhou���ɶ� */
		double *freedom_rhov;			/**< rhov���ɶ� */
		double *freedom_rhoE;			/**< rhoE���ɶ� */

		double *convar_rho_vol;			/**< �غ���rho�����˹���ֽڵ��ֵ */
		double *convar_rhou_vol;		/**< �غ���rhou�����˹���ֽڵ��ֵ */
		double *convar_rhov_vol;		/**< �غ���rhov�����˹���ֽڵ��ֵ */
		double *convar_rhoE_vol;		/**< �غ���rhoE�����˹���ֽڵ��ֵ */
		double *convar_rho_edge;		/**< �غ���rho�ڱ߸�˹���ֽڵ��ֵ */
		double *convar_rhou_edge;		/**< �غ���rhou�ڱ߸�˹���ֽڵ��ֵ */
		double *convar_rhov_edge;		/**< �غ���rhov�ڱ߸�˹���ֽڵ��ֵ */
		double *convar_rhoE_edge;		/**< �غ���rhoE�ڱ߸�˹���ֽڵ��ֵ */

		double *rhs_volume_rho;			/**< rho������ֲв� */
		double *rhs_volume_rhou;		/**< rhou������ֲв� */
		double *rhs_volume_rhov;		/**< rhov������ֲв� */
		double *rhs_volume_rhoE;		/**< rhoE������ֲв� */

		double *rhs_edge_rho;			/**< rho���߻��ֲв� */
		double *rhs_edge_rhou;			/**< rhou���߻��ֲв� */
		double *rhs_edge_rhov;			/**< rhov���߻��ֲв� */
		double *rhs_edge_rhoE;			/**< rhoe���߻��ֲв� */

		double *lfflux_coeff;			/**< LFͨ��ϵ�� */

		double *fedge_rho;				/**< f�ڱ��ϵ�ֵrho */
		double *fedge_rhou;				/**< f�ڱ��ϵ�ֵrhou */
		double *fedge_rhov;				/**< f�ڱ��ϵ�ֵrhov */
		double *fedge_rhoE;				/**< f�ڱ��ϵ�ֵrhoE */
               
		double *gedge_rho;				/**< g�ڱ��ϵ�ֵrho */
		double *gedge_rhou;				/**< g�ڱ��ϵ�ֵrhou */
		double *gedge_rhov;				/**< g�ڱ��ϵ�ֵrhov */
		double *gedge_rhoE;				/**< g�ڱ��ϵ�ֵrhoE */

		double *lfflux_rho;				/**< lfͨ��rho */
		double *lfflux_rhou;			/**< lfͨ��rhou */
		double *lfflux_rhov;			/**< lfͨ��rhov */
		double *lfflux_rhoE;			/**< lfͨ��rhoE */

		// �����ɶ�
		double *freedom_rho_old;			/**< �����ɶ� */
		double *freedom_rhou_old;			/**< �����ɶ� */
		double *freedom_rhov_old;			/**< �����ɶ� */
		double *freedom_rhoE_old;			/**< �����ɶ� */

		double *ddt;					/**< GPU�ϴ��ʱ�䲽���ı��� */
		double *residual;				/**< �в� */

		int *airfoil;
		double *cl;
		double *cd;

	protected:
		size_t _int_pitch;				/**< ������������ */
		size_t _double_pitch;			/**< �������������� */

	// ����
	public:
		/** ���캯�� */
		CCUDAArrays();
		
		/**
		 * ��GPU�Ϸ����ڴ�
		 * @param[in] num integer ����Ԫ����
		 * @param[in] wall_num integer ������浥Ԫ����
		 */
		void allocateMemory(int num, int wall_num);

		/**
		 * ��GPU�ϵ������ʼ��
		 * @param[in] num integer ����Ԫ����
		 */
		//void memsetArrays(int num);

		/** @return �����������䳤�� */
		size_t getIntPitch(void) const;

		/** @return ˫�����������䳤�� */
		size_t getDoublePitch(void) const;
		
		/** �������� */
		~CCUDAArrays();
};