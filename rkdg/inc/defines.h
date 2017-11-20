/**
 * ��������еĳ������ߺ궨��
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */
 
#pragma once 

 /** 
 * ���������εĳ���
 */
#ifndef TRIANGLE_CONSTANTS
	#define TRIANGLE_CONSTANTS
	#define TRIANGLE_VERTICES 3  // �����ζ������
	#define TRIANGLE_EDGES 3     // �����α���
	#define BASIS_FUNCTIONS 6    // ����������
	#define VOLUME_GPOINTS  7    // ����ָ�˹���ֽڵ����
	#define EDGE_GPOINTS  3		 // �߻��ָ�˹���ֽڵ����
	#define BASIS_FUNCTION_COEFFS 6  // ������ϵ������
#endif

// �����غ�������
#ifndef CONSERVATIVE_VARS
#define CONSERVATIVE_VARS 4
#endif

// ���������������
#ifndef RUNGE_KUTTA_STEPS
#define RUNGE_KUTTA_STEPS 3
#endif

// ���û��泣��
#ifndef RESIDUAL_VARS
	#define RESIDUAL_VARS 2  // ����ı�������
	#define BUFFER_PORTIONS 2 // �������м�����
#endif

// �ļ������߻���������󳤶�
#ifndef MAX_LENGTH
	#define MAX_LENGTH 1024 
#endif

/** 
 * ���嵥Ԫ����
 */
enum cell_tags{
	CELL_INTERIOR = 0,  // ��Ԫ���ڲ���Ԫ
	CELL_FARFIELD  = 1,  // Զ���߽�����
	CELL_OUTFLOW  = 2,  // ����
	CELL_REFLECTION = 3  // ����߽�
};

enum ver_tags{
	VERTICE_INTERIOR = 0,
	VERTICE_BOUNDARY = 1
};