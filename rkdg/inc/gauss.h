/**
 * ��˹����Ԫ��Ԫ������
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */

#pragma once
#include "cppstdheaders.h"

using namespace std;

extern "C"{

/**
 * ����Ԫ�ص�Gauss��ȥ�����򵥵����Է�����
 *
 * @param[in]  N : ��˾�������/����
 * @param[in]  A : ��˾���(��һά������ʽ�洢)
 * @param[in]  B : �Ҷ�����,��ⷽ������Ϻ��Ž�����
 */
void Gauss(int N, double* A, double* B);

}