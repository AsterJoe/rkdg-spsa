/**
 * ���ļ��ǳ����������ͷ�ļ�
 * ͨ���������ܹ�ʵ�ִ��ļ��ж�ȡ����������
 * @author tlanyan<tag.yuan@gmail.com>
 * @link http://tlanyan.me
 * @copyright Copyright &copy; 2013-2015 tlanyan
 * =============================================
 * @attention ����������Ƶ�ʹ�ñ��ļ�����Դ�Ļ���ҵ��Ŀ�ģ���Ψһ��Ҫ���Ǳ���������Ϣ����Ȩ˵���Լ���ʹ��ע�⡣
 */

#pragma once

#include "cppstdheaders.h"
#include "myexception.h"
using namespace std;


extern "C"{
	/**
	 * �����е�����: �޷�ʶ����У�ע���У����У��������У� ��֧�ֵ�������
	 * ���в�֧�ֵ�������ָ��������������Ǹ��������ڳ����в�ʹ�á�
	 */
	enum cl_type{
		CL_TYPE_UNKNOWN   = 0,
		CL_TYPE_COMMENT   = 1,
		CL_TYPE_EMPTY     = 2,
		CL_TYPE_CONFITEM  = 3,
		CL_TYPE_UNSUPPORT = 4
	};

}

/**
 * �ļ������࣬���ඨ����Щ���ÿ��Դ��ļ��л�ȡ����ֵ
 * ͨ������������ļ����ƽ��������ļ�
 * �����ļ������Ľ��������_message��
 * @example �����÷�����
 * // �ȶ�����Щ����Ҫ���ļ��ж�ȡ
 * string items[] = {"ia", "ib", "ic"}
 * // ʹ��Ĭ�Ϲ��캯������������
 * CConfig conf("path/to/configfile", items, 3) // ����3������3������������ļ��ж�ȡ
 * if ( !conf.parseConfigFile() )
 *	// �Խ������ɹ��Ĵ���
 *  // ������ļ��ж�ȡ��ֵ
 * conf.printConfigItems();
 * // ��ȡ�ļ�������Ϣ
 * conf.getMessage()
 */
class CConfig{
// ����
protected:
	/**
	 * ���ڴ�ӡ����Ϣ
	 */
	string _message;

	/** 
	 * ������飬��Ǵ�����������Ƿ����ļ��г�ʼ����
	 */
	bool *item_flag;

	/** 
	 * ��Ƿ��ţ�ָʾ�Ƿ�����������ɹ���ʼ��
	 */
	bool all_items_parsed;

public:
    /**
     * �����ļ���
     */
    string config_file;

	/**
	 * ֧�ֵ�����������
	 * ��γ�ʼ���������鿴 $Ĭ�Ϲ��캯��$ ��˵��
	 */
	map<string, string> config_items;

	/** 
	 * ������ĸ���
	 */
	int item_num;
	 
// ����
protected:
	 /**
     * ���캯��
	 * �ù��캯�������ã�����ֱ�ӵ��øú���
     */
    CConfig(){};

	 /** 
	  * �����ļ���ȡ֮�����Ƿ������������Ѿ��������ļ��еõ��˳�ʼ��
	  */
	virtual void afterParseConfigFile();

public:
	/**
     * Ĭ�Ϲ��캯��
	 * @param string filename �����ļ���
	 * @param string[] items ��Ҫ���ļ��ж�ȡ��������
	 * @param n ������ĸ���
	 */
	CConfig(const string& filename, string items[], int n);

	/**
	 * ���������ļ�
	 * ����ļ������ڣ����׳��쳣
	 * �����ɹ��Ľ������ͨ������ getMessage() ��������
	 * @return boolean ���������������ȫ���������ļ��г����򷵻�true, ���򷵻�false
	 */
	bool parseConfigFile(void);

	/**
	 * �����ļ�������
	 * @param integer lineno �к�
	 * @param string line ������
	 * @return integer ���������е�����
	 */
	int parseConfigLine(int lineno, const string& line);

	/** 
	 * ĳ���������Ƿ�Ϊ���е�������
	 * @param string item ���������Ƿ���֧�ֵ���������
	 * @return integer ����ҵ�������򷵻����������򷵻�-1 
	 */
	int configItemIndex(const string& item);

	/**
	 * ��ӡ����Ϣ,�밴����Ҫ�Զ���˺���
	 */
	void printMessage(void);

	/**
	 * �����ȡ��ֵ
	 */
	void printConfigItems(void);

	/**
	 * ��ȡ�ļ��������
	 */
	string& getMessage(void){ return _message; };

	/** 
	 * ��ӡ���ô��󣬲���ʾ�����˳�
	 */
	void printQuitMessage(void);

    /**
     * ��������
     */
    ~CConfig();

};
