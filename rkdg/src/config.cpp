#include "../inc/config.h"

// �������Ĭ�Ϲ��캯��
CConfig::CConfig(const string& filename, string items[], int n):
_message(""),
config_file(filename),
item_flag(NULL),
item_num(n),
all_items_parsed(true)
{
	item_flag = new bool[n];
	
	for ( int i=0; i<n; ++i )
	{
		item_flag[i] = false;
		config_items[items[i]] = "";
	}
}

void CConfig::printMessage(void)
{
	cout<<_message<<endl;
}

bool CConfig::parseConfigFile(void)
{

	ifstream fin(config_file.c_str());
	if (!fin)
		throw CMyException("\nFailed to open configure file "+config_file);

	int line_number = 0;
	string tmp_line;

	_message += "\n======================================\n";
	_message += "    Parse file   "+config_file;
	_message += " \n======================================\n\n";

	while(getline(fin, tmp_line))
	{
		++ line_number;
		parseConfigLine(line_number,tmp_line);
	}

	afterParseConfigFile();

	// �ж��Ƿ��������������ļ��ж�ȡ���
	_message += "\n======================================\n";
	_message += "   Parse file finished!\n";
	_message += "======================================\n";

	fin.close();

	return all_items_parsed;
}

int CConfig::parseConfigLine(int lineno, const string& line)
{
	char buffer[256];

	// ȥ�����׿ո�
	string::size_type first_index = line.find_first_not_of(" \t");

	// ������к�ע����
	if (first_index==string::npos)
		return CL_TYPE_EMPTY;
	if (line.at(first_index)=='#')
		return CL_TYPE_COMMENT;

	// �ҳ�����������
	string separator = ",; #\t";
	string::size_type mid_index = line.find_first_of(" =\t", first_index);
	// ���и���������������ȴû�и�������ֵ
	if (mid_index==string::npos)
		return CL_TYPE_UNKNOWN;
	string item(line, first_index, mid_index-first_index);
	int index = configItemIndex(item);
	if ( -1==index )
	{
		sprintf(buffer, "note: unsurpported item found in line %d : %s\n", lineno, item.c_str());
		_message += string(buffer);
		return CL_TYPE_UNSUPPORT;
	}

	// ȥ��������ֵ��ǰ���ո���ߵ���
	string::size_type another_index = line.find_first_not_of(" =\t", mid_index);
	// ���и���������������ȴû�и�������ֵ
	if (another_index==string::npos)
		return CL_TYPE_UNKNOWN;


	// �ҵ�����ֵ����ֹλ��
	string::size_type final_index = line.find_first_of(separator, another_index);
	string value;
	if ( string::npos==final_index )
	{
		value = line.substr(another_index);
	} 
	else if ( final_index==another_index )
	{
		// û�и�������ֵ
		return CL_TYPE_UNKNOWN;
	}
	else
	{
		value = line.substr(another_index,final_index-another_index);
	}

	config_items[item] = value;
	item_flag[index] = true;
	return CL_TYPE_CONFITEM;
}


int CConfig::configItemIndex(const string& item)
{
	map<string,string>::iterator iter;
	int i;
	for ( i=0,iter=config_items.begin(); iter!=config_items.end(); ++iter,++i )
	{
		if ( iter->first==item )
			return i;
	}

	return -1;
}

void CConfig::printConfigItems(void)
{
	map<string,string>::iterator iter;
	cout<<"\n======================================\n";
	cout<<"    configures in file "+config_file+":\n";
	cout<<"======================================\n";
	int i=0;
	for ( iter=config_items.begin(); iter!=config_items.end(); ++iter, ++i )
	{
		if ( item_flag[i] )
			cout<<iter->first<<" : "<<iter->second<<endl;
	}
	cout<<"\n======================================\n";
}


void CConfig::afterParseConfigFile(void)
{
	map<string,string>::iterator iter;
	int i;
	string tmp;
//	bool unfinished = false;
	tmp += "\n===============  WARNING  ==============\n\n";
	for ( i=0, iter=config_items.begin(); iter!=config_items.end(); ++iter, ++i )
	{
		if ( !item_flag[i] )
		{
//			unfinished = true;
			tmp += "configure item: "+iter->first+" is not be initialized.\n";
			all_items_parsed = false;
		}
	}

	if ( !all_items_parsed )
		_message += tmp;
}

void CConfig::printQuitMessage(void)
{
	printMessage();
	cout<<"Some configure items are not provided initial value."<<endl;
	cout<<"press ENTER to exit."<<endl;
	getchar();
	exit(-1);
}

CConfig::~CConfig()
{
	if ( item_flag!=NULL )
	{
		delete[] item_flag;
	}
}