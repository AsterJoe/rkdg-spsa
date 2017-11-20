#include "../inc/unstructuredgrid.h"

CUnstructuredGrid::CUnstructuredGrid():
_vertice_num(0),
_edge_num(0),
_triangle_num(0),
_ghost_triangle_num(0),
_cell_num(0),
tri_neighbour(NULL),
tri_sharedEdge(NULL),
tri_flag(NULL)
{}

void CUnstructuredGrid::initializeGrid(void)
{
	// ��Ҫ�������ļ����ȡ��������
	string conf_items[] = {
		"trianglenum",
		"verticenum", 
		"edgenum", 
		"trivertice_filename",
		"trineigh_filename",
		"triedge_filename", 
		"vertice_filename", 
		"edgevertice_filename",
		"output_filename"
	};
	
	CConfig grid_conf(config_file, conf_items, 9);
	if ( !grid_conf.parseConfigFile() )
		grid_conf.printQuitMessage();
	
	parseGhostNum(grid_conf.config_items["trineigh_filename"]);	
	
	_triangle_num = atoi(grid_conf.config_items["trianglenum"].c_str());
	_vertice_num  = atoi(grid_conf.config_items["verticenum"].c_str());
	_edge_num = atoi(grid_conf.config_items["edgenum"].c_str());
	output_filename = grid_conf.config_items["output_filename"];
	_cell_num = _triangle_num + _ghost_triangle_num;

	// ����ռ�
	vertice.resize((_vertice_num+_ghost_triangle_num));
	edge.resize(_edge_num);
	tri_vertice.resize(3*_cell_num);
	tri_edge.resize(3*_cell_num);
	tri_neighbour  = new int[3*_cell_num];
	tri_sharedEdge = new int[3*_cell_num];
	tri_flag       = new int[_cell_num];
	ver_flag       = new int[_vertice_num];

	readTriangleVertice(grid_conf.config_items["trivertice_filename"]);
	readTriangleNeighbour(grid_conf.config_items["trineigh_filename"]);
	readTriangleEdge(grid_conf.config_items["triedge_filename"]);
	readVertice(grid_conf.config_items["vertice_filename"]);
	readEdgeVertice(grid_conf.config_items["edgevertice_filename"]);

	initGhostGrid();
	
	initSharedEdge();

	getUpperNode2D();
}

void CUnstructuredGrid::initializeTriangleInfos(void)
{
	int i, j;
	double x[3], y[3];

	for ( i=0; i<_triangle_num; ++i )
	{
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getY();
		}

		triangle_infos.initInformation(i, x, y);
	}

	for ( ; i<_cell_num; ++i )
	{
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getY();
		}

		triangle_infos.barycenter[i].setVertice((x[0]+x[1]+x[2])/3.0, (y[0]+y[1]+y[2])/3.0);
	}
}

void CUnstructuredGrid::parseGhostNum(const string& trineigh_filename)
{
	ifstream fin(trineigh_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open neighbour file: "+trineigh_filename);
	
	int s(0);
	string tmp_str;
	while ( fin>>tmp_str )
	{
		if ( tmp_str=="-1" )
			++ s;
	}

	_ghost_triangle_num = s;

	fin.close();
}

void CUnstructuredGrid::readTriangleVertice(const string& trivertice_filename)
{
	ifstream fin(trivertice_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open triangle vertice file: "+trivertice_filename);
	
	string tmp_string;
	getline(fin, tmp_string);
	
	// ��ȡ��ɵ�Ԫ�Ķ�����
	for ( int i=0; i<_triangle_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_VERTICES; ++j )
			fin>>tri_vertice[TRIANGLE_VERTICES*i+j];
			
		tri_flag[i] = CELL_INTERIOR;
	}

	fin.close();
}

void CUnstructuredGrid::readTriangleNeighbour(const string& trineigh_filename)
{
	ifstream fin(trineigh_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open triangle neighbour file: "+trineigh_filename);

	string tmp_string;
	getline(fin, tmp_string);

	// ��ȡ��Ԫ�����ڵ�Ԫ���
	for ( int i=0; i<_triangle_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_EDGES; ++j )
			fin>>tri_neighbour[i+j*_cell_num];
	}

	fin.close();
}

void CUnstructuredGrid::readTriangleEdge(const string& triedge_filename)
{
	ifstream fin(triedge_filename.c_str());
	
	if ( !fin )
		throw CMyException("Failed to open triangle edge file: "+triedge_filename);

	string tmp_string;
	getline(fin, tmp_string);

	// ��ȡ����߱��
	for ( int i=0; i<_triangle_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_EDGES; ++j )
			fin>>tri_edge[i*TRIANGLE_EDGES+j];
	}

	fin.close();
}

void CUnstructuredGrid::readVertice(const string& vertice_filename)
{
	ifstream fin(vertice_filename.c_str());	
	if ( !fin )
		throw CMyException("Failed to open vertice file: "+vertice_filename);

	string tmp_string;
	getline(fin, tmp_string);
	
	double x, y;
	// ��ȡ��������
	for ( int i=0; i<_vertice_num; ++i )
	{
		fin>>x>>y;
		vertice[i].setVertice(x,y);
		ver_flag[i] = VERTICE_INTERIOR;
	}

	fin.close();
}

void CUnstructuredGrid::readEdgeVertice(const string& edgevertice_filename)
{
	ifstream fin(edgevertice_filename.c_str());
	if ( !fin )
		throw CMyException("Failed to open edge vertice file: "+edgevertice_filename);

	string tmp_string;
	getline(fin, tmp_string);

	int start, terminal;
	// ��ȡ��������
	for ( int i=0; i<_edge_num; ++i )
	{
		fin>>start>>terminal;
		edge.at(i).setEdge(start,terminal);
	}

	fin.close();
}

void CUnstructuredGrid::initGhostGrid(void)
{
	int i, j, k, e;
	double x1, x2, x3, y1, y2, y3, xt, yt;
	
	int ten = _triangle_num;
	int tvn = _vertice_num;

	set<int> wall_temp;
	set<int> far_temp;
	for ( e=0; e<_triangle_num; ++e )
	{
		if ( tri_neighbour[e]==-1 )
		{
			// ��Ԫ��ź͵�����
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			if (x1 * x1 + y1 * y1 < 4) {
				wall_temp.insert(j);
				wall_temp.insert(k);
				wall_elem.push_back(e);
			} else {
				far_temp.insert(j);
				far_temp.insert(k);
			}
			ver_flag[j] = VERTICE_BOUNDARY;
			ver_flag[k] = VERTICE_BOUNDARY;
			// �������������
			xt = x2 + x3 - x1;
			yt = y2 + y3 - y1;

			vertice[tvn].setVertice(xt,yt);

			// ���µ�Ԫ�б�
			tri_vertice[ten*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten*TRIANGLE_VERTICES+1] = k;
			tri_vertice[ten*TRIANGLE_VERTICES+2] = j;
			
			// ������ⵥԪ���ھӹ�ϵ
			tri_neighbour[ten] = e;
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 0;

			// �������ڵ�Ԫ�ھ��б�
			tri_neighbour[e] = ten;

			++ tvn;
			++ ten;
		}
		else if ( tri_neighbour[e+_cell_num]==-1 )
		{
			// ��Ԫ��ź͵�����
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			if (x1 * x1 + y1 * y1 < 4) {
				wall_temp.insert(i);
				wall_temp.insert(k);
				wall_elem.push_back(e);
			} else {
				far_temp.insert(i);
				far_temp.insert(k);
			}
			ver_flag[i] = VERTICE_BOUNDARY;
			ver_flag[k] = VERTICE_BOUNDARY;

			// ������������
			xt = x1 + x3 - x2;
			yt = y1 + y3 - y2;
			
			vertice[tvn].setVertice(xt,yt);

			// ���µ�Ԫ�б�
			tri_vertice[ten*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten*TRIANGLE_VERTICES+1] = i;
			tri_vertice[ten*TRIANGLE_VERTICES+2] = k;
			
			// ������ⵥԪ���ھӹ�ϵ
			tri_neighbour[ten] = e;
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 1;

			// �������ڵ�Ԫ�ھ��б�
			tri_neighbour[e+_cell_num] = ten;

			++ tvn;
			++ ten;
		}
		else if ( tri_neighbour[e+2*_cell_num]==-1 )
		{
			// ��Ԫ��ź͵�����
			i = tri_vertice[e*TRIANGLE_VERTICES];
			j = tri_vertice[e*TRIANGLE_VERTICES+1];
			k = tri_vertice[e*TRIANGLE_VERTICES+2];

			x1 = vertice[i].getX();
			x2 = vertice[j].getX();
			x3 = vertice[k].getX();
			y1 = vertice[i].getY();
			y2 = vertice[j].getY();
			y3 = vertice[k].getY();

			if (x1 * x1 + y1 * y1 < 4) {
				wall_elem.push_back(e);
				wall_temp.insert(i);
				wall_temp.insert(j);
			} else {
				far_temp.insert(i);
				far_temp.insert(j);
			}
			ver_flag[i] = VERTICE_BOUNDARY;
			ver_flag[j] = VERTICE_BOUNDARY;

			// ������������
			xt = x1 + x2 - x3;
			yt = y1 + y2 - y3;
			vertice[tvn].setVertice(xt,yt);

			// ���µ�Ԫ�б�
			tri_vertice[ten*TRIANGLE_VERTICES] = tvn;
			tri_vertice[ten*TRIANGLE_VERTICES+1] = i;
			tri_vertice[ten*TRIANGLE_VERTICES+2] = j;

			// ������ⵥԪ���ھӹ�ϵ
			tri_neighbour[ten] = e;
			tri_neighbour[ten+_cell_num] = -1;
			tri_neighbour[ten+2*_cell_num] = -1;
			tri_sharedEdge[ten] = 2;
			
			// �������ڵ�Ԫ�ھ��б�
			tri_neighbour[e+2*_cell_num] = ten;

			++ tvn;
			++ ten;
		}
	}

	set<int>::iterator iterCur = wall_temp.begin();
	set<int>::iterator iterEnd = wall_temp.end();
	for (; iterCur != iterEnd; ++iterCur)
	{
		wall.push_back(*iterCur);
	}
	iterCur = far_temp.begin();
	iterEnd = far_temp.end();
	for (; iterCur != iterEnd; ++iterCur)
	{
		farfiled.push_back(*iterCur);
	}
	ofstream a("output/wallelem.dat");
	airfoil = new int[wall_elem.size()];
	for (int i = 0; i < wall_elem.size(); i++) {
		airfoil[i] = wall_elem[i];
	}
}

int CUnstructuredGrid::getGhostTriangleNumber(void) const
{
	return _ghost_triangle_num;
}

int CUnstructuredGrid::getEdgeNumber() const
{
	return _edge_num;
}

int CUnstructuredGrid::getTriangleNumber() const
{
	return _triangle_num;
}

int CUnstructuredGrid::getVerticeNumber() const
{
	return _vertice_num;
}

int CUnstructuredGrid::getCellNumber() const
{
	return _cell_num;
}

void CUnstructuredGrid::markBoundaryTriangles( void )
{
	double bcx, bcy;

	for ( int i=_triangle_num; i<_cell_num; ++i )
	{
		// ��������
		bcx = triangle_infos.getBarycenter(i).getX();
		bcy = triangle_infos.getBarycenter(i).getY();

		//if ( bcx<-15.0 ) // ������
		//{
		//	tri_flag[i] = CELL_FARFIELD;
		//} 
		//else if ( bcx>15.0 ) // �ҷ�ȥ��
		//{
		////	_triangles[i].setFlag(bnd_outflow);
		//	tri_flag[i] = CELL_FARFIELD;
		//}
		//else if ( bcy<-15.0 ) // �·�����
		//{
		//	tri_flag[i] = CELL_FARFIELD;
		//}
		//else if ( bcy>15.0 ) // �Ϸ�����
		//{
		//	tri_flag[i] = CELL_FARFIELD;
		//}
		if ((bcx * bcx + bcy * bcy) > 4) {
			tri_flag[i] = CELL_FARFIELD;
		}
		else
		{
			tri_flag[i] = CELL_REFLECTION; // ����߽�
		}
	}
}

void CUnstructuredGrid::testTrianglesAntiwise() const
{
	int i, j;
	double x[3], y[3];
	double delt;

	for ( i=0; i<_triangle_num; ++i )
	{
		// ����
		for ( j=0; j<TRIANGLE_VERTICES; ++j )
		{
			x[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getX();
			y[j] = vertice[tri_vertice[i*TRIANGLE_VERTICES+j]].getY();
		}
		// ��������ʽ�ж��Ƿ�����
		delt = (x[1]-x[0])*(y[2]-y[1]) - (y[1]-y[0])*(x[2]-x[1]);

		if ( delt<0 )
		{
			cout<<"The order of vertices in  "<<i+1<<"th triangle is not antiwise."<<endl;
		}
	}
}

void CUnstructuredGrid::outputGrid(void) const
{
	// ���������Ԫ��Ϣ
	ofstream fout(output_filename.c_str());
	if ( !fout )
	{
		cout<<"\nFailed to open grid output file: "<<output_filename<<endl;
		cout<<"Grid output will be omitted."<<endl;
		return;
	}
	
	cout<<"\n========================================\nGrid information:"<<endl;
	cout<<"triangle number: "<<_triangle_num<<endl;
	cout<<"vertice number:  "<<_vertice_num<<endl;
	cout<<"edge number:     "<<_edge_num<<endl;
	cout<<"========================================"<<endl;
		
	fout<<"TITLE=RKDG"<<endl;
	fout<<"VARIABLES="<<"X"<<" , "<<"Y"<< endl;
	fout<<"ZONE N= "<<_vertice_num<<" , "<<"E= "<<_triangle_num<<" , "<<"F=FEPOINT"<<" , "<<"ET=TRIANGLE"<<endl;

	fout.precision(10);
	int i;

	// ���������
	for ( i=0; i<_vertice_num; ++i )
	{
		fout<<vertice.at(i).getX()<<"    "<<vertice.at(i).getY()<<endl;
	}

	// �����Ԫ����
	for ( i=0; i<_triangle_num; ++i )
	{
		// �����б���Ǵ�1��ʼ���ʶ���Ҫ+1
		for ( int j=0; j<TRIANGLE_VERTICES; ++j )
			fout<<tri_vertice[TRIANGLE_VERTICES*i+j]+1<<"    ";
			
		fout<<endl;
	}

	fout.close();
}

void CUnstructuredGrid::outputGridWithGhostCells(const string& filename) const
{
	ofstream fout(filename.c_str(), ios::out);
	if ( !fout )
	{
		cout<<"Failed to open ghost mesh file: "<<filename<<" and output will be omitted."<<endl<<endl;
		return;
	}
	int vnum = _vertice_num + _ghost_triangle_num;

	fout<<"TITLE=mesh-with-ghost-cells"<<endl;
	fout<<"VARIABLES="<<"x"<<" , Y"<<endl;
	fout<<"ZONE N="<<vnum<<" , E="<<_cell_num<<" , F=FEPOINT , ET=TRIANGLE"<<endl;

	int i;
	for ( i=0; i<vnum; ++i )
	{
		fout<<vertice[i].getX()<<"    "<<vertice[i].getY()<<endl;
	}

	for ( i=0; i<_cell_num; ++i )
	{
		for ( int j=0; j<TRIANGLE_VERTICES; ++j )
			fout<<tri_vertice[i*TRIANGLE_VERTICES+j]+1<<"    ";
		
		fout<<endl;
	}

	fout.close();
}

void CUnstructuredGrid::initSharedEdge()
{
	int neigh_index, edge_index;

	for ( int e=0; e<_triangle_num; ++e )
	{
		for ( int k=0; k<TRIANGLE_EDGES; ++k )
		{
			neigh_index = tri_neighbour[e+k*_cell_num];

			// ����õ�ԪΪ���ⵥԪ���������ⵥԪ���ɹ��򣬱ض���0�ű�
			if ( neigh_index>=_triangle_num )
			{
				tri_sharedEdge[e+k*_cell_num] = 0;
			} 
			else
			{
				edge_index = tri_edge[e*TRIANGLE_EDGES+k];
				if ( edge_index==tri_edge[neigh_index*TRIANGLE_EDGES] )
				{
					tri_sharedEdge[e+k*_cell_num] = 0;
				} 
				else if ( edge_index==tri_edge[neigh_index*TRIANGLE_EDGES+1] )
				{
					tri_sharedEdge[e+k*_cell_num] = 1;
				}
				else if ( edge_index==tri_edge[neigh_index*TRIANGLE_EDGES+2] )
				{
					tri_sharedEdge[e+k*_cell_num] = 2;
				}
				else
				{
					cout<<"\nedge relation error:"<<endl;
					cout<<"detail: triangle "<<e+1<<"th's  "<<k+1<<"th edge is the share edge of "<<neigh_index+1<<"th triangle, but indice of "<<neigh_index+1<<"th edge is not correct."<<endl;
					cout<<"press ENTER to exit."<<endl;
					getchar();
					exit(-1);
				}

			}
		}
	}
}

void CUnstructuredGrid::getUpperNode2D() {
	double y_loc;
	for (int i = 0; i < wall.size(); i++)
	{
		y_loc = vertice[wall[i]].getY();
		if (y_loc >= 0) {
			upper_node.push_back(wall[i]);
		}
		if (y_loc <= 0) {
			lower_node.push_back(wall[i]);
		}
	}
	/*for (int i = 0; i < wall.size(); i++) {
		y_loc = vertice[wall.at(i)].getY();
		if (y_loc >= 0) {
			upper_node.push_back(wall.at(i));
		}
		if (y_loc <= 0) {
			lower_node.push_back(wall.at(i));
		}
		}*/
	cout<<"print over"<<endl;

	sortUpperNode2D(&upper_node);
	sortUpperNode2D(&lower_node);
}

void CUnstructuredGrid::sortUpperNode2D(vector<int>* vec) {
	int min, temp;
	for (int i = 0; i < vec->size(); i++) {
		min = i;
		for (int j = i + 1; j < vec->size(); j++) {
			if (vertice[vec->at(min)].getX() > vertice[vec->at(j)].getX()) {
				min = j;
			}
		}
		temp = vec->at(min);
		vec->at(min) = vec->at(i);
		vec->at(i) = temp;
	}
}

double CUnstructuredGrid::calArea() {
	double area(0.0);
	int upperNum, lowerNum;
	int node_a, node_b;

	upperNum = upper_node.size();
	lowerNum = lower_node.size();
	for(int i = 0; i < upperNum-1; ++i){
		node_a = upper_node[i];
		node_b = upper_node[i+1];
		area = area + 0.5 * (vertice[node_a].getY() + vertice[node_b].getY()) * (vertice[node_b].getX() - vertice[node_a].getX());
	}
	//Step 2: Lower side (y <= 0 or maybe y > 0 near trailing edge)
	for(int i = 0; i < lowerNum-1; ++i){
		node_a = lower_node[i];
		node_b = lower_node[i+1];
		area = area + 0.5 * (vertice[node_a].getY() + vertice[node_b].getY()) * (vertice[node_a].getX() - vertice[node_b].getX());
	}
	return area;
}

CUnstructuredGrid::~CUnstructuredGrid()
{
	// �ͷ���Դ
	delete []tri_neighbour;
	delete []tri_sharedEdge;
	delete []tri_flag;
}