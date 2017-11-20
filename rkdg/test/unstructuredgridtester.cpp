#include "unstructuredgridtester.h"

CUnstructuredGridTester::CUnstructuredGridTester()
{

}

template <typename T>
void CUnstructuredGridTester::output2DVar( int rnum, int cnum, const T *var, const string &output )
{
	ofstream fout(output.c_str());
	if (!fout)
	{
		cout<<"文件： "<<output<<" 打开失败，文件输出过程将被忽略"<<endl;
		return;
	}

	fout.width(15);
	fout.precision(12);
	for ( int i=0; i<rnum; ++i )
	{
		for ( int j=0; j<cnum; ++j )
		{
			fout<<var[i+j*rnum]<<"   ";
		}
		fout<<endl;
	}

	fout.close();
}

void CUnstructuredGridTester::testNeighbour( const string& output/*="output/neighbour.dat"*/ )
{
	output2DVar(getCellNumber(), TRIANGLE_EDGES, tri_neighbour, output);
}

void CUnstructuredGridTester::testSharedEdge( const string& output/*="output/sharededge.dat"*/ )
{
	output2DVar(getCellNumber(), TRIANGLE_EDGES, tri_sharedEdge, output);
}

void CUnstructuredGridTester::testTriangleFlag( const string& output/*="output/triangleflag.dat"*/ )
{
	output2DVar(getCellNumber(), 1, tri_flag, output);
}

void CUnstructuredGridTester::testArea( const string& output/*="output/area.dat"*/ )
{
	output2DVar(getCellNumber(), 1, triangle_infos.area, output);
}

void CUnstructuredGridTester::testPerimeter( const string& output/*="output/perimeter.dat"*/ )
{
	output2DVar(getCellNumber(), 1, triangle_infos.perimeter, output);
}


void CUnstructuredGridTester::testOuterNormalVector( const string& output/*="output/outernormalvector.dat"*/ )
{
	output2DVar(getCellNumber(), TRIANGLE_EDGES*2, triangle_infos.outer_normal_vector, output);
}

void CUnstructuredGridTester::testVolumeBfValue( void )
{
	
	int num = getCellNumber();
	char buffer[10];

	for ( int i=0; i<BASIS_FUNCTIONS; ++i )
	{
		string filename = "output/volumebfvalue-f";

		_itoa(i, buffer,10);
		filename += buffer;
		filename += ".dat";

		output2DVar(num, VOLUME_GPOINTS, &triangle_infos.vol_bf_value[i*VOLUME_GPOINTS*num], filename);
	}
}

void CUnstructuredGridTester::testMassCoeff( const string& output/*="output/masscoeff.dat"*/ )
{
	output2DVar(getCellNumber(), BASIS_FUNCTIONS, triangle_infos.mass_coeff, output);
}


void CUnstructuredGridTester::testVolumeBdfValue( void )
{

	int num = getCellNumber();
	char buffer[10];

	for ( int i=0; i<2*BASIS_FUNCTIONS; ++i )
	{
		string filename = "output/volumebdfvalue-f";
		_itoa(i, buffer, 10);

		filename += buffer;
		filename += ".dat";

		output2DVar(num, VOLUME_GPOINTS, &triangle_infos.vol_bdf_value[i*VOLUME_GPOINTS*num], filename);

	}
}


void CUnstructuredGridTester::testEdgeBfValue( void )
{
	int num = getCellNumber();
	char buffer[10];

	for ( int i=0; i<BASIS_FUNCTIONS; ++i )
	{
		
		for ( int j=0; j<TRIANGLE_EDGES; ++j )
		{
			string filename = "output/edgebfvalue-f";

			_itoa(i, buffer,10);
			filename += buffer;

			filename += "-e";
			_itoa(j, buffer,10);
			filename += buffer;
			filename += ".dat";

			output2DVar(num, EDGE_GPOINTS, &triangle_infos.edge_bf_value[j*EDGE_GPOINTS*num+i*TRIANGLE_EDGES*EDGE_GPOINTS*num], filename);
		}
	}
}

void CUnstructuredGridTester::testVolumeGaussWeight( const string& output/*="output/volumegaussweight.dat"*/ )
{
	output2DVar(getCellNumber(), VOLUME_GPOINTS, triangle_infos.vol_gauss_weight, output);
}

void CUnstructuredGridTester::testEdgeGaussWeight( const string& output/*="output/edgegaussweight.dat"*/ )
{
	output2DVar(getCellNumber(), TRIANGLE_EDGES*EDGE_GPOINTS, triangle_infos.edge_gauss_weight, output);
}


CUnstructuredGridTester::~CUnstructuredGridTester()
{

}

