#include "/home/emir/Рабочий стол/VS/libFolder/my_math.h"

using namespace std;

int main()
{
	
	matrix m;
	std::ifstream IF("data.txt");
    	if(!IF) return -1;
    	IF>>m;
    	m = m.transpose();
    	m.matr.pop_back();
        m = m.transpose();

	//matrix m = {{1,2,3},{2,5,6},{1,8,9},{1,1,1}};
	//matrix m = {{1,6,3},{2,5,6},{12,8,9},{1,19,1}};
	PCA M(m);
	matrix T, P, E;
	try
	{
		std::cout<< "The given matrix:\n"<<M.mat<<std::endl;
		std::cout<<"Centering:\n"<<M.center()<<std::endl;
		
		std::tuple<matrix, matrix, matrix> NIPALS_RESULT = M.NIPALS(12);
		T = std::get<0>(NIPALS_RESULT);
		P = std::get<1>(NIPALS_RESULT);
		E = std::get<2>(NIPALS_RESULT);
		
		std::cout<<"X:\n"<<M.scaling()<<std::endl;
		std::cout<<"T:\n"<<T<<std::endl;
		std::cout<<"P:\n"<<P<<std::endl;
		std::cout<<"E:\n"<<E<<std::endl;
		std::cout<< "X = T * P' + E:\n" << T * (P.transpose()) + E <<std::endl;
		
		std::cout<<"Deviation:\n"<<M.deviation(4);
		std::cout<<"Leverage:\n"<<M.leverage(4);
		
		std::pair<double, double> p = M.dispersion(4);
		std::cout<<"\nERV = " << p.first<<"\tTRV = "<<p.second<<std::endl;	
	}
	catch (exception const& ex) {
		cerr << "Exception: " << ex.what() <<endl;
		return -1;
	}
	catch(...){
		cout<<"Unexpected error."<<endl;
	}
	return 0;
}

/*
vector<vector<double>> v = { {1, 2, 3, 4},{3, 4, 5, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}};
	vector<vector<double>> v1 = { {1, 1, 1, 1},{1, 1, 1, 1}, {2, 2, 2, 2}, {3, 3, 3, 3}, {4, 4, 4, 4}};
	vector<vector<double>> v2 = { {0, 0, 0, 0},{1, 2, 3, 4},{3, 4, 5, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}};
	vector<vector<double>> v3 = {{1,2},{3,4},{5,6}};
	vector<vector<double>> v4 = {{1,3,5},{2,4,6}};
	vector<vector<double>> v5 = {{1, 1,1}, {3,1,5}, {1,4,1}};
	vector<vector<double>> vec = {{1,0,0}};
	vector<vector<double>> vec1 = {{0},{0},{1}};
	vector<vector<double>> vec2 = {{0},{0},{0}};
	vector<vector<double>> vec3 = {{1,2,3}};
	vector<vector<double>> v6 = {{1,2,3},{3,4,9},{5,6,15}};
	vector<vector<double>> v7 = {{1,2},{3,4},{5,6}, {5, 10}};
	
	double mass[5] = {1,2,3,4,5};
	
matrix m;
	matrix m1(v1), m2(v2), m3(v), m4(v3), m5(v4), m6(vec), m7(vec1), m8(v5), m9(vec2), m10(v6), m11(v7);
	identity iden(3);
	diagonal diag(5, mass);
	upper up(v);
	lower down(v);
	symmetric sym(v);
	symmetric symm({ {1, 2, 3, 4},{3, 4, 5, 4}, {1, 2, 3, 4}, {1, 2, 3, 4}});
	matrix MMM = {{1,1,1},{2,2,2},{1.5, 3.3, -90}};
	std::cout<<symm<<std::endl<<sym;
	//std::cout<<m8.inverse() * m8;
	/*std::cout<<m8.det()<<std::endl;
	std::cout<<diag.det()<<std::endl;
	std::cout<<iden.det()<<std::endl;
	std::cout<<m6.scalar(m7)<<std::endl;
	std::cout<<m4.Fronebius()<<std::endl;
	std::cout<<diag.Fronebius()<<std::endl;
	std::cout<<m6.norm()<<' '<<m7.norm()<<std::endl;
	std::cout<<m6.maxnorm()<<' '<<m7.maxnorm()<<std::endl;
	*/
	/*
	std::cout<<m6.transpose()<<std::endl;
	std::cout<<m6.angle(m7)<<std::endl;
	std::cout<<m7.angle(m6)<<std::endl;
	std::cout<<m7.angle(m7)<<std::endl;
	std::cout<<m6.angle(m6)<<std::endl;
	*/
	//std::cout<<sym<<endl;
    	/*
    	std::cout<<m<<endl;
    	std::cout<<iden<<endl;
    	std::cout<<diag<<endl;
    	std::cout<<up<<endl;
    	std::cout<<down<<endl;
    	std::cout<<sym<<endl;
    	std::cout<<"sym + sym\t"<<sym + sym<<endl;
    	std::cout<<"matrix + matrix\t"<<m1 + m2<<endl;
    	std::cout<<"matrix - matrix\t"<<m1 - m2<<endl;
    	std::cout<<"matrix * matrix\t"<<m3 * m3<<endl;
    	std::cout<<"matrix * matrix\t"<<m4*m5<<endl;
    	std::cout<<"matrix * iden\t"<<m3*iden<<endl;
    	std::cout<< "matrix.Adam(matrix)\t"<<m4.Adam(m4)<<endl;
    	std::cout<< "symmetric.Adam(symmetric)\t" <<sym.Adam(sym)<<endl;
    	std::cout<< "symmetric.Adam(up)\t" <<sym.Adam(up)<<endl;
    	std::cout<< "identity.Adam(up)\t" <<iden.Adam(up)<<endl;
    	std::cout<< "down.Adam(identity)\t" <<down.Adam(iden)<<endl;
    	std::cout<< "7 * diag\t" <<diag.mul(7)<<endl;
    	std::cout<< "7 * matrix\t" <<m1.mul(7)<<endl;
    	std::cout<<"vector * vector\t"<<m6*m7<<endl;
    	*/
    	/*std::cout<< m4.rank()<<std::endl;
    	std::cout<< m10.rank()<<std::endl;
    	std::cout<< m8.rank()<<std::endl;
    	std::cout<<m8<<std::endl;
    	std::cout<<m8.inverse()<<std::endl;
    	*/
    	/*std::ofstream OF("fileOUT.txt");
    	if(!OF) return -1;
    	OF<<m7;
    	std::ifstream IF("fileIN.txt");
    	if(!IF) return -1;
    	//OF>>m7;
    	//std::cin>>m;
    	std::cout<<m4;
    	OF<<m4;
    	IF>>m;
    	std::cout<<m;
    	*/
