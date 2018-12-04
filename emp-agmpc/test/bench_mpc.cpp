#include <emp-tool/emp-tool.h>
#include "emp-agmpc/emp-agmpc.h"
#include "../json/json.hpp"
#include "../Eigen/Dense"

using namespace std;
using namespace emp;
using json = nlohmann::json;
using namespace Eigen;

const string circuit_file_location = macro_xstr(EMP_CIRCUIT_PATH);

const static int nP = 2;
int party, port;
int rows = 100;
int cols = 10;
double rho = 0.01;
int VALUE_LENGTH = 16;
int EXPONENT_LENGTH = 16;

//Convert double input into bit representation with 16 bit decimal precision and 16 bit value 
//0th element is exponent, 1st element is value
vector<short> convert_double(int value_length, int expnt_length, double input) {
    double abs = std::abs(input);
    double lo = pow(2, value_length - 2);
    double up = pow(2, value_length - 1);
    short p = 0;
    while (abs > 0. && abs < lo) {abs *=2; --p;}
    while (abs >= up) {abs /=2; ++p;}
    vector<short> result(2);
    result[0] = p;
    short val;
    val = (abs * (input > 0 ? 1: -1));
    result[1] = val;
    return result; 
}


vector<vector<short>*> readMatrix(string file_name, double rho, int rows, int cols) {
    cout << "Reading matrix" << endl;
    std::ifstream i("../"+file_name);
    json j;
    i >> j;
    MatrixXd data_matrix(rows, cols);
    VectorXd y(rows);
    vector<double> x_data = j["x"];
    vector<double> y_data = j["y"];
    for (int i = 0; i < rows; i++) {
        y[i] = y_data[i];
        for (int j = 0; j < cols; j++) {
            data_matrix(i, j) = x_data[i * cols + j];
        }
    }


    //Compute values of matrices
    MatrixXd transpose = data_matrix.transpose();

    
    MatrixXd XTX = transpose * data_matrix;
    MatrixXd identity = MatrixXd::Identity(cols, cols);
    MatrixXd rho_identity = rho * identity;
    MatrixXd XTX_rhoI = XTX + rho_identity;
    MatrixXd inverse = XTX_rhoI.inverse();
    MatrixXd XTy = transpose * y;

    cout << "XTX + rho I inverse" << endl;
    cout << inverse << endl;
    cout << "XTy" << endl;
    cout << XTy << endl;
    
    vector<short> * inverse_float = new vector<short>[cols * cols];
  
    int r = inverse.rows();
    int c = inverse.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            inverse_float[i * c + j] = convert_double(VALUE_LENGTH, EXPONENT_LENGTH, inverse(i, j));
        }
    }  

    cout << "Got XTX inverse" << endl;

    vector<short> * XTy_float = new vector<short>[cols];
    r = XTy.rows();
    c = XTy.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            XTy_float[i * c + j] = convert_double(VALUE_LENGTH, EXPONENT_LENGTH, XTy(i, j));
        }
    }
    cout << "Got XTy " << endl;
    
    
    vector<vector<short> *> values(2);
    values[0] = inverse_float;
    values[1] = XTy_float;
    return values;
}
void bench_once(NetIOMP<nP> * ios[2], ThreadPool * pool, string filename) {
	cout << "Party " << party << endl;
	if(party == 1)cout <<"CIRCUIT:\t"<<filename<<endl;
	//string file = circuit_file_location+"/"+filename;
	CircuitFile cf(filename.c_str());

	auto start = clock_start();
	CMPC<nP>* mpc = new CMPC<nP>(ios, pool, party, &cf);
	ios[0]->flush();
	ios[1]->flush();
	double t2 = time_from(start);
//	ios[0]->sync();
//	ios[1]->sync();
	if(party == 1)cout <<"Setup:\t"<<party<<"\t"<< t2 <<"\n"<<flush;

	start = clock_start();
	mpc->function_independent();
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	if(party == 1)cout <<"FUNC_IND:\t"<<party<<"\t"<<t2<<" \n"<<flush;
//	cout << "Input size " << cf.n1 + cf.n2 << endl;
	start = clock_start();
	mpc->function_dependent();
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	if(party == 1)cout <<"FUNC_DEP:\t"<<party<<"\t"<<t2<<" \n"<<flush;


//	bool *in = new bool[cf.n1+cf.n2]; bool *out = new bool[cf.n3];
	cout << "Output size " << cf.n3 << endl;
	string file_name = "data" + to_string(party) + ".json";
//	vector<Float*> values = readMatrix(file_name, rho, rows, cols);
        
	cout << "Read matrix" << endl;
	vector<vector<short>*> values = readMatrix(file_name, rho, rows, cols);
	vector<short>*  inverse_float  = values[0];
	vector<short>*  XTy_float = values[1];
 	
	// Convert json file to bits
	bool* in = new bool[cf.n1 + cf.n2]; bool* out = new bool[cf.n3];
	int index = 0;
	cout << "Looking at the bits of XTX" << endl;
	for (int i = 0; i < cols * cols; i++) {
		vector<short> inverse_i  = inverse_float[i];
                short expnt = inverse_i[0];
		short val = inverse_i[1];
		for (int i = VALUE_LENGTH - 1; i >= 0; i--) {
			in[index++] = (val >> i) & 1;
		}
		for (int i = EXPONENT_LENGTH - 1; i >= 0; i--) {
			in[index++] = (expnt >> i) & 1;	
		}
	
	} 
	for (int i = 0; i < cols; i++) {
		vector<short> XTy_i = XTy_float[i];
		short expnt = XTy_i[0];
		short val = XTy_i[1];
               	for (int i = VALUE_LENGTH - 1; i >= 0; i--) {
			in[index++] = (val >> i) & 1;
		}
		 for (int i = EXPONENT_LENGTH - 1; i >= 0; i--) {
			in[index++] = (expnt >> i) & 1;	
		}
	}
      

	
	// Test to see if the inputs actually work
	cout << "Index " << index << endl;
	memset(in, false, cf.n1+cf.n2);
//	string res = "";
	int check_index = 0;
	for (int i = 0; i < cols * cols; i++) {
		
		short val = 0;
		for (int j = 0; j < VALUE_LENGTH; j++) {
			int x = in[check_index++] ? 1 : 0;
			val = (val << 1) | x;
		}
		short expnt  = 0;
		for (int j = 0; j < EXPONENT_LENGTH; j++) {
			int x = in[check_index++] ? 1 : 0;
			expnt = (expnt << 1) | x;
		}
		cout << "Exponent: " << expnt << "Value: " << val << " Decimal: " << (double) (val *  pow(2, expnt)) << endl;
	}
	start = clock_start();
	mpc->online(in, out);
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	// Turn output into actual floats
//	uint64_t band2 = io.count();
//	if(party == 1)cout <<"bandwidth\t"<<party<<"\t"<<band2<<endl;
	if(party == 1)cout <<"ONLINE:\t"<<party<<"\t"<<t2<<" \n"<<flush;
	index = cf.n3 - (cols * (EXPONENT_LENGTH + VALUE_LENGTH));
	cout << "Starting at output index " << index << endl;
	vector<double> result(cols);
	for (int i = 0; i < cols; i++) {
		short expnt = 0;
		for (int j = 0; j < EXPONENT_LENGTH; j++) {
			expnt = expnt | out[index++];
			expnt = expnt << 1;
		}
		short val = 0;
		for (int k = 0; k < VALUE_LENGTH; k++) {
			val = val | out[index++];
			val = val << 1;
		}
	
		result[i] = (double) (val *  pow(2, expnt));
	}
	cout << "Printing Weights" << endl;
	for (int i = 0; i < cols; i++) {
		cout << result[i] << endl;
	}

	
		string res = "";
		for (int i = 0; i < cf.n3; i++) {
			res += (out[i] ? "1":"0");
		}
		cout << "Result " << res << endl;
	
	
	delete mpc;
}
int main(int argc, char** argv) {
	int func = 3;
	parse_party_and_port(argv, &party, &port);
	cout << "Party " << party << "Port " << port << endl;
	if(party > nP)return 0;

	NetIOMP<nP> io(party, port);
	cout << "Did at least something work" << endl;
#ifdef LOCALHOST
	NetIOMP<nP> io2(party, port+2*(nP+1)*(nP+1)+1);
#else
	NetIOMP<nP> io2(party, port+2*(nP+1));
#endif
	NetIOMP<nP> *ios[2] = {&io, &io2};
	ThreadPool pool(2*(nP-1)+2);	
	cout << "HELLO" << endl;
//	for(int i = 0; i < 10; ++i)	
	if(func == 0)	
	bench_once(ios, &pool, circuit_file_location+"AES-non-expanded.txt");
//	for(int i = 0; i < 10; ++i)
	if(func == 1)	
	bench_once(ios, &pool, circuit_file_location+"sha-1.txt");
//	for(int i = 0; i < 10; ++i)
	if(func == 2)	
	bench_once(ios, &pool, circuit_file_location+"sha-256.txt");
	if(func == 3)
	bench_once(ios, &pool, circuit_file_location+"sort.txt");	
//	bench_once(ios, &pool, "/home/wangxiao/git/emp-toolkit/constantmpc/circ.txt");
	return 0;
}
