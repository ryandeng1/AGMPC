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
int cols = 5;
double rho = 0.01;
int VALUE_LENGTH = 20;
int EXPONENT_LENGTH = 20;

vector<Float*> readMatrix(string file_name, double rho, int rows, int cols) {
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
    //cout << "Matrix" << endl;
    //cout << data_matrix << endl;
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
    cout << "WHAT IS GOING ON" << endl;
    Float *inverse_float  = new Float[cols * cols];
    int r = cols;
    int c = cols;
    cout << "Rows " << r << "Columns " << c << endl;
  

    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            inverse_float[i * c + j] = Float(VALUE_LENGTH, EXPONENT_LENGTH, inverse(i, j));
         
        }
    }  

    cout << "Got XTX inverse" << endl;

    Float* XTy_float = new Float[cols];
    r = XTy.rows();
    c = XTy.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            XTy_float[i * c + j] = Float(VALUE_LENGTH, EXPONENT_LENGTH, XTy(i, j));
        }
    }
    cout << "Got XTy " << endl;

    vector<Float*> values(2);
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
	cout << "Hello" << endl;
	CMPC<nP>* mpc = new CMPC<nP>(ios, pool, party, &cf);
	ios[0]->flush();
	ios[1]->flush();
	cout << "Flushed" << endl;
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
	cout << "Did I complete function dependent?" << endl;
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);
	if(party == 1)cout <<"FUNC_DEP:\t"<<party<<"\t"<<t2<<" \n"<<flush;


	bool *in = new bool[cf.n1+cf.n2]; bool *out = new bool[cf.n3];
	cout << "Output size " << cf.n3 << endl;
	string file_name = "data" + to_string(party) + ".json";
	vector<Float*> values = readMatrix(file_name, rho, rows, cols);
	cout << "Read matrix" << endl;
	Float* XTX_rhoI = values[0];
	Float* XTy = values[1];
	int input_index = 0;
	cout << "HERE" << endl;
	for (int i = 0; i < cols * cols; i++) {
		Float float_val = XTX_rhoI[i];
		Integer val = float_val.value;
		Integer expnt = float_val.expnt;
		for (int j = 0; j < VALUE_LENGTH; j++) {
			in[input_index] = val[j].reveal<bool>();
			input_index++;
		}
		for (int j = 0; j < EXPONENT_LENGTH; j++) {
			in[input_index] = expnt[j].reveal<bool>();
			input_index++;
		}
	}
	for (int i = 0; i < cols; i++) {
		Float float_val = XTy[i];
		Integer val = float_val.value;
		Integer expnt = float_val.expnt;
		for (int j = 0; j < VALUE_LENGTH; j++) {
			in[input_index] = val[j].reveal<bool>();
			input_index++;
		}
		for (int j = 0; j < EXPONENT_LENGTH; j++) {
			in[input_index] = expnt[j].reveal<bool>();
			input_index++;
		}
	}





	//memset(in, false, cf.n1+cf.n2);
	start = clock_start();
	mpc->online(in, out);
	ios[0]->flush();
	ios[1]->flush();
	t2 = time_from(start);

//	uint64_t band2 = io.count();
//	if(party == 1)cout <<"bandwidth\t"<<party<<"\t"<<band2<<endl;
	if(party == 1)cout <<"ONLINE:\t"<<party<<"\t"<<t2<<" \n"<<flush;
	Float* result = new Float[cols];
	bool* pointer = out;
	for (int i = 0; i < cols; i++) {
		Integer value = Integer(VALUE_LENGTH, pointer);
		Integer expnt = Integer(EXPONENT_LENGTH, pointer + VALUE_LENGTH);
		Float float_val = Float(VALUE_LENGTH, EXPONENT_LENGTH, 0.0);
		float_val.value = value;
		float_val.expnt = expnt;
		result[i] = float_val;
		pointer = pointer + VALUE_LENGTH + EXPONENT_LENGTH;
	}
	

	cout << "Weights" << endl;
	for (int i = 0; i < cols; i++) {
		cout << result[i].reveal<string>() << endl;
	}
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
