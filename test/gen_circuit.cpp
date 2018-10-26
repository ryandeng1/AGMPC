#include "emp-tool/emp-tool.h"
#include "json/json.hpp"

using namespace emp;
using namespace std;

NUM_BITS = 20;
NUM_DECIMAL = 40;
NUM_ROWS = 10;
NUM_COLS = 10;



/*
void ham(int n) {
	Integer a(n, 0, ALICE);
	Integer b(n, 0, BOB);
	Integer c = a^b;
	Integer d = c.hamming_weight();
	d.reveal<string>();
}

void mult(int n) {
	Integer a(n, 0, ALICE);
	Integer b(n, 0, BOB);
	Integer c = a*b;
	c.reveal<string>();
}
void modexp(int n1, int n2) {
	Integer a(n1, 0,  ALICE);
	Integer b(n2, 0,  BOB);
	Integer c(n1, 5, ALICE);
	Integer d = a.modExp(b, c);
}
void sort(int n) {
	Integer *A = new Integer[n];
	Integer *B = new Integer[n];
	for(int i = 0; i < n; ++i)
		A[i] = Integer(32, 0, ALICE);
	for(int i = 0; i < n; ++i)
		B[i] = Integer(32, 0, BOB);
	for(int i = 0; i < n; ++i)
		A[i] = A[i] ^ B[i];
	sort(A, n);
	for(int i = 0; i < n; ++i)
		A[i].reveal<string>();

}
*/


Float*  matrix_mul(Float*  left_matrix, Float* right_matrix, int dim) {
	Float* result = new Float[dim * dim];
	for (int i = 0; i < dim * dim; i++) {
		result[i] = Float(40, 20, 0);
	}
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			for (int k = 0; k < dim; k++) {
				result[i * dim + j] = result[i * dim + j] + left_matrix[i * dim + k] * right_matrix[k * dim + j];
			}
		}
	}
	return result;
} 

Float* add_matrix(Float* left_matrix, Float* right_matrix, int dim) {
	Float* result = new Float[dim * dim];
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			result[i * dim + j]  = left_matrix[i * dim + j] + right_matrix[i * dim + j]
		}
	}
	return result;
}


Float soft_threshold(Float th, Float v) {
    if (v > th) {
        return v - th;
    }
    else if (v < -1 * th) {
        return v + th;
    }
    else {
        return 0;
    }
}

Float* soft_threshold_vec(Float th, Float* vec) {
	Float* result = new Float[dim * dim];
   	for (int i = 0; i < dim * dim; i++) {
   		result[i] = soft_threshold(th, vec[i]);
   	}
    
    return result;
}

vector<Float*> readMatrix(string file_name, double rho) {
    cout << "Reading matrix" << endl;
    cout << finish << endl;
    std::ifstream i(file_name);
    json j;
    i >> j;
    MatrixXd data_matrix(NUM_ROWS, NUM_COLUMNS);
    VectorXd y(NUM_ROWS);
    vector<double> x_data = j["x"];
    vector<double> y_data = j["y"];
    for (int i = 0; i < NUM_ROWS; i++) {
        y[i] = y_data[i];
        for (int j = 0; j < NUM_COLUMNS; j++) {
            data_matrix(i, j) = x_data[i * NUM_COLUMNS + j];
        }
    }


    //Compute values of matrices
    cout << "Matrix" << endl;
    cout << data_matrix << endl;
    MatrixXd transpose = data_matrix.transpose();

    
    MatrixXd XTX = transpose * data_matrix;
    MatrixXd identity = MatrixXd::Identity(NUM_ROWS, NUM_COLUMNS);
    MatrixXd rho_identity = rho * identity;
    MatrixXd XTX_rhoI = XTX + rho_identity;
    MatrixXd inverse = XTX_rhoI.inverse();
    MatrixXd XTy = transpose * y;

    Float* inverse_float = new Float[NUM_ROWS * NUM_COLS];
    int r = inverse.rows();
    int c = inverse.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            inverse_float[i * c + j] = Float(40, 20, inverse(i, j));
         
        }
    }  



    Float* XTy_float = new Float[NUM_ROWS];
    r = XTy.rows();
    c = XTy.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            XTy_float[i * c + j] = Float(40, 20, XTy(i, j));
        }
    }
    vector<Float*> values(2);
    values[0] = inverse_float;
    values[1] = XTy_float;
    return values;
}


Float* admm_local(Float* XXinv, Float* XTy, Float* u, Float* z, Float rho, Float l) {
	Float* neg_u = new Float[NUM_COLS];
	for (int i = 0; i < NUM_COLS; i++) {
		neg_u[i] = Float(40, 20, -1) * u[i];
	}

	Float* z_u = add_matrix(z, neg_u);
	for (int i = 0; i < NUM_COLS; i++) {
		z_u[i] = rho * z_u[i];
	}
	Float* second_term = add_matrix(XTy, z_u);
    Float* w = matrix_mul(XXinv, second_term)
    return w

}

vector<vector<Float*>> admm_coordinate(vector<Float*> w_list, vector<Float*> u_list, Float* z, Float rho, Float l, int nparties) {
	Float* w_avg = new Float[NUM_COLS];
	Float* u_avg = new Float[NUM_COLS];
	for (int i = 0; i < NUM_COLS; i++) {
		w_avg[i] = Float(40, 20, 0);
		u_avg[i] = Float(40, 20, 0);
	}
	for (int i = 0; i < nparties; i++) {
		w_avg = add_matrix(w_avg, w_list[i]);
		u_avg = add_matrix(w_avg, w_list[i]);
	}

	for (int i = 0; i < NUM_COLS; i++) {
		w_avg[i] = w_avg[i] / Float(40, 20, nparties);
		u_avg[i] = u_avg[i] / Float(40, 20, nparties);
	}

	Float th = l / (rho * nparties);
	Float* z_new = soft_threshold_vec(th, add_matrix(w_avg, u_avg));
	Float* z_new_neg = new Float[NUM_COLS];
	for (int i = 0; i < NUM_COLS; i++) {
		z_new_neg[i] = Float(40, 20, -1) * z_new;
	}
	vector<Float*> new_ulist(nparties);
	for (int i = 0; i < nparties; i++) {
		Float* u = u_list[i];
		Float* w = w_list[i];
		Float* new_u = add_matrix(add_matrix(u, w). z_new_neg);
		new_ulist[i] = new_u;
	}

	vector<vector<Float*>> result(2);
	vector<Float*> z_new_vec(1);
	z_new_vec[0] = z_new;
	result[0] = new_ulist;
	result[1] = z_new_vec;
	return result;
}




Float* admm(vector<Float*> XXinv_cache, vector<Float*> XTy_cache, int admm_iter, Float rho, Float l, int nparties) {
	vector<Float*> w_list(nparties);
	vector<Float*> u_list(nparties);
	Float* z = new Float[NUM_COLS];
	for (int i = 0; i < nparties; i++) {
		w_list[i] = new Float[NUM_ROWS];
		u_list[i] = new Float[NUM_ROWS];
		for (int j = 0; j < NUM_ROWS; j++) {
			w_list[i][j] = Float(40, 20, 0);
			u_list[i][j] = Float(40, 20, 0);
		}
	}

	for (int i = 0; i < NUM_COLS; i++) {
		z[i] = Float(40, 20, 0);
	}

	for (int i = 0; i < admm_iter; i++) {
		for (j = 0; j < nparties; j++) {
			w_list[j] = admm_local(XXinv_cache[j], XTy_cache[j], u_list[j], z, rho, l);
		}

		vector<vector<Float*>> vals = admm_coordinate(w_list, u_list, z, rho, l);
		u_list = vals[0];
		z = vals[1][0];
	}

	return z;
}













int main(int argc, char** argv) {
	setup_plain_prot(true, "sort.txt");
	int dim = 10;
	int npartiea = 4;
	int admm_iter = 10;
	Float rho(40, 20, 0.01);
	Float l(40, 20, 0.08);
	string file_name = "data.json";
	vector<Float*> vals = readMatrix(file_name, rho);
	Float* XXinv = vals[0];
	Float* XTy = vals[1];
	vector<Float*> XXinv_cache = {XXinv};
	vector<Float*> XTy_cache = {XTy};
	Float* z = admm(XXinv_cache, XTy_cache, admm_iter, rho, l, nparties);
	cout << "Printing weights" << endl;
	for (int i = 0; i < NUM_COLS; i++) {
		cout << z[i] << endl;
	}



	
	finalize_plain_prot();
	return 0;
}	
