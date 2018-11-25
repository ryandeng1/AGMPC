#include "emp-tool/emp-tool.h"
#include "json/json.hpp"
#include "Eigen/Dense"
using namespace emp;
using namespace std;
using json = nlohmann::json;
using namespace Eigen;
const int NUM_BITS = 20;
const int NUM_DECIMAL = 40;
const int DIMENSION = 10;



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


Float*  matrix_mul(Float*  left_matrix, Float* right_matrix, int left_rows, int right_cols, int left_cols) {
	Float* result = new Float[left_rows * right_cols];
	for (int h = 0; h < left_rows * right_cols; h++) {
		result[h] = Float(40, 20, 0);
	}
	for (int i = 0; i < left_rows; i++) {
		for (int j = 0; j < right_cols; j++) {
			for (int k = 0; k < left_cols; k++) {
				result[i * right_cols + j] = result[i * right_cols + j] + left_matrix[i * left_cols + k] * right_matrix[k * right_cols + j];
			}
		}
	}
	return result;
} 

Float* add_matrix(Float* left_matrix, Float* right_matrix, int rows, int cols) {
	Float* result = new Float[rows * cols];
	for (int i = 0; i < rows; i++) {
		for (int j = 0; j < cols; j++) {
			result[i * cols + j]  = left_matrix[i * cols + j] + right_matrix[i * cols + j];
		}
	}
	return result;
}



Float soft_threshold(Float th, Float v) {
    Bit v_greater_th = v.greater(th);
    Bit neg_th_greater_v = (Float(40, 20, -1) * th).greater(v); 
    if (v_greater_th.reveal(PUBLIC)) {
        return v - th;
    }
    else if (neg_th_greater_v.reveal(PUBLIC)) {
        return v + th;
    }
    else {
        return Float(40, 20, 0);
    }
}

Float* soft_threshold_vec(Float th, Float* vec, int rows, int cols) {
	Float* result = new Float[rows * cols];
	for (int i = 0; i < rows * cols; i++) {
   		result[i] = soft_threshold(th, vec[i]);
   	}
    
    return result;
}

vector<Float*> readMatrix(string file_name, double rho, int rows, int cols) {
    cout << "Reading matrix" << endl;
    std::ifstream i(file_name);
    json j;
    i >> j;
    int dim = DIMENSION;
    MatrixXd data_matrix(rows, dim);
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
    MatrixXd identity = MatrixXd::Identity(dim, dim);
    MatrixXd rho_identity = rho * identity;
    MatrixXd XTX_rhoI = XTX + rho_identity;
    MatrixXd inverse = XTX_rhoI.inverse();
    MatrixXd XTy = transpose * y;

    cout << "XTX + rho I inverse" << endl;
    cout << inverse << endl;
    cout << "XTy" << endl;
    cout << XTy << endl;

    Float* inverse_float = new Float[cols * cols];
    int r = inverse.rows();
    int c = inverse.cols();
    for (int i = 0; i < r; i++) {
        for (int j = 0; j < c; j++) {
            inverse_float[i * c + j] = Float(40, 20, inverse(i, j));
         
        }
    }  



    Float* XTy_float = new Float[cols];
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
    int dim = DIMENSION;
    Float* neg_u = new Float[dim];
    for (int i = 0; i < dim; i++) {
	neg_u[i] = Float(40, 20, -1) * u[i];
    }
    Float* z_u = add_matrix(z, neg_u, dim, 1);
    for (int i = 0; i < dim; i++) {
	z_u[i] = rho * z_u[i];
    }
    Float* second_term = add_matrix(XTy, z_u, dim, 1);
    
    Float* w = matrix_mul(XXinv, second_term, dim, 1, dim);
    cout << "ADMM LOCAL WEIGHT " << endl;
    for (int i = 0; i < dim; i++) {
	cout << w[i].reveal<string>() << endl;
    }
    return w;

}

vector<vector<Float*>> admm_coordinate(vector<Float*> w_list, vector<Float*> u_list, Float* z, Float rho, Float l, int nparties) {
	cout << "Nparties " << nparties << endl;
	int dim = DIMENSION;
	Float* w_avg = new Float[dim];
	Float* u_avg = new Float[dim];
	for (int i = 0; i < dim; i++) {
		w_avg[i] = Float(40, 20, 0);
		u_avg[i] = Float(40, 20, 0);
	}
	for (int i = 0; i < nparties; i++) {
		w_avg = add_matrix(w_avg, w_list[i], dim, 1);
		u_avg = add_matrix(u_avg, u_list[i], dim, 1);
	}
	double nparties_d = 1.0 / nparties;
	for (int i = 0; i < dim; i++) {
		w_avg[i] = w_avg[i] * Float(40, 20, nparties_d);
		u_avg[i] = u_avg[i] * Float(40, 20, nparties_d);
	}
	Float th = l / (rho * Float(40, 20,nparties));
	Float* w_u_avg_sum = add_matrix(w_avg, u_avg, dim, 1);
	Float* z_new = soft_threshold_vec(th, w_u_avg_sum, dim, 1);
	Float* z_new_neg = new Float[dim];
	for (int i = 0; i < dim; i++) {
		z_new_neg[i] = Float(40, 20, -1) * z_new[i];
	}
	vector<Float*> new_ulist(nparties);
	for (int i = 0; i < nparties; i++) {
		Float* u = u_list[i];
		Float* w = w_list[i];
		Float* new_u = add_matrix(add_matrix(u, w, dim, 1), z_new_neg, dim, 1);
		new_ulist[i] = new_u;
	}

	vector<vector<Float*>> result(2);
	vector<Float*> z_new_vec(1);
	z_new_vec[0] = z_new;
	result[0] = new_ulist;
	result[1] = z_new_vec;
	/*
	cout << "Printing z at the end of admm coordinate" << endl;
	for (int i = 0; i < dim; i++) {
		cout << z_new[i].reveal<string>() << endl;
	}
	*/
	return result;
}




Float* admm(vector<Float*> XXinv_cache, vector<Float*> XTy_cache, int admm_iter, Float rho, Float l, int nparties) {
	vector<Float*> w_list(nparties);
	vector<Float*> u_list(nparties);
        int dim = DIMENSION;
	Float* z = new Float[dim];
	
	for (int i = 0; i < nparties; i++) {
		w_list[i] = new Float[dim];
		u_list[i] = new Float[dim];
		for (int j = 0; j < dim; j++) {
			w_list[i][j] = Float(40, 20, 0);
			u_list[i][j] = Float(40, 20, 0);
		}
	}

	for (int i = 0; i < dim; i++) {
		z[i] = Float(40, 20, 0);
	}
	cout << "Inited everything" << endl;
	for (int i = 0; i < admm_iter; i++) {
		for (int j = 0; j < nparties; j++) {
			cout << "Did I make it here?" << endl;
			w_list[j] = admm_local(XXinv_cache[j], XTy_cache[j], u_list[j], z, rho, l);
		}

		vector<vector<Float*>> vals = admm_coordinate(w_list, u_list, z, rho, l, nparties);
		u_list = vals[0];
		z = vals[1][0];
	}

	return z;
}













int main(int argc, char** argv) {
	setup_plain_prot(true, "sort.txt");
	int rows = 100;
	int cols = 10;
	int nparties = 2;
	int admm_iter = 10;
	Float rho(40, 20, 0.01);
	double rho_double = 0.01;
	Float l(40, 20, 0.008);
	vector<Float*> XXinv_cache(nparties);
	vector<Float*> XTy_cache(nparties);
	for (int i = 0; i < nparties; i++) {
		string file_name = "data" + to_string(i + 1) + ".json";
		vector<Float*> vals = readMatrix(file_name, rho_double, rows, cols);
		Float* XXinv = vals[0];
		Float* XTy = vals[1];
		XXinv_cache[i] = XXinv;
		XTy_cache[i] = XTy;
	}
	//string file_name = "data.json";
	//cout << "Did I make it here?" << endl;
	Float* z = admm(XXinv_cache, XTy_cache, admm_iter, rho, l, nparties);
	cout << "Printing weights" << endl;
	for (int i = 0; i < cols; i++) {
		cout << z[i].reveal<string>()  << endl;
	}



	
	finalize_plain_prot();
	return 0;
}	
