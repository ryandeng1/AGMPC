#include "emp-tool/emp-tool.h"
#include "json/json.hpp"
#include "Eigen/Dense"
using namespace emp;
using namespace std;
using json = nlohmann::json;
using namespace Eigen;
const int DIMENSION = 2;
const int NUM_EXPNT = 16;
const int NUM_VALUE = 16;




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
	Integer c = Integer(2 * n, 0, PUBLIC); //a*b;
	c = a * b;
	c.reveal<string>();
}
void modexp(int n1, int n2) {
	Integer a(n1, 0,  ALICE);
	Integer b(n2, 0,  BOB);
	Integer c(n1, 5, ALICE);
	Integer d = a.modExp(b, c);
}


void add(int n) {
	Integer a(n, 0, ALICE);
	Integer b(n, 0, BOB);
	Integer c = a + b;
	c.reveal<string>();
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



Float*  matrix_mul(Float*  left_matrix, Float* right_matrix, int left_rows, int right_cols, int left_cols) {
	Float* result = new Float[left_rows * right_cols];
	for (int h = 0; h < left_rows * right_cols; h++) {
		result[h] = Float(NUM_VALUE, NUM_EXPNT, 0);
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
    Bit neg_th_greater_v = (Float(NUM_VALUE, NUM_EXPNT, -1) * th).greater(v); 
    if (v_greater_th.reveal(PUBLIC)) {
        return v - th;
    }
    else if (neg_th_greater_v.reveal(PUBLIC)) {
        return v + th;
    }
    else {
        return Float(NUM_VALUE, NUM_EXPNT, 0);
    }
}

Float* soft_threshold_vec(Float th, Float* vec, int rows, int cols) {
	Float* result = new Float[rows * cols];
	for (int i = 0; i < rows * cols; i++) {
   		result[i] = soft_threshold(th, vec[i]);
   	}
    
    return result;
}
/*
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

*/
Float* admm_local(Float* XXinv, Float* XTy, Float* u, Float* z, Float rho, Float l) {
    int dim = DIMENSION;
    Float* neg_u = new Float[dim];
    for (int i = 0; i < dim; i++) {
	neg_u[i] = Float(NUM_VALUE, NUM_EXPNT, -1) * u[i];
    }
    Float* z_u = add_matrix(z, neg_u, dim, 1);
    for (int i = 0; i < dim; i++) {
	z_u[i] = rho * z_u[i];
    }
    Float* second_term = add_matrix(XTy, z_u, dim, 1);
    
    Float* w = matrix_mul(XXinv, second_term, dim, 1, dim);
    /*
    cout << "ADMM LOCAL WEIGHT " << endl;
    for (int i = 0; i < dim; i++) {
	cout << w[i].reveal<string>() << endl;
    }
    */
    return w;

}

vector<vector<Float*>> admm_coordinate(vector<Float*> w_list, vector<Float*> u_list, Float* z, Float rho, Float l, int nparties) {
	int dim = DIMENSION;
	Float* w_avg = new Float[dim];
	Float* u_avg = new Float[dim];
	for (int i = 0; i < dim; i++) {
		w_avg[i] = Float(NUM_VALUE, NUM_EXPNT, 0);
		u_avg[i] = Float(NUM_VALUE, NUM_EXPNT, 0);
	}
	for (int i = 0; i < nparties; i++) {
		w_avg = add_matrix(w_avg, w_list[i], dim, 1);
		u_avg = add_matrix(u_avg, u_list[i], dim, 1);
	}
	double nparties_d = 1.0 / nparties;
	for (int i = 0; i < dim; i++) {
		w_avg[i] = w_avg[i] * Float(NUM_VALUE, NUM_EXPNT, nparties_d);
		u_avg[i] = u_avg[i] * Float(NUM_VALUE, NUM_EXPNT, nparties_d);
	}
	Float th = l / (rho * Float(NUM_VALUE, NUM_EXPNT, nparties));
	Float* w_u_avg_sum = add_matrix(w_avg, u_avg, dim, 1);
	Float* z_new = soft_threshold_vec(th, w_u_avg_sum, dim, 1);
	Float* z_new_neg = new Float[dim];
	for (int i = 0; i < dim; i++) {
		z_new_neg[i] = Float(NUM_VALUE, NUM_EXPNT, -1) * z_new[i];
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
			w_list[i][j] = Float(NUM_VALUE, NUM_EXPNT, 0);
			u_list[i][j] = Float(NUM_VALUE, NUM_EXPNT, 0);
		}
	}

	for (int i = 0; i < dim; i++) {
		z[i] = Float(NUM_VALUE, NUM_EXPNT, 0);
	}
	//cout << "Inited everything" << endl;
	for (int i = 0; i < admm_iter; i++) {
		for (int j = 0; j < nparties; j++) {
			//cout << "Did I make it here?" << endl;
			w_list[j] = admm_local(XXinv_cache[j], XTy_cache[j], u_list[j], z, rho, l);
		}

		vector<vector<Float*>> vals = admm_coordinate(w_list, u_list, z, rho, l, nparties);
		u_list = vals[0];
		z = vals[1][0];
	}

	return z;
}



void main_func() {
	int cols = DIMENSION;
	int nparties = 2;
	int admm_iter = 10;
	Float rho(NUM_VALUE, NUM_EXPNT, 0.01);
	Float l(NUM_VALUE, NUM_EXPNT, 0.008);
	vector<Float*> XXinv_cache(nparties);
	vector<Float*> XTy_cache(nparties);
		
	for (int i = 0; i < nparties; i++) {
		Float* XXinv = new Float[cols * cols];
		for (int j = 0; j < cols * cols; j++) {
			if (i == 0) {
				XXinv[j] = Float(NUM_VALUE, NUM_EXPNT, 0, ALICE);
			}
			else {
				XXinv[j] = Float(NUM_VALUE, NUM_EXPNT, 0, BOB);
			}
		} 
		Float* XTy = new Float[cols];
		for (int h = 0; h < cols; h++) {
			if (i == 0) {
				XTy[h] = Float(NUM_VALUE, NUM_EXPNT, 0, ALICE);
			}
			else {
				XTy[h] = Float(NUM_VALUE, NUM_EXPNT, 0, BOB);
			}
		}
		XXinv_cache[i] = XXinv;
		XTy_cache[i] = XTy;
	}
	
	/*
	Float* inverse_one = new Float[cols * cols];
	for (int i = 0; i < cols * cols; i++) {
		inverse_one[i] = Float(40, 20, 0, ALICE);
	}

	Float* XTy_one = new Float[cols];
	for (int i = 0; i < cols; i++) {
		XTy_one[i] = Float(40, 20, 0, ALICE);
	}
	Float* inverse_two = new Float[cols * cols];
	for (int i = 0; i < cols * cols; i++) {
		inverse_two[i] = Float(40, 20, 0, BOB);
	}
	Float* XTy_two = new Float[cols];
	for (int i = 0; i < cols; i++) {
		XTy_two[i] = Float(40, 20, 0, BOB);
	}
	XXinv_cache[0] = inverse_one; XXinv_cache[1] = inverse_two;
	XTy_cache[0] = XTy_one; XTy_cache[1] = XTy_two; 
	*/
	Float* z = admm(XXinv_cache, XTy_cache, admm_iter, rho, l, nparties);
	
	for (int i = 0; i < DIMENSION; i++) {
		cout << z[i].reveal<string>() << endl;
	}
		
		
}



vector<Integer*> readMatrix(string file_name, int length) {
    cout << "Reading matrix" << endl;
    std::ifstream i(file_name);
    json j;
    i >> j;
    int rows = DIMENSION;
    int cols = DIMENSION;
    Integer* matrix1 = new Integer[rows * cols];
    Integer* matrix2 = new Integer[rows * cols];
    vector<double> x_data = j["x"];
    vector<double> y_data = j["y"];

    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
	    matrix1[i * cols + j] = Integer(length, x_data[i * cols + j], ALICE);
        }
    }
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
	    matrix2[i * cols + j] = Integer(length, y_data[i * cols + j], BOB);
        }
    }

    vector<Integer*> values(2);
    values[0] = matrix1;
    values[1] = matrix2;
    return values;
}


Integer* mat_mul(Integer* left, Integer* right, int left_row, int left_col, int right_row, int right_col) {
    Integer* result = new Integer[left_row * right_col];
    for (int i = 0; i < left_row * right_col; i++) {
	result[i] = Integer(64, 0, PUBLIC);
    }    
    for (int i = 0; i < left_row; i++) {
	for (int j = 0; j < right_col; j++) {
	    for (int k = 0; k < left_col; k++) {
		result[i * right_col + j] = result[i * right_col + j] +  left[i * left_col + k] * right[k * right_col + j];
	    }
	}
    }


    return result;
}





int main(int argc, char** argv) {
	setup_plain_prot(true, "sort.txt");
/*
	//string file_name = "data.json";
	//cout << "Did I make it here?" << endl;

	*/
	/*
        string file_name = "data.json";
        int len = 32;
	
	int dim = DIMENSION;
	Integer* left = new Integer[dim * dim];//matrices[0];
	Integer* right = new Integer[dim * dim]; //matrices[1];
	for (int i = 0; i < dim; i++) {
		for (int j = 0; j < dim; j++) {
			left[i * dim + j] = Integer(len, 0, ALICE);
			right[i * dim + j] = Integer(len, 0, BOB);	
		}
	}


	Integer* result = mat_mul(left, right, dim, dim, dim, dim);
	cout << "RESULT" << endl;
	for (int i = 0; i < dim; i++) {
	    for (int j = 0; j < dim; j++) {
		Integer x = result[i * dim + j];
		x.reveal<string>();
	    }
	}
	*/
	add(32);
	finalize_plain_prot();
}	
