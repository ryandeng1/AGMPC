#include "emp-tool/emp-tool.h"
using namespace emp;
using namespace std;
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


int main(int argc, char** argv) {
//	setup_plain_prot(true, "sort.txt");
//	sort(128);	
//	mult(2048);
//	ham(1<<10);
//	finalize_plain_prot ();
	int dim = 3;
	Float left_matrix[dim*dim];
	Float right_matrix[dim*dim];
	for (int i = 0; i < dim * dim; i++) {
		left_matrix[i] = Float(40, 20, i);
		right_matrix[i] = Float(40, 20, i);	
	}
	
	Float* result  = matrix_mul(left_matrix, right_matrix, dim);
	for (int i = 0; i < dim * dim; i++) {
		result[i].reveal<double>();
	}
	
	finalize_plain_prot();
}	
