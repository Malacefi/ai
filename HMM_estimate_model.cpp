// Emelie Eriksson, Pontus Brink
#include <iostream>
#include <vector>
#include <stdlib.h>
using namespace std;

class Matrix {
  public:
  	int rows, cols;
	long double ** matrix;
	Matrix(int rows, int cols){
		this->rows= rows;
		this->cols = cols;
		this->matrix = new long double*[rows];
		for (int i = 0; i<rows; ++i){
			this->matrix[i] = new long double[cols];
		}
	}
	Matrix(const Matrix& m){
		this->rows = m.rows;
		this->cols = m.cols;
		this->matrix = new long double*[rows];
		for (int i =0; i<this->rows; ++i){
			this->matrix[i] = new long double[cols];
			for (int j = 0; j<this->cols; ++j){
				this->matrix[i][j] = m.matrix[i][j];
			}
		}
	}

	long double* operator[](int col){
			return matrix[col];
	}

	Matrix operator=(const Matrix& m2){
		for (int i =0; i<this->rows; ++i){
			for (int j = 0; j<this->cols; ++j){
				long double value = m2.matrix[i][j];
				this->matrix[i][j] = value;
			}
		}
		return *this;
	}

	Matrix operator-(Matrix m2){
		Matrix result(this->rows, this->cols);
		for (int i =0; i<this->rows; ++i){
			for (int j = 0; j<this->cols; ++j){
				if (this->matrix[i][j] > m2[i][j]){
					result[i][j] = this->matrix[i][j] - m2[i][j];
				}
				else{
					result[i][j] = m2[i][j]-this->matrix[i][j];
				}
				//cout << result[i][j] << " ";
			}
			//cout << "\n";
		}
		return result;
	}

	long double sum(void){
		long double sum = 0.0;
		for (int i = 0; i< this->rows; ++i){
			for (int j = 0; j<this->cols; ++j){
				sum += this->matrix[i][j];
				//cout << sum << "\n";
			}
		}
		return sum;
	}

	Matrix operator*(Matrix m2){
		if (this->cols != m2.rows){
			throw invalid_argument("Invalid matrix dimensions\n");
		}
		Matrix result(this->rows, m2.cols);
		long double temp;
		for (int i = 0; i<this->rows; ++i){
			for (int j = 0; j<m2.cols; ++j){
				// Now we are at a element in result matrix
				temp = 0.0;
				for (int k = 0; k<this->cols; ++k){
					temp += this->matrix[i][k]*m2[k][j];
				}
				result[i][j] = temp;
			}
		}
		return result;
	}

	void print(void){
		cout << this->rows << " " << this->cols << " ";
		for (int i = 0; i<this->rows; ++i){
			for (int j = 0; j<this->cols; ++j){
				cout << this->matrix[i][j] << " ";
			}
		
		}
		cout << "\n";
	}

};

Matrix alpha(Matrix A, Matrix B, Matrix pi, Matrix omissionSequence){
	Matrix alpha(omissionSequence.cols, A.rows);
	for (int i = 0; i<A.rows; ++i){
		alpha[0][i] = B[i][int(omissionSequence[0][0])] * pi[0][i];
	}

	for(int t = 1; t<omissionSequence.cols; ++t){
		for (int i = 0; i<A.rows;++i){
			long double sum = 0;
			for (int j = 0; j<A.rows; ++j){
				sum += A[j][i] * alpha[t-1][j];
			}
			alpha[t][i] = B[i][int(omissionSequence[0][t])] * sum;
		}
	}
	return alpha;
}


Matrix beta(Matrix A, Matrix B, Matrix omissionSequence){
	Matrix beta(omissionSequence.cols, A.rows);
	for (int i = 0; i<A.rows; ++i){
		beta[omissionSequence.cols-1][i] = 1;
	}
	for (int t = omissionSequence.cols-2; t>=0; --t){
		for (int i = 0; i<A.rows; ++i){
			long double sum = 0;
			for (int j = 0; j<A.rows; ++j){
				sum += beta[t+1][j] * B[j][int(omissionSequence[0][t+1])] * A[i][j];
			}
			beta[t][i] = sum;
		}
	}
	return beta;
}

void digamma(Matrix A, Matrix B, Matrix alpha, Matrix beta, Matrix omissionSequence, vector<Matrix>& digamma){
	for (int t = 0; t<omissionSequence.cols-1; ++t){
		Matrix toAdd(A.rows, A.cols);
		//cout << "running digamma: " << t << "\n";
		for (int i = 0; i<A.rows; ++i){
			for (int j = 0; j<A.rows; ++j){
				long double numerator = alpha[t][i] * A[i][j] * B[j][int(omissionSequence[0][t+1])] * beta[t+1][j];
				//cout << "Crashed at numerator calc\n";
				long double denominator = 0;
				for (int k = 0; k<A.rows; ++k){
					denominator += alpha[omissionSequence.cols-1][k];
				}
				toAdd[i][j] = numerator/denominator;
			}
		}
		//cout << "DIGAMMA DUN\n";
		digamma.push_back(toAdd);
		//cout << "digamma pushed\n";
	}
}


Matrix gamma(int omissionLength, vector<Matrix> &digamma){
	Matrix gamma(omissionLength-1, digamma.at(0).rows);
	for (int t = 0; t<gamma.rows; ++t){
		for (int i = 0; i<gamma.cols; ++i){
			long double sum = 0;
			for (int j = 0; j<gamma.cols; ++j){
				sum += digamma.at(t)[i][j];
			}
			gamma[t][i] = sum;
		}
	}
	return gamma;
}

void updateLambda(Matrix &A, Matrix &B, Matrix &pi, Matrix omissionSequence, vector<Matrix> &digamma, Matrix gamma){
	for (int i = 0; i<A.rows; ++i){
		for (int j = 0; j<A.cols; ++j){
			long double numerator = 0;
			long double denominator = 0;
			for(int t = 0; t<omissionSequence.cols-1; ++t){
				numerator += digamma.at(t)[i][j];
				denominator += gamma[t][i];
			}
			A[i][j] = numerator / denominator;
		}
	}
	for(int j = 0; j<A.rows; ++j){
		for(int k = 0; k<B.cols; ++k){
			long double numerator = 0;
			long double denominator = 0;
			for(int t = 0; t<omissionSequence.cols-1; ++t){
				if (omissionSequence[0][t] == k){
					numerator += gamma[t][j];
				}
				denominator += gamma[t][j];
			}
			B[j][k] = numerator/denominator;
		}
	}
	for (int i = 0; i<pi.cols; ++i){
		pi[0][i] = gamma[0][i];
	}

}

bool converges(Matrix A, Matrix oldA, Matrix B, Matrix oldB, Matrix pi, Matrix oldpi){
	long double threshold = 0.001;
	//cout << "DIFFERENCE: " << (A - oldA).sum()  << "\n";

	if ((A - oldA).sum() > threshold){
		return true;
	}
	return false;
}

Matrix viterbi(Matrix A, Matrix B, Matrix pi, Matrix omissionSequence){
	Matrix t1(A.rows, omissionSequence.cols);
	Matrix t2(A.rows, omissionSequence.cols);
	for (int i = 0; i<A.rows; ++i){
		t1[i][0] = pi[0][i] * B[i][int(omissionSequence[0][0])];
		t2[i][0] = 0;
	}
	long double maxTxA = 0;
	int maxState = 0;
	for (int i = 1; i<omissionSequence.cols; ++i){
		for (int j = 0; j<A.rows; ++j){
			maxTxA = 0;
			maxState = 0;
			for (int k = 0; k<A.rows; ++k){
				if (maxTxA < t1[k][i-1] * A[k][j]){
					maxTxA = t1[k][i-1] * A[k][j];
					maxState = k;
				}
			}
			t1[j][i] = B[j][int(omissionSequence[0][i])] * maxTxA;
			t2[j][i] = maxState;
		}
	}
	Matrix X(1, omissionSequence.cols);
	long double maxInT = 0;
	int index = 0;
	for (int i = 0; i<A.rows; ++i){
		if (maxInT < t1[i][omissionSequence.cols-1]){
			maxInT = t1[i][omissionSequence.cols-1];
			index = i;
		}
	}
	X[0][omissionSequence.cols-1] = index;
	for (int i = omissionSequence.cols-1; i>0; --i){
		X[0][i-1] = t2[X[0][i]][i];
	}

	return X;
}

int main(){

	int maRows, maCols, mbRows, mbCols;
	// Init ma.
	cin >> maRows >> maCols;
	Matrix ma(maRows, maCols);
	for (int i = 0; i<maRows; ++i){
		for (int j = 0; j<maCols; ++j){
			cin >> ma[i][j]; 
		}
	}

	// Init mb.
	
	cin >> mbRows >> mbCols;
	Matrix mb(mbRows, mbCols);
	for (int i = 0; i<mbRows; ++i){
		for (int j = 0; j<mbCols; ++j){
			cin >> mb[i][j]; 
		}
	}

	// Init pi
	int piRows, piCols;
	cin >> piRows >> piCols;
	Matrix pi(piRows, piCols);
	//cout << piRows << " " << piCols << "\n";
	for (int i = 0; i<piRows; ++i){
		//cout << "Survived row\n";
		for (int j = 0; j<piCols; ++j){
			cin >> pi[i][j]; 
		}
	}
	//cout << "Pi survived\n";
	int M;
	cin >> M;
	Matrix omissionSequence(1,M);
	for (int i = 0; i<M; ++i){
		cin >> omissionSequence[0][i];
	}



	
	Matrix al = alpha(ma, mb, pi, omissionSequence);
	Matrix be = beta(ma, mb, omissionSequence);
	vector<Matrix> dig;
	digamma(ma, mb, al, be, omissionSequence, dig);
	Matrix ga = gamma(omissionSequence.cols, dig);
	
	Matrix oldMa = ma;
	Matrix oldMb = mb;
	Matrix oldPi = pi;
	
	updateLambda(ma, mb, pi, omissionSequence, dig, ga);
	
	int max = 0;
	bool keepGoing = converges(ma, oldMa, mb, oldMb, pi, oldPi);
	while(keepGoing){
		dig.clear();
		al = alpha(ma, mb, pi, omissionSequence);
		be = beta(ma, mb, omissionSequence);
		digamma(ma, mb, al, be, omissionSequence, dig);
		ga = gamma(omissionSequence.cols, dig);
		oldMa = ma;
		oldMb = mb;
		oldPi = pi;
		//ma.print();
		updateLambda(ma, mb, pi, omissionSequence, dig, ga);
		//ma.print();
		++max;
		keepGoing = converges(ma, oldMa, mb, oldMb, pi, oldPi);
		

	}
	ma.print();
	mb.print();
	//cout << max << "\n";
}