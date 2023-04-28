#include <iostream>

using namespace std;

void matrix_plus_matrix(int n, double** a, double** a2) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = a[i][j] + a2[i][j];
        }
    }
}

void matrix_mult_matrix(int n, double** c, double** a, double** b) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            c[i][j] = 0;
            for (int k = 0; k < n; k++) {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

void matrix_mult_value(int n, double** a, double value) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = a[i][j] * value;
        }
    }
}

void matrix_mult_transp_vector(int n, double** a, double* u, double temp_v[]) {
    for (int i = 0; i < n; i++) {
        temp_v[i] = 0;
        for (int j = 0; j < n; j++) {
            temp_v[i] += a[i][j] * u[j];
        }
    }
}

void vector_mult_matrix(int n, double** a, double* v, double* v2) {
    for (int i = 0; i < n; i++) {
        v[i] = 0;
        for (int j = 0; j < n; j++) {
            v[i] += a[j][i] * v2[j];
        }
    }
}

double vector_mult_transp_vector(int n, double* v, double* u) {
    double value = 0;
    for (int i = 0; i < n; i++) {
        value += u[i] * v[i];
    }
    return value;
}

void transp_vector_mult_vector(int n, double* u, double* v, double** a) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = u[i] * v[j];
        }
    }
}

void transp_matrix(int n, double** a, double** transp_a) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            transp_a[i][j] = a[j][i];
        }
    }
}

void edging(int n, double** a, double** transp_a, double** temp_a, double** temp_a2, double* temp_ver) {
    double value = 0;

    vector_mult_matrix(n, a, temp_ver, a[n]);
    value = vector_mult_transp_vector(n, temp_ver, transp_a[n]);
    value = a[n][n] - value;
    a[n][n] = 1 / value; // правая нижняя часть в окаймленной матрице


    matrix_mult_transp_vector(n, a, transp_a[n], temp_ver);
    for (int i = 0; i < n; i++) {
        a[i][n] = -(temp_ver[i] / value); // верхняя правая чать в окаймленной матрице
    }

    transp_vector_mult_vector(n, temp_ver, a[n], temp_a);
    matrix_mult_matrix(n, temp_a2, temp_a, a);
    matrix_mult_value(n, temp_a2, (1 / value));

    vector_mult_matrix(n, a, temp_ver, a[n]);
    for (int i = 0; i < n; i++) {
        a[n][i] = -(temp_ver[i] / value); // нижняя левая часть в окаймленной матрице
    }

    matrix_plus_matrix(n, a, temp_a2); // верхняя левая часть в окаймленной матрице
}

int main()
{
    setlocale(LC_ALL, "Russian");
    srand(time(NULL));

    const int n = 7;

    double temp_ver[n]{ 0 };
    double b[n]{ 0 };
    double** a = new double* [n];
    double** a_copy = new double* [n];
    double** transp_a = new double* [n];
    double** temp_a = new double* [n];
    double** temp_a2 = new double* [n];

    double c[n][n] = { {0.411, 0.421, -0.333, 0.313, -0.141, -0.381, 0.245},
                     {0.241, 0.705, 0.139, -0.409, 0.321, 0.0625, 0.101},
                     {0.123, -0.239, 0.502, 0.901, 0.243, 0.819, 0.321},
                     {0.413, 0.309, 0.801, 0.865, 0.423, 0.118, 0.183},
                     {0.241, -0.221, -0.243, 0.134, 1.274, 0.712, 0.423},
                     {0.281, 0.525, 0.719, 0.118, -0.974, 0.808, 0.923 },
                     {0.246, -0.301, 0.231, 0.813, -0.702, 1.223, 1.105} };

    double d[n] = { 0.096, 1.252, 1.024, 1.023, 1.155, 1.937, 1.673 };

    for (int i = 0; i < n; i++) {
        a[i] = new double[n];
        a_copy[i] = new double[n];
        transp_a[i] = new double[n];
        temp_a[i] = new double[n];
        temp_a2[i] = new double[n];
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            a[i][j] = c[i][j];
            a_copy[i][j] = a[i][j];
        }
    }

    for (int i = 0; i < n; i++) {
        b[i] = d[i];
    }

    transp_matrix(n, a, transp_a);

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << a[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << "Свободные члены: ";

    for (int i = 0; i < n; i++) {
        cout << b[i] << " ";
    }

    cout << endl << endl;


    a[0][0] = 1 / a[0][0];

    for (int i = 1; i < n; i++) {
        edging(i, a, transp_a, temp_a, temp_a2, temp_ver);
        cout << "шаг" << i << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << a[i][j] << " ";
            }
            cout << endl;
        }
        cout << endl;
    }

    matrix_mult_transp_vector(n, a, b, temp_ver);

    matrix_mult_matrix(n, temp_a, a, a_copy);

    cout << "A*A^(-1): " << endl;

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            cout << temp_a[i][j] << " ";
        }
        cout << endl;
    }

    cout << endl << endl;

    cout << "Xk: ";

    for (int i = 0; i < n; i++) {
        cout << temp_ver[i] << " ";
    }

    cout << endl << endl << "Подставим Xk в вектор x и вычислим погрешность(Ax - b): ";

    for (int i = 0; i < n; i++) {
        double s = 0;
        for (int j = 0; j < n; j++) {
            s += temp_ver[j] * a_copy[i][j];
        }
        cout << abs(b[i] - s) << " ";
    }

    cout << endl;
}
