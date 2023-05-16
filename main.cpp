#include <ccomplex>
#include <cmath>
#include <iostream>
#include <vector>


std::complex<double> W(int n, int p) {
    return {cos(-2 * M_PI * p / n), sin(-2 * M_PI * p / n)};
}

void fft_square(int n, std::vector<std::vector<std::complex<double>>> &ar) {
    if (n < 3) {
        if (n == 1) {
            return;
        }
        auto f00 = ar[0][0] + ar[0][1] + ar[1][0] + ar[1][1];
        auto f01 = ar[0][0] - ar[0][1] + ar[1][0] - ar[1][1];
        auto f10 = ar[0][0] + ar[0][1] - ar[1][0] - ar[1][1];
        auto f11 = ar[0][0] - ar[0][1] - ar[1][0] + ar[1][1];
        ar[0][0] = f00; ar[0][1] = f01;
        ar[1][0] = f10; ar[1][1] = f11;
        return;
    }

    std::vector<std::vector<std::complex<double>>> tmp[4];
    for (auto &v : tmp) {
        v = std::vector<std::vector<std::complex<double>>>(n / 2,std::vector<std::complex<double>>(n / 2));
    }

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            tmp[(i & 1) + 2 * (j & 1)][i / 2][j / 2] = ar[i][j];  // 00 10 01 11
        }
    }

    fft_square(n / 2, tmp[0]);
    fft_square(n / 2, tmp[1]);
    fft_square(n / 2, tmp[2]);
    fft_square(n / 2, tmp[3]);

    for (int i = 0; i < n / 2; ++i) {
        for (int j = 0; j < n / 2; ++j) {
            auto wnp = W(n, i), wnm = W(n, j), wnpm = W(n, i + j);
            ar[i][j] = tmp[0][i][j] + tmp[1][i][j] * wnp + tmp[2][i][j] * wnm + tmp[3][i][j] * wnpm;
            ar[i + n / 2][j] = tmp[0][i][j] - tmp[1][i][j] * wnp + tmp[2][i][j] * wnm - tmp[3][i][j] * wnpm;
            ar[i][j + n / 2] = tmp[0][i][j] + tmp[1][i][j] * wnp - tmp[2][i][j] * wnm - tmp[3][i][j] * wnpm;
            ar[i + n / 2][j + n / 2] = tmp[0][i][j] - tmp[1][i][j] * wnp - tmp[2][i][j] * wnm + tmp[3][i][j] * wnpm;
        }
    }
}

void fft_cooley_tukey(int n, int m, std::vector<std::vector<std::complex<double>>> &ar) {
    if (n == m) {
        fft_square(n, ar);
        return;
    }

    std::vector<std::vector<std::complex<double>>> tmp[2];
    tmp[0] = tmp[1] = std::vector<std::vector<std::complex<double>>>(n, std::vector<std::complex<double>>(m / 2));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            tmp[j & 1][i][j / 2] = ar[i][j];
        }
    }

    fft_cooley_tukey(n, m / 2, tmp[0]);
    fft_cooley_tukey(n, m / 2, tmp[1]);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m / 2; ++j) {
            auto mult = W(m, j);
            ar[i][j] = tmp[0][i][j] + mult * tmp[1][i][j];
            ar[i][j + m / 2] = tmp[0][i][j] - mult * tmp[1][i][j];
        }
    }
}

int main() {
    freopen("input.txt", "r", stdin);
    int n, m;
    std::cin >> n >> m;
    std::vector<std::vector<std::complex<double>>> a(n, std::vector<std::complex<double>>(m));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            double t;
            std::cin >> t;
            a[i][j] = {t, 0};
        }
    }

    fft_cooley_tukey(n, m, a);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            std::cout << a[i][j].real() << " " << a[i][j].imag() << "     ";
        }
        std::cout << '\n';
    }
}
