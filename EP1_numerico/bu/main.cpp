#include <stdio.h>
#include <iostream>
#include <vector>
#include <math.h>       /* potenciação */
#include <matriz.h>
using namespace std;

double modulo(double numero) {
    if (numero < 0) {
        return -numero;
    } else {
        return numero;
    }
}

vector<double> cos_sen_Givens(double w_ik, double w_jk) {
    double tau;
    double cos;
    double sen;
    if (modulo(w_ik) > modulo(w_jk)) {
        tau = -(w_jk)/(w_ik);
        cos = 1/pow(1+tau*tau, 0.5);
        sen = cos*tau;
    } else {
        tau = -(w_ik)/(w_jk);
        sen = 1/pow(1+tau*tau, 0.5);
        cos = sen*tau;
    }
    vector<double> cos_sen = {cos, sen};
    return cos_sen;
}

void rotacao_Givens(matriz &W, int i, int j, double cos, double sen) {
    double aux;
    for (int r = 0; r <  W.get_n_colunas(); r++) {
        aux = cos*W.get_elemento(i, r) - sen*W.get_elemento(j, r);
        W.mudar_elemento(sen*W.get_elemento(i, r) + cos*W.get_elemento(j, r), j, r);
        W.mudar_elemento(aux, i, r);
    }
}

vector<double> resolver_por_QR(matriz &W, matriz &b) {
    int i;
    int m = W.get_n_colunas();
    int n = W.get_n_linhas();
    double cos, sen, aux;
    vector<double> solucao_inversa;
    vector<double> solucao;
    for (int k = 0; k < m; k++) {
        for (int j = n - 1; j > k; j--) {
            i = j - 1;
            if (W.get_elemento(j, k)!= 0) {
                cos = cos_sen_Givens(W.get_elemento(i, k), W.get_elemento(j, k)).at(0);
                sen = cos_sen_Givens(W.get_elemento(i, k), W.get_elemento(j, k)).at(1);
                rotacao_Givens(W, i, j, cos, sen);
                rotacao_Givens(b, i, j, cos, sen);
            }
        }
    }
    for (int k = m - 1; k >= 0; k--) {
        aux = 0;
        for (int j = k; j < m - 1; j++) {
            aux += W.get_elemento(k, j)*solucao_inversa.at(m-1-j);
        }
        solucao_inversa.insert(solucao_inversa.begin(), (b.get_elemento(k, 0) - aux)/W.get_elemento(k, k));
    }
    for (int i = 0; i < m; i++) {
        solucao.push_back(solucao_inversa.at((m-1)-i));
    }
    return solucao;
}

int main() {
    vector<double> teste1 = {1, 2, 3};
    vector<double> teste2 = {1, 1, 2};
    vector<double> teste3 = {6};
    vector<double> teste4 = {4};
    vector<double> teste7 = {1, 1, 1};
    vector<double> teste8 = {3};
    vector<vector<double> > teste5 = {teste1, teste2, teste7};
    vector<vector<double> > teste6 = {teste3, teste4, teste8};
    matriz teste_matriz1(teste5);
    matriz teste_matriz2(teste6);
    cout << resolver_por_QR(teste_matriz1, teste_matriz2).at(0) << '\n';
    cout << resolver_por_QR(teste_matriz1, teste_matriz2).at(1) << '\n';
    cout << resolver_por_QR(teste_matriz1, teste_matriz2).at(2) << '\n';
}
