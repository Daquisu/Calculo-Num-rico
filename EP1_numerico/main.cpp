#include <stdio.h>
#include <time.h>       /* time */
#include <stdlib.h>     /* srand, rand */
#include <iostream>
#include <deque>
#include <math.h>       /* potenciação */
#include <Matriz.h>
using namespace std;

double modulo(double numero) { // função módulo padrão
    if (numero < 0) {
        return -numero;
    } else {
        return numero;
    }
}

deque<double> cos_sen_Givens(double w_ik, double w_jk) { // calcula o cosseno e o seno de Givens. Retorna uma deque onde o primeiro elemenento é o cos e o segundo é o sen
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
    deque<double> cos_sen = {cos, sen};
    return cos_sen;
}

void rotacao_Givens(Matriz &matriz_W, int i, int j, double cos, double sen) {  // efetua a rotacao de Givens em duas linhas, dados cos e sen de Givens
    double aux;
    for (int r = 0; r <  matriz_W.get_n_colunas(); r++) {
        aux = cos*matriz_W.get_elemento(i, r) - sen*matriz_W.get_elemento(j, r);
        matriz_W.mudar_elemento(sen*matriz_W.get_elemento(i, r) + cos*matriz_W.get_elemento(j, r), j, r);
        matriz_W.mudar_elemento(aux, i, r);
    }
}

deque<deque<double> > resolver_por_QR(Matriz &matriz_W, Matriz &matriz_A) { // resolve varios sistemas lineares por QR
    int i;
    int p = matriz_A.get_n_colunas();
    int n = matriz_W.get_n_linhas();
    int m = matriz_W.get_n_colunas();
    double cos, sen, aux;
    Matriz solucao(m, p);
    for (int k = 0; k < m; k++) {
        for (int j = n - 1; j > k; j--) {
            i = j - 1;
            if (matriz_W.get_elemento(j, k)!= 0) {
                cos = cos_sen_Givens(matriz_W.get_elemento(i, k), matriz_W.get_elemento(j, k)).at(0);
                sen = cos_sen_Givens(matriz_W.get_elemento(i, k), matriz_W.get_elemento(j, k)).at(1);
                rotacao_Givens(matriz_W, i, j, cos, sen);
                rotacao_Givens(matriz_A, i, j, cos, sen);
            }
        }
    }

    for (int k = m - 1; k >= 0; k--) {
        for (int j = 0; j < p; j++) {
            aux = 0;
            for (int i = k; i <= m - 1; i++) {
                aux += matriz_W.get_elemento(k, i)*solucao.get_elemento(i, j);
            }
            solucao.mudar_elemento(((matriz_A.get_elemento(k, j) - aux)/matriz_W.get_elemento(k, k)), k, j);
        }
    }
    return solucao.get_a();
}

Matriz achar_W_e_H(Matriz &matriz_W, Matriz & matriz_A) {
    const int itmax = 100;
    const float epsilon = 0.00001;
    double erro = epsilon + 1;
    int n_iteracoes = 0;
    double somatoria_quadrados = 0;
    srand(time(0));
    int n_linhas_W = matriz_W.get_n_linhas();
    int n_colunas_W = matriz_W.get_n_colunas();
    int n_colunas_A = matriz_A.get_n_colunas();
    int n_linhas_A = matriz_A.get_n_linhas();
    Matriz* matriz_H;
    double erro_atual;
    double erro_anterior;

    for (int i = 0; i < n_linhas_W; i++) {
        for (int j = 0; j < n_colunas_W; j++) {
            matriz_W.mudar_elemento(rand(), i, j);
        }
    }
    while (n_iteracoes < itmax && erro > epsilon) {;
        Matriz matriz_A_modificada = matriz_A;
        for (int j = 0; j < n_colunas_W; j++) {
            somatoria_quadrados = 0;
            for (int i = 0; i < n_linhas_W; i++) {
                somatoria_quadrados += matriz_W.get_elemento(i, j)*matriz_W.get_elemento(i, j);
            }
            somatoria_quadrados = pow(somatoria_quadrados, 0.5);
            for (int i = 0; i < n_linhas_W; i++) {
                matriz_W.mudar_elemento(matriz_W.get_elemento(i, j)/somatoria_quadrados, i, j);
            }
        }
        matriz_H = new Matriz(resolver_por_QR(matriz_W, matriz_A_modificada));
        int n_linhas_H = matriz_H->get_n_linhas();
        int n_colunas_H = matriz_H->get_n_colunas();
        for (int j = 0; j < n_colunas_H; j ++) {
            for (int i = 0; i < n_linhas_H; i ++) {
                if (matriz_H->get_elemento(i, j) < 0) {
                    matriz_H->mudar_elemento(0, i, j);
                }
            }
        }

        Matriz matriz_A_transposta = matriz_A.transposta();
        Matriz matriz_H_transposta = matriz_H->transposta();
        Matriz* matriz_W_transposta = new Matriz(resolver_por_QR(matriz_H_transposta, matriz_A_transposta));
        Matriz* matriz_W_nova = new Matriz((matriz_W_transposta->transposta()).get_a());
        delete matriz_W_transposta;
        for (int j = 0; j < n_colunas_W; j ++) {
            for (int i = 0; i < n_linhas_W; i ++) {
                if (matriz_W_nova->get_elemento(i, j) < 0) {
                    matriz_W_nova->mudar_elemento(0, i, j);
                }
            }
        }
        for (int j = 0; j < n_colunas_W; j ++) {
            for (int i = 0; i < n_linhas_W; i ++) {
                matriz_W.mudar_elemento(matriz_W_nova->get_elemento(i, j), i, j);
            }
        }
        somatoria_quadrados = 0;
        for (int i = 0; i < n_linhas_A; i++) {
            for (int j = 0; j < n_colunas_A; j++) {
                somatoria_quadrados += (matriz_A.get_elemento(i, j) - (matriz_W*(*matriz_H)).get_elemento(i, j))*(matriz_A.get_elemento(i, j) - (matriz_W*(*matriz_H)).get_elemento(i, j));
            }
        }
/*
        if (n_iteracoes == 0) {
            erro_atual = somatoria_quadrados;
        } else {
            erro_anterior = erro_atual;
            erro_atual = somatoria_quadrados;
            erro = modulo(erro_atual - erro_anterior);
        }
*/


        n_iteracoes++;
    }
    return *matriz_H;
}

int main() {
// segunda tarefa

    deque<double> teste1 = {0.3, 0.6, 0};
    deque<double> teste2 = {0.5, 0, 1};
    deque<double> teste3 = {0.4, 0.8, 0};
    //deque<double> teste4 = {0.6, 0);
    //deque<double> teste5 = {0, 1);
    //deque<double> teste6 = {0.8, 0};
    //deque<double> teste7 = {0.5, 1, 0};
    //deque<double> teste8 = {0.5, 0, 1};
    deque<deque<double> > testeA = {teste1, teste2, teste3};
    //deque<deque<double> > testeW = {teste4, teste5, teste6};
    //deque<deque<double> > testeH = {teste7, teste8};
    Matriz matriz_A(testeA);

    Matriz matriz_W(3, 2);

    // matriz_W.show();
    //cout << endl;

    //matriz matriz_W(testeW);
    //matriz matriz_H(testeH);

    Matriz matriz_H = achar_W_e_H(matriz_W, matriz_A);
    cout << "Resultados utilizando 100 iteracoes" << endl;
    cout << "matriz_H calculada:" << endl;
    matriz_H.show();
    cout << endl;
    cout << "matriz_W calculada:" << endl;
    matriz_W.show();
    cout << endl;
    cout << "W*H:" << endl;
    (matriz_W*matriz_H).show();
    //cout << endl;

    //matriz_W.show();



    // testes para o resolver_por_QR
/*
    deque<double> teste1 = {2, 1, 1, -1, 1};
    deque<double> teste2 = {0, 3, 0, 1, -2};
    deque<double> teste3 = {0, 0, 2, 2, -1};
    deque<double> teste4 = {0, 0, -1, 1, 2};
    deque<double> teste5 = {0, 0, 0, 3, 1};
    deque<deque<double> > teste6 = {teste1, teste2, teste3, teste4, teste5};
    deque<double> b1 = {8};
    deque<double> b2 = {4};
    deque<double> b3 = {6};
    deque<double> b4 = {4};
    deque<double> b5 = {8};
    deque<deque<double> > bzao = {b1, b2, b3, b4, b5};
    Matriz b(bzao);
    Matriz testezao(teste6);
    Matriz(resolver_por_QR(testezao, b)).show();
*/

    // abaixo o teste a)
/*
    int cte = 64;
    Matriz teste_matriz(cte, cte);
    Matriz teste_b     (cte, 1  );
    for (int i = 0; i < cte; i++) {
        teste_b.mudar_elemento(1, i, 0);
        for (int j = 0; j < cte; j++) {
            if (i == j) {
                teste_matriz.mudar_elemento(2, i, j);
            }
            if (modulo(i-j) == 1) {
                teste_matriz.mudar_elemento(1, i, j);
            }
        }
    }
    teste_matriz.show();
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show();
*/

    // teste b
    /*
    Matriz teste_matriz(20, 17);
    Matriz teste_b(20, 1);
    for (int i = 0; i < 20; i++) {
        teste_b.mudar_elemento(i, i, 0);
        for (int j = 0; j < 17; j++) {
            if (modulo(i-j) <= 4) {
                teste_matriz.mudar_elemento(1.0/((i+1)+(j+1)-1), i, j);
            }
        }
    }
    teste_matriz.show();
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show();
    */

    // teste c
    /*
    int cte = 64;
    Matriz teste_matriz(cte, cte);
    Matriz teste_b     (cte, 3  );
    for (int i = 0; i < cte; i++) {
        teste_b.mudar_elemento(1, i, 0);
        teste_b.mudar_elemento(i, i, 1);
        teste_b.mudar_elemento((2*i)-1, i, 2);
        for (int j = 0; j < cte; j++) {
            if (i == j) {
                teste_matriz.mudar_elemento(2, i, j);
            }
            if (modulo(i-j) == 1) {
                teste_matriz.mudar_elemento(1, i, j);
            }
        }
    }
    teste_matriz.show();
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show();
    */

    // teste d
    /*
    Matriz teste_matriz(20, 17);
    Matriz teste_b(20, 3);
    for (int i = 0; i < 20; i++) {
        teste_b.mudar_elemento(1, i, 0);
        teste_b.mudar_elemento(i, i, 1);
        teste_b.mudar_elemento((2*i)-1, i, 2);
        for (int j = 0; j < 17; j++) {
            if (modulo(i-j) <= 4) {
                teste_matriz.mudar_elemento(1.0/((i+1)+(j+1)-1), i, j);
            }
        }
    }
    teste_matriz.show();
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show();
    */
    // teste operador =
    /*
    Matriz a(2, 2);
    a.mudar_elemento(4, 0, 0);
    Matriz b = a;
    a.mudar_elemento(2, 0, 0);
    a.show();
    cout << endl;
    b.show();
    */
}
