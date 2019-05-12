#include <stdio.h>
#include <time.h>       /* time */
#include <stdlib.h>     /* srand, rand */
#include <iostream>
//#include <vector>
#include <math.h>       /* potenciação */
#include <Matriz.h>
#include <vector>
#include <string.h>
using namespace std;

double modulo(double numero) { // função módulo padrão
    if (numero < 0) {
        return -numero;
    } else {
        return numero;
    }
}

vector<double> cos_sen_Givens(double w_ik, double w_jk) { // calcula o cosseno e o seno de Givens. Retorna uma vector onde o primeiro elemenento é o cos e o segundo é o sen
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

void rotacao_Givens(Matriz &matriz_W, int i, int j, double cos, double sen) {  // efetua a rotacao de Givens em duas linhas, dados cos e sen de Givens
    double aux;
    for (int r = 0; r <  matriz_W.get_n_colunas(); r++) {
        aux = cos*matriz_W.get_elemento(i, r) - sen*matriz_W.get_elemento(j, r);
        matriz_W.mudar_elemento(sen*matriz_W.get_elemento(i, r) + cos*matriz_W.get_elemento(j, r), j, r);
        matriz_W.mudar_elemento(aux, i, r);
    }
}

vector<vector<double> > resolver_por_QR(Matriz &matriz_W, Matriz &matriz_A) { // resolve varios sistemas lineares por QR
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

Matriz resolver_por_QR_m(Matriz &matriz_W, Matriz &matriz_A) { // resolve varios sistemas lineares por QR
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
    return solucao;
}



vector <double> norma_euclidiana_por_coluna_raiz(Matriz matriz){
    vector <double> vetor_norma(matriz.get_n_colunas());
    for (int i = 0; i < matriz.get_n_linhas(); i++){
        for (int j = 0; j < matriz.get_n_colunas(); j++){
            vetor_norma.at(j) += matriz.get_elemento(i,j)*matriz.get_elemento(i,j);
        }
    }
    for (int j = 0; j < vetor_norma.size(); j++) {
            vetor_norma.at(j) = pow(vetor_norma.at(j), 0.5);
    }
    return vetor_norma;
}


Matriz achar_W_e_H(Matriz &matriz_W, Matriz &matriz_A) {
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
    double erro_atual;
    double erro_anterior;
    double aux;
    int t1 = time(0);
    for (int i = 0; i < n_linhas_W; i++) {
        for (int j = 0; j < n_colunas_W; j++) {
            matriz_W.mudar_elemento(rand(), i, j);
        }
    }
    vector <double> norma;

    while (1) {
        Matriz matriz_A_modificada = matriz_A;
        norma = norma_euclidiana_por_coluna_raiz(matriz_W);
        for (int i = 0; i < n_linhas_W; i++) {
            for (int j = 0; j < n_colunas_W; j++) {
                matriz_W.mudar_elemento(matriz_W.get_elemento(i, j)/norma.at(j), i, j);
            }
        }
        Matriz matriz_H = resolver_por_QR_m(matriz_W, matriz_A_modificada);  // 6 a 7 segundos
        int n_linhas_H = matriz_H.get_n_linhas();
        int n_colunas_H = matriz_H.get_n_colunas();
        for (int i = 0; i < n_linhas_H; i ++) {
            for (int j = 0; j < n_colunas_H; j ++) {
                if (matriz_H.get_elemento(i, j) < 0) {
                    matriz_H.mudar_elemento(0, i, j);
                }
            }
        }
        Matriz matriz_H_transposta = matriz_H.transposta();
        Matriz matriz_A_transposta = matriz_A.transposta();
        matriz_W.mudar_a(resolver_por_QR_m(matriz_H_transposta, matriz_A_transposta).transposta().get_a()); // 6 a 7 segundos
        for (int i = 0; i < n_linhas_W; i ++) {
            for (int j = 0; j < n_colunas_W; j ++) {
                if (matriz_W.get_elemento(i, j) < 0) {
                    matriz_W.mudar_elemento(0, i, j);
                }
            }
        }
        //cout << "Tempo uma iteracao dentro do while: " << time(0) - t1 << endl;

/*
        somatoria_quadrados = 0;
        Matriz matriz_WH = matriz_W*(matriz_H);
        for (int i = 0; i < n_linhas_A; i++) {
            for (int j = 0; j < n_colunas_A; j++) {
                aux = matriz_A.get_elemento(i, j) - matriz_WH.get_elemento(i, j);
                somatoria_quadrados += aux*aux;
            }
        }

        if (n_iteracoes == 0) {
            erro_atual = somatoria_quadrados;
            erro = erro_atual;
        } else {
            erro_anterior = erro_atual;
            erro_atual = somatoria_quadrados;
            erro = modulo(erro_atual - erro_anterior);
        }
*/

        if (n_iteracoes%10 == 0){
            cout << "Iteracao " << n_iteracoes << endl;
            //cout << "Erro = " << erro << endl << endl;
        }

        n_iteracoes++;
        if (n_iteracoes >= itmax || erro <= epsilon) {
            cout << "Tempo uma iteracao while: " << time(0) - t1 << endl;
            return matriz_H;
        }
    }
}

Matriz armazenamento(char arq[]){
    FILE *fp;
    fp = fopen(arq,"r");
    int linha = 784;
    vector <int> dado (8000000);
    long int tamanho = 0;
    while(!feof(fp)){
        fscanf(fp,"%d",&dado[tamanho]);
        tamanho++;
    }
    fclose(fp);
    long int coluna = tamanho/linha;
    long int k = 0;
    Matriz *matriz_aux = new Matriz(linha,coluna);
    for (int i=0;i<linha;i++){
        for(int j=0; j<coluna; j++){
            matriz_aux->mudar_elemento((dado.at(k))/255.0,i,j);
            ++k;
        }
    }
    return *matriz_aux;
}

vector <double> armazenamento_index(char arq[]){
    FILE *fp;
    fp = fopen(arq,"r");
    vector <double> dado (10000);
    long int tamanho = 0;
    while(!feof(fp)){
        fscanf(fp,"%lf",&dado[tamanho]);
        tamanho++;
    }
    fclose(fp);
    return dado;
}

Matriz matriz_de_treinamento(char arq[], int ndig_treino)//Leitura do arquivo criando a matriz de treinamento
{
    Matriz matriz_aux = armazenamento(arq);
    Matriz matriz_a = matriz_aux.pegarcoluna(ndig_treino);
    return matriz_a;
}

vector<Matriz> classificador(int ndig_treino,int p){ //Lista de classificadores 0 a 9, na ordem
    vector<Matriz> treinamento;
    int n = 784;
    for (int i = 0; i < 10; i++){
        char buff[50];
        string s = to_string(i);
        char const *pchar = s.c_str();
        strcpy(buff,"train_dig");
        strcat(buff,pchar);
        strcat(buff,".txt");
        Matriz matriz_A = matriz_de_treinamento(buff,ndig_treino); //criação da matriz A
        Matriz matriz_W = Matriz(n, p);
        Matriz matriz_H = Matriz(achar_W_e_H(matriz_W,matriz_A));
        treinamento.push_back(matriz_W);
    }
    return treinamento;
}

vector <double> norma_euclidiana_por_coluna(Matriz matriz){
    vector <double> vetor_norma(matriz.get_n_colunas());
    for (int i = 0; i < matriz.get_n_linhas(); i++){
        for (int j = 0; j < matriz.get_n_colunas(); j++){
            vetor_norma.at(j) += matriz.get_elemento(i,j)*matriz.get_elemento(i,j);
        }
    }
    return vetor_norma;
}

vector <vector <double> > teste_imagens (vector <Matriz> w){ //teste dos classificadores
    vector <vector <double> > vetores;
    vector <double> indice;
    Matriz teste = armazenamento("test_images.txt");
    for (unsigned int i = 0; i < w.size(); i++){
        Matriz matriz = w.at(i);
        cout << i << endl;
        Matriz teste_copia = teste;
        Matriz matriz_H(resolver_por_QR(matriz,teste_copia));
        int t1 = time(0);
        Matriz matriz_error = teste - w.at(i)*matriz_H;
        cout << "Tempo multiplicacao: " << time(0) - t1 << endl;
        vector <double> norma = norma_euclidiana_por_coluna(matriz_error);
        //for (int j = 0; j < norma.size(); j++){
        //    cout << norma.at(j);
       // }
        //matriz_H.show();
        //cout << " " << endl;
        if (i == 0){
            vetores.push_back(norma);
            for (unsigned int k = 0; k < norma.size(); k++){
                indice.push_back(0);
            }
            vetores.push_back(indice);
        }
        else{
            for (unsigned int k = 0; k < norma.size(); k++){
                if (vetores.at(0).at(k) > norma.at(k)){
                    vetores.at(0).at(k) = norma.at(k);
                    vetores.at(1).at(k) = i;
                }
            }
        }
    }
    return vetores;
}

int main() {
// tarefa principal
//1)Treinamento

vector <Matriz> renshu;
renshu = classificador(1000, 10);
//2)Teste
vector <vector <double> > tesuto;
tesuto = teste_imagens(renshu);

//3)Comparar com o index
vector <double> sakuin;
sakuin = armazenamento_index("test_index.txt");
//3.a)Calcular o percentual total de acertos
double counter = 0;
for (unsigned int i = 0; i < sakuin.size(); i++){
    if (tesuto.at(1).at(i) == sakuin.at(i)) counter++;
}
double percentual = 100*counter/sakuin.size();

cout << "O percentual de acerto eh " << percentual << "%" << endl;

//3.b)percentual de acerto de cada digito
for (int k = 0; k < 10; k++){
    double cont = 0;
    double total = 0;
    for (unsigned int i = 0; i < sakuin.size(); i++){
        if(sakuin.at(i) == k){
            ++total;
            if (tesuto.at(1).at(i) == sakuin.at(i)) {
                cont++;
            }
        }
    }
    cout << "A ocorrencia correta de " << k << " foi igual a " << cont << " no total de " << total << " ocorrencias"<< endl;
    cout << "O percentual relativo: " << 100*cont/total << "%    Percentual do total: "  <<  100*cont/sakuin.size() << "%"<< endl;
}


//for (unsigned int i = 0; i < sakuin.size(); i++){
//    cout << tesuto.at(1).at(i) << sakuin.at(i) << endl;
//}

// segunda tarefa
/*
    vector<double> teste1 = {0.3, 0.6, 0};
    vector<double> teste2 = {0.5, 0, 1};
    vector<double> teste3 = {0.4, 0.8, 0};
    vector<double> teste4 = {0.6, 0};
    vector<double> teste5 = {0, 1};
    vector<double> teste6 = {0.8, 0};
    vector<double> teste7 = {0.5, 1, 0};
    vector<double> teste8 = {0.5, 0, 1};
    vector<vector<double> > testeA = {teste1, teste2, teste3};
    vector<vector<double> > testeW = {teste4, teste5, teste6};
    vector<vector<double> > testeH = {teste7, teste8};
    Matriz matriz_A(testeA);
    Matriz matriz_W(3, 2);
    cout << "Segunda tarefa:" << endl;
    Matriz matriz_H = achar_W_e_H(matriz_W, matriz_A);
    cout << "Resultados com 100 iteracoes:" << endl;
    cout << "matriz_H calculada:" << endl;
    matriz_H.show();
    cout << endl;
    cout << "matriz_W calculada:" << endl;
    matriz_W.show();
    cout << endl;
    cout << "W*H:" << endl;
    (matriz_W*matriz_H).show();
    cout << endl;
*/

    // testes para o resolver_por_QR
/*
    vector<double> teste1 = {2, 1, 1, -1, 1};
    vector<double> teste2 = {0, 3, 0, 1, -2};
    vector<double> teste3 = {0, 0, 2, 2, -1};
    vector<double> teste4 = {0, 0, -1, 1, 2};
    vector<double> teste5 = {0, 0, 0, 3, 1};
    vector<vector<double> > teste6 = {teste1, teste2, teste3, teste4, teste5};
    vector<double> b1 = {8};
    vector<double> b2 = {4};
    vector<double> b3 = {6};
    vector<double> b4 = {4};
    vector<double> b5 = {8};
    vector<vector<double> > bzao = {b1, b2, b3, b4, b5};
    Matriz b(bzao);
    Matriz testezao(teste6);
    Matriz(resolver_por_QR(testezao, b)).show();
*/

    // abaixo o teste a)
/*
    cout << "Teste a:" << endl;
    cout << "Matriz W do Teste a: " << endl;

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
    teste_matriz.show(1, 0);

    cout << "Resultados do teste a: " << endl;
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show(6, 5);
*/

    // teste b
/*
    cout << "Teste b" << endl;
    cout << "Matriz W do Teste b:" << endl << endl;
    Matriz teste_matriz(20, 17);
    Matriz teste_b(20, 1);
    for (int i = 0; i < 20; i++) {
        teste_b.mudar_elemento(i+1, i, 0);
        for (int j = 0; j < 17; j++) {
            if (modulo(i-j) <= 4) {
                teste_matriz.mudar_elemento(1.0/((i+1)+(j+1)-1), i, j);
            }
        }
    }
    teste_matriz.show(6, 5);
    cout << endl << "Resultado teste b:" << endl;
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show(9, 5);
*/

    // teste c
/*
    int cte = 64;
    Matriz teste_matriz(cte, cte);
    Matriz teste_b     (cte, 3  );
    for (int i = 0; i < cte; i++) {
        teste_b.mudar_elemento( 1,         i, 0);
        teste_b.mudar_elemento( i+1,       i, 1);
        teste_b.mudar_elemento(2*(i+1)-1, i, 2);
        for (int j = 0; j < cte; j++) {
            if (i == j) {
                teste_matriz.mudar_elemento(2, i, j);
            }
            if (modulo(i-j) == 1) {
                teste_matriz.mudar_elemento(1, i, j);
            }
        }
    }
    cout << "Teste c" << endl;
    cout << "Matriz W do Teste C:" << endl;
    teste_matriz.show(1, 0);
    cout << endl << "Resultados do item c: " << endl;
    Matriz(resolver_por_QR(teste_matriz, teste_b)).show(9, 5);
*/

    // teste d
/*
    Matriz teste_matriz(20, 17);
    Matriz teste_b(20, 3);
    for (int i = 0; i < 20; i++) {
        teste_b.mudar_elemento(1, i, 0);
        teste_b.mudar_elemento(i+1, i, 1);
        teste_b.mudar_elemento(2*(i+1)-1, i, 2);
        for (int j = 0; j < 17; j++) {
            if (modulo(i-j) <= 4) {
                teste_matriz.mudar_elemento(1.0/((i+1)+(j+1)-1), i, j);
            }
        }
    }
    cout << "Teste d" << endl;
    cout << "Matriz W do Teste D:" << endl;
    teste_matriz.show(7, 5);
    cout << endl << "Resultados do item d: " << endl;
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
