#include "Matriz.h"
#include <iomanip>      // std::setw

Matriz::Matriz() { // declara

}

Matriz::Matriz( vector<vector<double> > lista_de_listas) { // cria uma matriz por meio de uma lista de listas
    a = lista_de_listas;
    n_linhas = a.size();
    n_colunas = a.at(0).size();
}

Matriz::Matriz(int n_linhas, int n_colunas) { // cria uma matriz com elementos nulos, definindo apenas o número de linhas e de colunas
    vector<double> aux;
    this->n_linhas = n_linhas;
    this->n_colunas = n_colunas;
    for (int j = 0; j < n_colunas; j++) {
        aux.push_back(0);
    }
    for (int i = 0; i < n_linhas; i++) {
        a.push_back(aux);
    }
}

Matriz::~Matriz() {
    //dtor
}

int Matriz::get_n_colunas() {
    return n_colunas;
}

int Matriz::get_n_linhas() {
    return n_linhas;
}

vector<vector<double> > Matriz::get_a() {
    return a;
}

double Matriz::get_elemento(int i, int j) {
    return a.at(i).at(j);
}


void Matriz::show(int digitos_total, int digitos_depois_virgula) { // os dois parametros sao para a formatacao do printf. as matrizes ficam alinhadas
    for (int i = 0; i < n_linhas; i++) {
        for (int j = 0; j < n_colunas; j++) {
            // cout << a.at(i).at(j) << "     ";
            printf ("%*.*f ", digitos_total, digitos_depois_virgula, a.at(i).at(j)); // primeiro digito eh o total de digitos. o segundo eh digitos apos a virgula
        }
        printf ("\n");
    }
}

void Matriz::show() { // parametros padroes
    for (int i = 0; i < n_linhas; i++) {
        for (int j = 0; j < n_colunas; j++) {
            // cout << a.at(i).at(j) << "     ";
            printf ("%10.5f ", a.at(i).at(j));
        }
        printf ("\n");
    }
}

void Matriz::mudar_elemento(double a_ij, int i, int j) { // set_elemento basicamente
    a.at(i).at(j) = a_ij;
}

Matriz Matriz::transposta() {
    int linhas = n_colunas;
    int colunas = n_linhas;
    Matriz matriz_retorno(linhas, colunas);
    for (int i = 0; i < linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            matriz_retorno.mudar_elemento(a.at(j).at(i), i, j);
        }
    }
    return matriz_retorno;
}

void Matriz::mudar_a(vector<vector<double> > a) {
    this->a = a;
}

Matriz Matriz::pegarcoluna(int colunas){
    Matriz matriz_parte(n_linhas,colunas);
    for (int i = 0; i < n_linhas; i++) {
        for (int j = 0; j < colunas; j++) {
            matriz_parte.mudar_elemento(a.at(i).at(j), i, j);
        }
    }
    return matriz_parte;
}

Matriz Matriz::operator *(Matriz &outra) { // multiplicação de matriz
    vector<vector<double> > futura_matriz;
    if (n_colunas == outra.get_n_linhas()) {
        for (int i = 0; i < n_linhas; i++) {
            vector<double> vetor_linha;
            futura_matriz.push_back(vetor_linha);
            for (int j = 0; j < outra.get_n_colunas(); j++) {
                double temp = 0;
                for (int k = 0; k < n_colunas; k++) {
                    temp += a.at(i).at(k) * outra.a.at(k).at(j);
                }
                futura_matriz.at(i).push_back(temp);
            }
        }
    }
    return Matriz(futura_matriz);
}

Matriz Matriz::operator =(Matriz const&outra) { // atribuição de matriz
    Matriz m_return = Matriz(outra.n_linhas, outra.n_colunas);
    m_return.a = outra.a;
    return m_return;
}

Matriz Matriz::operator -(Matriz const&outra) { // subtração de matrizes
    vector<vector<double> > futura_matriz;
    for (int i = 0; i < outra.n_linhas; i++) {
        vector<double> vetor_linha;
        futura_matriz.push_back(vetor_linha);
        for (int j = 0; j < outra.n_colunas; j++) {
            futura_matriz.at(i).push_back(a.at(i).at(j)-outra.a.at(i).at(j));
        }
    }
    return Matriz(futura_matriz);
}



