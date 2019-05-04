#include "Matriz.h"

Matriz::Matriz() { // declara

}

Matriz::Matriz( deque<deque<double> > lista_de_listas) { // cria uma matriz por meio de uma lista de listas
    a = lista_de_listas;
    n_linhas = a.size();
    n_colunas = a.at(0).size();
}

Matriz::Matriz(int n_linhas, int n_colunas) { // cria uma matriz com elementos nulos, definindo apenas o número de linhas e de colunas
    deque<double> aux;
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

deque<deque<double> > Matriz::get_a() {
    return a;
}

double Matriz::get_elemento(int i, int j) {
    return a.at(i).at(j);
}


void Matriz::show() { // printa a matriz, pode ser melhorado pra ficar mais bonito mas não é prioridade
    for (int i = 0; i < n_linhas; i++) {
        for (int j = 0; j < n_colunas; j++) {
            cout << a.at(i).at(j) << "     ";
        }
        cout << '\n';
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

Matriz Matriz::operator *(Matriz &outra) { // multiplicação de matriz
    deque<deque<double> > futura_matriz;
    if (n_colunas == outra.get_n_linhas()) {
        for (int i = 0; i < n_linhas; i++) {
            deque<double> vetor_linha;
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
    deque<deque<double> > futura_matriz;
    for (int i = 0; i < outra.n_linhas; i++) {
        deque<double> vetor_linha;
        futura_matriz.push_back(vetor_linha);
        for (int j = 0; j < outra.n_colunas; j++) {
            futura_matriz.at(i).push_back(outra.a.at(i).at(j));
        }
    }
    return Matriz(futura_matriz);
}

