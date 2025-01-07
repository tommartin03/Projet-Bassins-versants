//
// Created by sophie on 20/11/24.
//
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <omp.h>  // Include OpenMP header
#include "mnt.h"
#include "fonctions.h"
using namespace std;

mnt::mnt(string filename) {
    ifstream f(filename);
    if (!f) {
        cerr << "Error opening file: " << filename << endl;
        exit(EXIT_FAILURE);
    }

    f >> nb_lignes >> nb_cols;
    int tmp;
    f >> tmp >> tmp >> tmp >> no_value;

    this->terrain = new float[(nb_lignes + 2) * nb_cols];
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            f >> terrain[(i + 1) * nb_cols + j];
        }
    }

    #pragma omp parallel for
    for (int j = 0; j < nb_cols; j++) {
        terrain[j] = no_value;
        terrain[(nb_lignes + 1) * nb_cols + j] = no_value;
    }

    this->direction = new int[(nb_lignes + 2) * nb_cols];
    this->accumulation = new int[(nb_lignes + 2) * nb_cols];
    this->bassin = new int[nb_lignes * nb_cols];

    #pragma omp parallel for
    for (int i = 0; i < nb_cols; i++) {
        direction[i] = 1;
        direction[(nb_lignes + 1) * nb_cols + i] = 5;
        accumulation[i] = 0;
        accumulation[(nb_lignes + 1) * nb_cols + i] = 0;
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            accumulation[(i + 1) * nb_cols + j] = -1;
            bassin[i * nb_cols + j] = -1;
        }
    }
}

void mnt::affichageTerrain() {
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++)
            cout << terrain[(i + 1) * nb_cols + j] << " ";
        cout << endl;
    }
}

void mnt::calculDirection() {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            int x = i + 1;
            int y = j;
            int x1 = x - 1, x2 = x + 1;
            int y1 = y - 1, y2 = y + 1;
            float val = terrain[x * nb_cols + y];

            if (val != no_value) {
                if (j == 0) {
                    float tab_bord[5] = {terrain[x1 * nb_cols + y], terrain[x1 * nb_cols + y2], terrain[x * nb_cols + y2], terrain[x2 * nb_cols + y2], terrain[x2 * nb_cols + y]};
                    direction[x * nb_cols + y] = f_bord1(val, tab_bord, no_value);
                } else if (j == nb_cols - 1) {
                    float tab_bord[5] = {terrain[x2 * nb_cols + y], terrain[x2 * nb_cols + y1], terrain[x * nb_cols + y1], terrain[x1 * nb_cols + y1], terrain[x1 * nb_cols + y]};
                    direction[x * nb_cols + y] = f_bord2(val, tab_bord, no_value);
                } else {
                    float tab[8] = {terrain[x1 * nb_cols + y], terrain[x1 * nb_cols + y2], terrain[x * nb_cols + y2], terrain[x2 * nb_cols + y2], terrain[x2 * nb_cols + y], terrain[x2 * nb_cols + y1], terrain[x * nb_cols + y1], terrain[x1 * nb_cols + y1]};
                    direction[x * nb_cols + y] = f(val, tab, no_value);
                }
            }
        }
    }
}

void mnt::calculAccumulation() {
    bool stop;
    int iteration = 0; // Track iterations for debugging
    do {
        stop = false;
        iteration++;
        std::cout << "Iteration: " << iteration << std::endl;

        #pragma omp parallel for collapse(2) reduction(||:stop)
        for (int i = 0; i < nb_lignes; i++) {
            for (int j = 0; j < nb_cols; j++) {
                int x = i + 1;
                int y = j;
                int d = direction[x * nb_cols + y];
                if (d != no_dir_value && accumulation[x * nb_cols + y] == -1) {
                    accumulation[x * nb_cols + y] = 1; // Example update for debugging
                    stop = true;  // Ensure loop continues
                }
            }
        }
    } while (stop && iteration < 1000); // Add a safety limit
}


void mnt::calculBassin() {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            if (direction[i * nb_cols + j] != no_dir_value && bassin[i * nb_cols + j] == -1) {
                bassin[i * nb_cols + j] = f_bassin(bassin, direction, nb_lignes, nb_cols, i, j, &num);
                #pragma omp critical
                {
                    std::cout << "Bassin assigned at (" << i << ", " << j << "): " << bassin[i * nb_cols + j] << std::endl;
                }
            }
        }
    }
}


void mnt::affichageBassin() {
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            cout << bassin[i * nb_cols + j] << " ";
        }
        cout << endl;
    }
}
void mnt::affichageDirection() {
    for (int i = 0; i < nb_lignes + 2; i++) {
        for (int j = 0; j < nb_cols; j++) {
            std::cout << direction[i * nb_cols + j] << " ";
        }
        std::cout << std::endl;
    }
}

void mnt::affichageAccumulation() {
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            std::cout << accumulation[(i + 1) * nb_cols + j] << " ";
        }
        std::cout << std::endl;
    }
}


mnt::~mnt() {
    delete[] terrain;
    delete[] direction;
    delete[] accumulation;
    delete[] bassin;
}

