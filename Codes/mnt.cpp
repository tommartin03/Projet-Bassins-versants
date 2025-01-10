
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <omp.h>
#include "mnt.h"
#include "fonctions.h"

using namespace std;

mnt::mnt(std::string filename) {
    std::ifstream f(filename);
    if (!f) { // vérification de l'ouverture du fichier
        std::cerr << "Error: Unable to open file " << filename << std::endl;
        exit(EXIT_FAILURE);
    }

    f >> nb_lignes >> nb_cols;
    int tmp;
    f >> tmp >> tmp >> tmp >> no_value;

    // initialisation dynamique des tableaux
    terrain = new float[(nb_lignes + 2) * nb_cols];
    direction = new int[(nb_lignes + 2) * nb_cols];
    accumulation = new int[(nb_lignes + 2) * nb_cols];
    bassin = new int[nb_lignes * nb_cols];

    
    #pragma omp parallel for collapse(2) // utitilisation de omp parallel for pour paralleliser les boucles et collapse(2) pour les boucles imbriquées 
    for (int i = 0; i < nb_lignes + 2; i++) { // initialisation du terrain, de la direction et de l'accumulation dans la même boucle
        for (int j = 0; j < nb_cols; j++) {
            terrain[i * nb_cols + j] = (i == 0 || i == nb_lignes + 1) ? no_value : 0; 
            direction[i * nb_cols + j] = (i == 0 || i == nb_lignes + 1) ? no_dir_value : -1;
            accumulation[i * nb_cols + j] = (i == 0 || i == nb_lignes + 1) ? 0 : -1;
        }
    }

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nb_lignes; i++) { // initialisation du bassin 
        for (int j = 0; j < nb_cols; j++) {
            bassin[i * nb_cols + j] = -1;
        }
    }


    for (int i = 0; i < nb_lignes; i++) { // chargement des valeurs du terrain (séquentiel)
        for (int j = 0; j < nb_cols; j++) {
            f >> terrain[(i + 1) * nb_cols + j];
        }
    }
}

void mnt::affichageTerrain() {
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            std::cout << terrain[(i + 1) * nb_cols + j] << " ";
        }
        std::cout << std::endl;
    }
}

void mnt::calculDirection() {
    #pragma omp parallel for  // initialisation des premières et dernières lignes pour la direction
    for (int j = 0; j < nb_cols; j++) {
        direction[j] = 1;
        direction[(nb_lignes + 1) * nb_cols + j] = 5; 
    }

    #pragma omp parallel for collapse(2) // calul des direction pour la grille principale
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            int x = i + 1;
            int y = j;
            float val = terrain[x * nb_cols + y];

            if (val != no_value) {
                if (j == 0) {
                    float tab_bord[5] = {
                        terrain[(x - 1) * nb_cols + y], // calcul des valeurs des bords pour chaque cellule de la grille avec des inidice simplifiés
                        terrain[(x - 1) * nb_cols + (y + 1)],
                        terrain[x * nb_cols + (y + 1)],
                        terrain[(x + 1) * nb_cols + (y + 1)],
                        terrain[(x + 1) * nb_cols + y]
                    };
                    direction[x * nb_cols + y] = f_bord1(val, tab_bord, no_value);
                } else if (j == nb_cols - 1) {
                    float tab_bord[5] = {
                        terrain[(x + 1) * nb_cols + y],
                        terrain[(x + 1) * nb_cols + (y - 1)],
                        terrain[x * nb_cols + (y - 1)],
                        terrain[(x - 1) * nb_cols + (y - 1)],
                        terrain[(x - 1) * nb_cols + y]
                    };
                    direction[x * nb_cols + y] = f_bord2(val, tab_bord, no_value);
                } else {
                    float tab[8] = {
                        terrain[(x - 1) * nb_cols + y],
                        terrain[(x - 1) * nb_cols + (y + 1)],
                        terrain[x * nb_cols + (y + 1)],
                        terrain[(x + 1) * nb_cols + (y + 1)],
                        terrain[(x + 1) * nb_cols + y],
                        terrain[(x + 1) * nb_cols + (y - 1)],
                        terrain[x * nb_cols + (y - 1)],
                        terrain[(x - 1) * nb_cols + (y - 1)]
                    };
                    direction[x * nb_cols + y] = f(val, tab, no_value);
                }
            }
        }
    }
}

void mnt::calculAccumulation() {
    int dir_bord1[5] = {5, 6, 7, 8, 1};
    int dir_bord2[5] = {1, 2, 3, 4, 5};
    int dir_general[8] = {5, 6, 7, 8, 1, 2, 3, 4};

    unsigned stop;
    int iteration_count = 0;
    const int max_iterations = 1000;

    do {
        stop = 0;
        iteration_count++;

 
        #pragma omp parallel for collapse(2) reduction(+ : stop) // parallelisation de la boucle de propagation, reduction de la variable stop pour éviter les conflits de lecture et d'écriture 
        for (int i = 0; i < nb_lignes; i++) {
            for (int j = 0; j < nb_cols; j++) {
                int x = i + 1;
                int y = j;
                int d = direction[x * nb_cols + y];

                if (d != no_dir_value && accumulation[x * nb_cols + y] == -1) { // si la direction est valide et que l'accumulation n'est pas encore calculée
                    int res = -1;

                    if (j == 0) {
                        int tab_bord[10] = {
                            accumulation[(x - 1) * nb_cols + y], //
                            accumulation[(x - 1) * nb_cols + (y + 1)],
                            accumulation[x * nb_cols + (y + 1)],
                            accumulation[(x + 1) * nb_cols + (y + 1)],
                            accumulation[(x + 1) * nb_cols + y],
                            direction[(x - 1) * nb_cols + y],
                            direction[(x - 1) * nb_cols + (y + 1)],
                            direction[x * nb_cols + (y + 1)],
                            direction[(x + 1) * nb_cols + (y + 1)],
                            direction[(x + 1) * nb_cols + y]
                        };
                        res = f_acc(tab_bord, no_dir_value, dir_bord1, 5);
                    } else if (j == nb_cols - 1) {
                        int tab_bord[10] = {
                            accumulation[(x + 1) * nb_cols + y],
                            accumulation[(x + 1) * nb_cols + (y - 1)],
                            accumulation[x * nb_cols + (y - 1)],
                            accumulation[(x - 1) * nb_cols + (y - 1)],
                            accumulation[(x - 1) * nb_cols + y],
                            direction[(x + 1) * nb_cols + y],
                            direction[(x + 1) * nb_cols + (y - 1)],
                            direction[x * nb_cols + (y - 1)],
                            direction[(x - 1) * nb_cols + (y - 1)],
                            direction[(x - 1) * nb_cols + y]
                        };
                        res = f_acc(tab_bord, no_dir_value, dir_bord2, 5);
                    } else {
                        int tab[16] = {
                            accumulation[(x - 1) * nb_cols + y],
                            accumulation[(x - 1) * nb_cols + (y + 1)],
                            accumulation[x * nb_cols + (y + 1)],
                            accumulation[(x + 1) * nb_cols + (y + 1)],
                            accumulation[(x + 1) * nb_cols + y],
                            accumulation[(x + 1) * nb_cols + (y - 1)],
                            accumulation[x * nb_cols + (y - 1)],
                            accumulation[(x - 1) * nb_cols + (y - 1)],
                            direction[(x - 1) * nb_cols + y],
                            direction[(x - 1) * nb_cols + (y + 1)],
                            direction[x * nb_cols + (y + 1)],
                            direction[(x + 1) * nb_cols + (y + 1)],
                            direction[(x + 1) * nb_cols + y],
                            direction[(x + 1) * nb_cols + (y - 1)],
                            direction[x * nb_cols + (y - 1)],
                            direction[(x - 1) * nb_cols + (y - 1)]
                        };
                        res = f_acc(tab, no_dir_value, dir_general, 8);
                    }

                    if (res != -1) {
                        accumulation[x * nb_cols + y] = res;
                    } else {
                        stop++;
                    }
                }
            }
        }

        std::cout << "Iteration: " << iteration_count << ", Stop: " << stop << std::endl;
        //affichageAccumulation();

        if (iteration_count >= max_iterations) {
            std::cerr << "Error: Maximum iterations reached in calculAccumulation." << std::endl;
            break;
        }

    } while (stop > 0);
}

void mnt::calculBassin() {
    #pragma omp parallel for collapse(2)
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            int x = i + 1;
            int y = j;

            if (direction[x * nb_cols + y] != no_dir_value) {
                if (bassin[i * nb_cols + j] == -1) { // si le bassin n'est pas encore assigné
                    bassin[i * nb_cols + j] = f_bassin(
                        bassin, direction + nb_cols, nb_lignes, nb_cols, i, j, &num);
                }
            } else {
                bassin[i * nb_cols + j] = no_bassin_value; // si la direction est absente, on assigne une valeur spéciale
            }
        }
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

void mnt::affichageBassin() {
    for (int i = 0; i < nb_lignes; i++) {
        for (int j = 0; j < nb_cols; j++) {
            std::cout << bassin[i * nb_cols + j] << " ";
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
