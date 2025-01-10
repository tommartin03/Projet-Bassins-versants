//
// Created by sophie on 20/11/24.
//

#ifndef MNT_H
#define MNT_H
#include <stdlib.h>

class mnt {
private:
    float no_value; // absence de valeur dans le MNT initial (cf le fichier de données)
    int nb_lignes;
    int nb_cols;
    int no_dir_value = -9; // choisi pour initialiser l'absence de données dans la direction
    int no_bassin_value = -10;
    float *terrain;
    int *direction;
    int *accumulation;
    int *bassin;
    int num = 100; // début de l'étiquetage poru les bassins versants

public:
    mnt(std::string filename);
    void affichageTerrain();
    void calculDirection();
    void affichageDirection();
    void calculAccumulation();
    void affichageAccumulation();
    void calculBassin();
    void affichageBassin();
    ~mnt(void);
};




#endif //MNT_H
