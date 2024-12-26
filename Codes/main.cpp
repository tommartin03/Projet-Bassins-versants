#include <iostream>

#include "mnt.h"

int main(int argc, char* argv[]) {
    std::string filename = argv[1];
    mnt M1(filename);
    std::cout << "le terrain" << std::endl;
    M1.affichageTerrain();
    M1.calculDirection();
    std::cout << "les directions" << std::endl;
    M1.affichageDirection();
    M1.calculAccumulation();
    std::cout << "l'accumulation" << std::endl;
    M1.affichageAccumulation();
    M1.calculBassin();
    std::cout << "le bassin" << std::endl;
    M1.affichageBassin();
}
