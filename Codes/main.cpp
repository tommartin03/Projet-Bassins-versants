#include <iostream>
#include <omp.h>  // OpenMP Header
#include "mnt.h"

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <filename>" << std::endl;
        return 1;
    }

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

    return 0;
}

