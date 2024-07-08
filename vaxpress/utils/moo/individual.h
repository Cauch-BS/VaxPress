/* The individual is based on the original paper by Deb et al. (2002) titled 
"A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"*/
#ifndef INDIVIDUAL_H
#define INDIVIDUAL_H
#include <vector>
using namespace std;

struct individual {
    vector<double> objectives; // objective values
    int domination_count = 0; // number of individuals that dominate this individual
    vector<int> dominated_solutions; // indices of individuals that this individual dominates
    double distance = 0.0; // crowding distance
    int rank = 0; // rank of the individual (front number)

    individual(vector<double> obj) : objectives(move(obj)) {}
};

#endif // INDIVIDUAL_H