/*The algorithm is based on the
original paper by Deb et al. (2002) titled 
"A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II
*/

//The stucture for the code is as follows
#include <vector>
#include <algorithm>
#include "individual.hpp"
/*
1. Non-Dominated Sorting
We select population above a certain objective threshold. 
*/
#include "ndsort.cpp"
/* 
2. Crowding Distance
Within the population selected with non-dominated sorting, 
sort each front with crowding distance. 
This forms the new population
*/
#include "crowding.cpp"
/*
3. NSGA-II
The main function that combines the above two functions
Export this function to Python
*/
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

using namespace std;
vector<individual> nsga2(vector<individual>& population, int population_size) {
    /*pseudocode equivalent
    fronts = non_dominated_argsort(population)
    new_population = [] and i = 0
    while |new_population| + |fronts[i]| <= population_size:
        crowding_distance(population, fronts[i])
        new_population += fronts[i]
        i ++
    sort(fronts[i] by crowding_distance)
    new_population += fronts[i][:population_size - |new_population|]
    */
    vector<individual> new_population; // new_population = []
    new_population.reserve(population_size);
    vector<vector<int>> fronts = non_dominated_argsort(population); 
    // fronts = non_dominated_argsort(population)
    size_t i = 0; // i = 0
    while (static_cast<int>(new_population.size() + fronts[i].size())
            <= population_size) { // while |new_population| + |fronts[i]| <= population_size
        crowding_distance(population, fronts[i]); // crowding_distance(population, fronts[i])
        // new_population += fronts[i]
        for (int idx: fronts[i]) {
            new_population.push_back(population[idx]);
        }
        ++i ;
        if (i >= fronts.size()) {
            break;
        }
    }

    sort(fronts[i].begin(), fronts[i].end(), [&](int a, int b) {
        return crowd_compare(population[a], population[b]);
    }); // sort(fronts[i] by crowding_distance)

    for(size_t j = 0; j < population_size - new_population.size(); ++j) {
        new_population.push_back(population[fronts[i][j]]);
    }
    // new_population += fronts[i][:population_size - |new_population|]
    return new_population;
}

namespace py = pybind11;
using namespace py;

PYBIND11_MODULE(nsga2, m) {
    m.doc() = "NSGA-II algorithm";
    class_<individual>(m, "Individual")
        .def(py::init<vector<double>>())
        .def_readwrite("objectives", &individual::objectives)
        .def_readonly("domination_count", &individual::domination_count)
        .def_readonly("dominated_solutions", &individual::dominated_solutions)
        .def_readonly("distance", &individual::distance)
        .def_readonly("rank", &individual::rank);
    m.def("nsga2", &nsga2, "Run the NSGA-II algorithm");
}