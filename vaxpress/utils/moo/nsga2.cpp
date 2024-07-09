/**
 * @file nsga2.cpp
 * @brief Implementation of the NSGA-II algorithm.
 * 
 * The algorithm is based on the original paper by Deb et al. (2002) titled 
 * "A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II".
 * 
 * The structure of the code is as follows:
 * 1. Non-Dominated Sorting: Select population above a certain objective threshold.
 * 2. Crowding Distance: Sort each front with crowding distance within the selected population.
 * 3. NSGA-II: The main function that combines the above two functions.
 * 
 * This file also includes the necessary headers and defines the NSGA-II algorithm as a Python module.
 */

#include <vector>
#include <algorithm>
#include "individual.hpp"
#include "ndsort.cpp"
#include "crowding.cpp"
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/complex.h>

using namespace std;

/**
 * @brief Runs the NSGA-II algorithm.
 * 
 * @param population The initial population of individuals.
 * @param population_size The desired size of the new population.
 * @return The new population of individuals.
 */

//For debugging purposes
// #include <iostream>
// #include <iterator>

//For debugging purposes
// template<typename T>
// std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
//     os << "[";
//     std::copy(vec.begin(), vec.end(), std::ostream_iterator<T>(os, ", "));
//     if (!vec.empty()) {
//         os << "\b\b"; // Remove the last ", "
//     }
//     os << "]";
//     return os;
// }

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

    // //For debugging purposes
    // std::cout << "Fronts size: " << fronts.size() << endl;

    // //For debugging purposes
    // for (size_t i = 0; i < fronts.size(); ++i) {
    //    std::cout << "Front " << i + 1 << " : " << fronts[i] << endl;
    // }

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

    // //For debugging purposes
    // std::cout << "Individuals added from complete fronts : " << new_population.size() << endl;

    if (i < fronts.size()) {
        // calculte crowding distance for the remaining front
        crowding_distance(population, fronts[i]); // crowding_distance(population, fronts[i])

        sort(fronts[i].begin(), fronts[i].end(), [&](int a, int b) {
            return crowd_compare(population[a], population[b]);
        }); // sort(fronts[i] by crowding_distance)

        size_t remaining = population_size - new_population.size();
        size_t available = fronts[i].size();
        for(size_t j = 0; j < min(remaining, available); ++j) {
            new_population.push_back(population[fronts[i][j]]);

            // //For debugging purposes
            // std::cout<< "Individual added from remaining front : " 
            //     << population[fronts[i][j]].objectives
            //     << endl;
            // std::cout<< "Rank : " << population[fronts[i][j]].rank << endl;
            // std::cout<< "Crowding Distance : " << population[fronts[i][j]].distance << endl;
            // std::cout<< "Iterator at : " << j << endl;

        } // new_population += fronts[i][:population_size - |new_population|]
    }

    // //For debugging purposes
    // std::cout << "Individuals added from remaining fronts : " << new_population.size() << endl;

    return new_population;
}

namespace py = pybind11;
using namespace py;

PYBIND11_MODULE(nsga2, m) {
    m.doc() = "NSGA-II algorithm";
    class_<individual>(m, "ind")
        .def(py::init<vector<double>>())
        .def_readwrite("obj", &individual::objectives)
        .def_readonly("cnt", &individual::domination_count)
        .def_readonly("dominated_solutions", &individual::dominated_solutions)
        .def_readonly("distance", &individual::distance)
        .def_readonly("rank", &individual::rank);
    m.def("nsga2", &nsga2, "Run the NSGA-II algorithm");
}


// //For debugging purposes
// int main() {
//     // Create an initial population
//     std::vector<individual> initial_population = {
//         individual({1.0, 2.0}),
//         individual({2.0, 1.0}),
//         individual({0.5, 1.5}),
//         individual({1.5, 1.5}),
//         individual({1.2, 1.8}),
//         individual({1.8, 1.2}),
//         individual({1.9, 1.1})
//     };

//     // Set the desired new population size
//     int new_population_size = 4;

//     // Run NSGA-II
//     std::vector<individual> result = nsga2(initial_population, new_population_size);

//     for (size_t i = 0; i < result.size(); ++i) {
//         std::cout << "Individual " << i + 1 << ":" << std::endl;
//         std::cout << "  Objectives: " << result[i].objectives << std::endl;
//         std::cout << "  Rank: " << result[i].rank << std::endl;
//         std::cout << "  Crowding Distance: " << result[i].distance << std::endl;
//     }

//     return 0;
// }