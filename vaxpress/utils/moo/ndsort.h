/* An implementation of the fast non-dominated sorting algorithm
the algorithm has O(M*N^2) complexity where M is the number of parameters
and N is the size of the population. The algorithm is based on the
original paper by Deb et al. (2002) titled 
"A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"*/
#ifndef NDSORT_H
#define NDSORT_H

#include <vector>
#include "individual.h"

using namespace std;

bool dominates(const individual& a, const individual& b) {
    // complexity: O(M) where M is the number of objectives
    bool one_better = false; // assume a dominates b
    for (int i = 0; i < a.objectives.size(); ++i) {
        if (a.objectives[i] < b.objectives[i]) {
            return false; // a does not dominate b, return false
        } else if (a.objectives[i] > b.objectives[i]) {
            one_better = true; // a dominates b in at least one objective
        }
    }
    return one_better;
}

vector<vector<int>> non_dominated_argsort(vector<individual>& population) {
    // fronts is a list of indices of individuals in the population
    // fronts[0] is the first front, fronts[1] is the second front, and so on
    vector<vector<int>> fronts(1); 

    // generate front based on non-domination
    // complexity: O(M * N^2) where N is the size of the population and M is the number of objectives
    // space complexity: O(N^2)
    /* pseudocode equivalent
    for each p in population:
        for each q in population:
            if p dominates q:
                add q to p's dominated_solutions
            else if q dominates p:
                increment p's domination_count
        if p's domination_count == 0:
            p's rank = 1
            add p to front(0)
    */

    for (int i = 0; i < population.size(); ++i) { // for each p in population
        for (int j = 0; j < population.size(); ++j) { // for each q in population
            if (i == j) continue; // skip if p == q

            bool i_dominates_j = dominates(population[i], population[j]);
            bool j_dominates_i = dominates(population[j], population[i]);

            if (i_dominates_j && !j_dominates_i) { // if p dominates q
                population[i].dominated_solutions.push_back(j); // add q to p's dominated_solutions
            } else if (j_dominates_i && !i_dominates_j) { // if q dominates p
                population[i].domination_count++; // increment p's domination_count
            }
        }
        // generate front of non-dominated individuals (these individuas have rank 1)
        // if p's domination_count == 0, add p to front(0)
        if (population[i].domination_count == 0) {
                population[i].rank = 1; // p's rank = 1
                fronts[0].push_back(i); // add p to front(0)
        }
    }

    // generate subsequent fronts
    // complexity: O(fronts.size() * fronts[max].size()) 
    // This is less than O(N^2) where N is the size of the population
    // N = sum(fronts[i].size() for i in range(fronts.size())
    /* pseudocode equivalent
    i = 0
    while front(i) is not empty:
        next_front = []
        for each p in front(i):
            for each q in p's dominated_solutions:
                decrement q's domination_count
                if q's domination_count == 0:
                    q's rank = i + 1
                    add q to next_front
        i++
        add next_front to front(i)
    */
    int i = 0;
    while (!fronts[i].empty()) {
        vector<int> next_front;
        for (int j = 0; j < fronts[i].size(); ++j) { 
            // for each p in front(i)
            for (int k = 0; k < population[fronts[i][j]].dominated_solutions.size(); ++k) {
                int q = population[fronts[i][j]].dominated_solutions[k]; // for each q in p's dominated_solutions
                population[q].domination_count--; // decrement q's domination_count
                if (population[q].domination_count == 0) { // if q's domination_count == 0
                    population[q].rank = i + 1; // q's rank = i + 1
                    next_front.push_back(q); // add q to next_front
                }
            }
        }
        i++;
        fronts.push_back(next_front);
    }

    return fronts;
};

#endif // NDSORT_H
