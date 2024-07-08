/* An implementation of the crowd distance assignment algorithm
the algorithm has O(MNlog(N)) complexity where M is the number of parameters
and N is the size of the population. The algorithm is based on the
original paper by Deb et al. (2002) titled 
"A Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II"*/
#ifndef CROWDING_H
#define CROWDING_H

#include <vector>
#include <algorithm>
#include <limits>
#include "individual.h"

using namespace std;

void crowding_distance(vector<individual>& population, 
                        vector<int>& front){
    /*pseudocode equivalent
    l = |front|
    for each i in front:
        crowding_distance[i] = 0
    for each m in objectives:
        front = sort(front by objective [m])
        crowding_distance[front[0]] = crowding_distance[front[l-1]] = infinity
        for i in range(1, l-1):
            crowding_distance[front[i]] += (
                population[front[i+1]].objectives[m] - 
                population[front[i-1]].objectives[m]
            ) / (max(objectives[m]) - min(objectives[m])
    */
    int n = front.size(); // l = |front|
    // complexity O(M * Nlog(N))
    for (int m = 0; m < population[0].objectives.size(); ++m) { // for each m in objectives
        sort(front.begin(), front.end(), [&](int a, int b) {
            return population[a].objectives[m] < population[b].objectives[m];
        }); // sort front by objective [m] 
        //complexity O(Nlog(N))
        // (python equivalent: front.sort(key=lambda x: population[x].objectives[m]) )
        population[front[0]].distance = population[front[n-1]].distance = numeric_limits<double>::infinity(); 
        // crowding_distance[front[0]] = crowding_distance[front[l-1]] = infinity
        for (int i = 1; i < n-1; ++i) { // for i in range(1, l-1)
            population[front[i]].distance += (population[front[i+1]].objectives[m] - population[front[i-1]].objectives[m]) / 
            // divide by normalizating constant
            (population[front[n-1]].objectives[m] // max(objectives[m])
            - population[front[0]].objectives[m]); // min(objectives[m])
        }
        /* pseudocode equivalent is
        crowding_distance[front[i]] += (
            population[front[i+1]].objectives[m] - 
            population[front[i-1]].objectives[m]
        ) / (max(objectives[m]) - min(objectives[m])
        */
    }
}

/* after all population members of the front are assigned a crowding distance, 
sort the front by crowding distance. A solution with a smaller value of crowding
distance is i some sense more similar to other solutins (i.e. more crowded)
To generate a pareto front, we need to select the solutions which are spaced far apart,
i.e. the solutions with the largest crowding distance
*/

bool crowd_compare(const individual& a, const individual& b) {
    bool prefer_a_to_b = false;
    if (a. rank < b. rank){
        prefer_a_to_b = true;
    } else if (a.rank == b.rank) {
        if (a.distance > b.distance){
            prefer_a_to_b = true;
        }
    }
    return prefer_a_to_b;
}

#endif // CROWDING_H