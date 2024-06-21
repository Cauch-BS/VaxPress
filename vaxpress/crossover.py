#add support for crossover
#add support for Reproduction using crossover
import numpy as np
from time import time

class BinaryCrossOver:
    def __init__(self, parent1: list[str] , parent2: list[str]):
        """Initializes a recombination object with a list of sequences of codons"""
        self.sequence1 = parent1[:]
        self.sequence2 = parent2[:]

        assert len(self.sequence1) == len(self.sequence2), "Recombined Sequences must have the same length"

        self.length = len(self.sequence1) 

    def single_point_crossover(self, crossover_prob: float) -> tuple[list]:
        assert crossover_prob < 0.5, "Crossover probability must be less than 0.5"
        self.child1, self.child2 = self.sequence1[:], self.sequence2[:]
        pos_crossover = max(1, min(np.random.binomial(self.length, crossover_prob), self.length - 1))
        #troubleshooting code
        print(f"Crossver position: {pos_crossover}")
        for i in range(pos_crossover):
            self.child1[i], self.child2[i] = self.sequence2[i], self.sequence1[i]
        return self.child1, self.child2

    def two_point_crossover(self, crossover_prob: float) -> tuple[list]:
        assert crossover_prob < 0.5, "Crossover probability must be less than 0.5"
        self.child1, self.child2 = self.sequence1[:], self.sequence2[:]
        start = np.random.randint(0, self.length)
        diff = max(1, min(np.random.binomial(self.length, crossover_prob), self.length - 1))
        if start + diff < self.length:
            end = start + diff
        elif start + diff >= self.length:
            if start > diff:
                start, end = start - diff, start
            else:
                if start < self.length // 2:
                    start, end = 0, start
                elif start >= self.length // 2:
                    start, end = start, self.length - 1
                else:
                    raise ValueError("Invalid start and end points")
        #troubleshooting code
        print(f"Crossver position: {start} and {end}")
        for i in range(start, end):
            self.child1[i], self.child2[i] = self.sequence2[i], self.sequence1[i]

        return self.child1, self.child2

    def binomial_crossover(self, crossover_prob: float) -> tuple[list]:
        assert crossover_prob < 0.5, "Crossover probability must be less than 0.5"
        self.child1, self.child2 = self.sequence1[:], self.sequence2[:]
        n_crossover = np.random.binomial(self.length, crossover_prob)
        n_crossover = max(1, min(n_crossover, self.length))

        crossover_choices = np.random.choice(self.length, n_crossover, replace = False)
        #troubleshooting code
        print(f"Crossver positions: {crossover_choices}")
        for i in crossover_choices:
            self.child1[i], self.child2[i] = self.sequence2[i], self.sequence1[i]

        return self.child1, self.child2

if __name__ == "__main__":
    #Run Tests for Binary Crossover
    print("Running Tests for Binary Crossover")
    parent1 = ['GGG']*5000
    parent2 = ['AAA']*5000
    crossover = BinaryCrossOver(parent1, parent2)
    current_time = time()
    print("Running Tests for Single Point Crossover...")
    #for tests in range(1000):
    single_point = crossover.single_point_crossover(0.2)
    print(f"Test successfully finished in {time() - current_time}")
    current_time = time()
    print("Running Tests for Dual Point Crossover...")
    #for tests in range(1000):
    two_point = crossover.two_point_crossover(0.2)
    print(f"Test successfully finished in {time() - current_time}")
    current_time = time()
    print("Running Tests for Binomial Crossover...")
    #for tests in range(1000):
    binomial_point = crossover.binomial_crossover(0.2)
    print(f"Test successfully finished in {time() - current_time}")
    print("Done")
