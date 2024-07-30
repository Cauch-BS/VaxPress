import nsga2
import time
from concurrent import futures

def run_moo(self):
    self.show_configuration()

    timelogs = [time.time()]
    n_survivors = self.execopts.n_survivors
    last_winddown = 0
    error_code = 0

    with futures.ProcessPoolExecutor(max_workers=self.n_processes) as executor:
        for i in range(self.execopts.n_iterations):
            iter_no = i + 1
            n_parents = len(self.population)

            try:
                self.next_generation(i)
            except StopIteration:
                break

            total_scores, scores, metrics, foldings, pairingprobs = self.seqeval.evaluate_moo(
                                                self.flatten_seqs, executor)
            
            for score in total_scores:
                if score is None:
                    # Termination due to errors from one or more scoring functions
                    error_code = 1
                    break
            total_scores = [nsga2.ind(scores) for scores in total_scores]
            survivor_indices = nsga2.nsga2(total_scores, n_survivors)
            survivors = [self.population[i] for i in survivor_indices]
            survivor_foldings = [foldings[i] for i in survivor_indices]
            survivor_bpps = [pairingprobs[i] for i in survivor_indices]  # noqa: F841
            self.best_scores.append(total_scores[survivor_indices[0]])
    return 
