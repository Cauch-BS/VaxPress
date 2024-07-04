#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <algorithm>
#define GET_ACGU_NUM(x) ((x=='A'? 0 : (x=='C'? 1 : (x=='G'? 2 : (x=='U'?3: 4)))))

using namespace std;

vector<tuple<string, int, int>> rabin_karp_repeated_substrings(const string& s, int min_length) {
    int n = s.length();
    int base = 4;
    int mod = 2147483647;
    vector<tuple<string, int, int>> result;
    vector<vector<bool>> is_substring(n, vector<bool>(n, false));

    auto calculate_hash = [&](const string& sub) {
        long long h = 0;
        for (char c : sub) {
            h = (h * base + GET_ACGU_NUM(c)) % mod;
        }
        return h;
    };

    auto power_mod = [&](int base, int exp, int modulus) {
        long long result = 1;
        base = base % modulus;
        while (exp > 0) {
            if (exp % 2 == 1) {
                result = (result * base) % modulus;
            }
            base = (base * base) % modulus;
            exp /= 2;
        }

        return result;
    };
    
    auto roll_hash = [&](long long prev_hash, char left_char, char right_char, long long exp_result) {
        long long new_hash = (prev_hash * base - GET_ACGU_NUM(left_char) * exp_result + GET_ACGU_NUM(right_char)) % mod;
        new_hash = new_hash < 0 ? new_hash + mod: new_hash;
        return new_hash;
    };

    unordered_map <long long, int> seen;
    long long current_hash = calculate_hash(s.substr(0, min_length));
    seen[current_hash] = 0;

    for (int i = 1; i < n - min_length; ++i) {
        long long exp_result = power_mod(base, min_length, mod);
        current_hash = roll_hash(current_hash, s[i - 1], s[i + min_length - 1], exp_result);
        
        if (seen.find(current_hash) != seen.end()) {
            string substring = s.substr(i, min_length);
            int prev_pos = seen[current_hash];
            result.emplace_back(substring, prev_pos, i);
            seen[current_hash] = i;            
        }else{
        seen[current_hash] = i;
        }
    }
    

    return result;
}

float returnRepeatsPenalty(const string& seq, int min_length){
    int length = seq.length();
    vector<tuple<string, int, int>> repeats = rabin_karp_repeated_substrings(seq, min_length);
    float penalty = 0;
    float epsilon = 1.0 / 256;
    for (const auto& tuple : repeats) {
        float score = 8/(get<2>(tuple)-get<1>(tuple));
        float norm_score = (1.0 + epsilon) /(1.0 + epsilon  - score);
        penalty += norm_score;
        cout << "Score: " << norm_score << endl;
    }
    penalty /= length;
    return penalty;
}


int main() {
    string input_sequence = "CAGAAATCAATTCTTTCAGAAATCAATTGGTACCTTACTGAATTATCGATTTTCTGTTTTCGTCCTACAAATACTTTAATGGGGTGCGGCAGGTAGTTATTGCCCATTGTACTCAGCAGAAATCAATTTACCGGTGTCGCGGCGTAGGCCAAGCCCCAACATAGGATCTTCCTT";
    int min_length = 8;
    
    vector <tuple<string, int, int>> repeats = rabin_karp_repeated_substrings(input_sequence, min_length);
    
    float penalty = returnRepeatsPenalty(input_sequence, min_length);

    for (const auto& tuple : repeats) {
        cout << "Substring: " << get<0>(tuple) << " Last Seen: " << get<1>(tuple) << " Now: " << get<2>(tuple) << endl;
    }

    cout << "Repeats penalty: " << penalty << endl;

    return 0;
}


