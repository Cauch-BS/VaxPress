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
                result = (result * static_cast<long long>(base)) % modulus;
            }
            base = (static_cast<long long>(base) * base) % modulus;
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
            cout << "Repeated subtring " << substring << " at position " << i << endl;
            
        }else{
        seen[current_hash] = i;
        }
    }
    

    return result;
}

vector<tuple<string, int, int, int, int>> groupConsecutive(const vector<tuple<string, int, int>>& input, int min_length) {
    vector<tuple<string, int, int, int, int>> result;

    // Start with the first tuple
    auto current = input[0];
    string current_str = get<0>(current);
    int last = get<1>(current);
    int last_end = last;
    int now = get<2>(current);
    int end = now;

    for (size_t i = 1; i < input.size(); ++i) {
        if (get<2>(input[i]) == end + 1) {
            // Extend the current range
            current_str += get<0>(input[i])[min_length - 1];
            end = get<2>(input[i]);
        } else {
            last_end = last + end - now + min_length - 1;
            end = end + min_length - 1;
            // Add the current group to the result
            result.emplace_back(current_str, last, last_end, now, end);
            // Start a new group
            current_str = get<0>(input[i]); 
            auto rbegin = make_reverse_iterator(input.begin() + i);
            auto rend = make_reverse_iterator(input.begin());

            auto last_it = find_if(
                        rbegin, rend, 
                        [&current_str](const tuple<string, int, int>& t) {
                            return get<0>(t) == current_str;
                        });

            last = last_it != rend ? - distance(last_it, rend) : get<1>(input[i]) ;
            last_end = last;
            now = get<2>(input[i]);
            end = get<2>(input[i]);
        }
    }

    last_end = last + end - now + min_length - 1;
    end = end + min_length - 1;
    // Add the last group to the result
    result.emplace_back(current_str, last, last_end, now, end);

    return result;
}


int main() {
    string input_sequence = "CAGAAATCAATTCTTTCAGAAATCAATTGGTACCTTACTGAATTATCGATTTTCTGTTTTCGTCCTACAAATACTTTAATGGGGTGCGGCAGGTAGTTATTGCCCATTGTACTCAGCAGAAATCAATTTACCGGTGTCGCGGCGTAGGCCAAGCCCCAACATAGGATCTTCCTT";
    int min_length = 8;
    vector<tuple<string, int, int, int, int>> result = groupConsecutive(
        rabin_karp_repeated_substrings(
            input_sequence, min_length
        ), min_length
    );

    for (const auto& tuple : result) {
        cout << "(" << get<0>(tuple) 
             << ", " << get<1>(tuple) 
             << ", " << get<2>(tuple) 
             << ", " << get<3>(tuple) 
             << ", " << get<4>(tuple) 
             << ")" << endl;
    }

    return 0;
}
