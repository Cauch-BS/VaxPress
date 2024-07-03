#include <iostream>
#include <string>
#include <unordered_map>
#include <vector>
#include <tuple>
#include <algorithm>

std::vector<std::tuple<std::string, int>> rabin_karp_repeated_substrings(const std::string& s, int min_length, int max_length) {
    int n = s.length();
    int truncate_n = std::min(n, max_length)  ;
    int base = 4;
    int mod = 2147483647;
    std::vector<std::tuple<std::string, int>> result;
    std::vector<std::vector<bool>> is_substring(n, std::vector<bool>(n, false));

    auto calculate_hash = [&](const std::string& sub) {
        long long h = 0;
        for (char c : sub) {
            h = (h * base + c) % mod;
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
        
        long long new_hash = (prev_hash * base - left_char * exp_result + right_char) % mod;
        new_hash = new_hash < 0 ? new_hash + mod: new_hash;
        return new_hash;
    };

    for (int length = truncate_n - 1; length >= min_length - 1; --length) {
        std::unordered_map <long long, int> seen;
        long long current_hash = calculate_hash(s.substr(0, length + 1));
        seen[current_hash] = 0;

        for (int i = 1; i <= n - length - 1; ++i) {
            long long exp_result = power_mod(base, length + 1, mod);
            current_hash = roll_hash(current_hash, s[i - 1], s[i + length], exp_result);
            std::string substring = s.substr(i, length + 1);
            
            if (seen.find(current_hash) != seen.end()) {
                std::string substring = s.substr(i, length + 1);
                int prev_pos = seen[current_hash];
                if (s.substr(prev_pos, length + 1) == substring){
                    result.push_back(std::make_tuple(substring, i));
                    std::cout << "Repeated subtring " << substring << " at position " << i << std::endl;
                }
            }else{
            seen[current_hash] = i;
            }
        }
    }

    return result;
}

int main() {
    std::string input_sequence = "CAGAAATCAATTCTTTCAGAAATCAATTGGTACCTTACTGAATTATCGATTTTCTGTTTTCGTCCTACAAATACTTTAATGGGGTGCGGCAGGTAGTTATTGCCCATTGTACTCAGCAGAAATCAATTTACCGGTGTCGCGGCGTAGGCCAAGCCCCAACATAGGATCTTCCTT";
    int min_length = 8;
    int max_length = 100;
    std::vector<std::tuple<std::string, int>> result = rabin_karp_repeated_substrings(input_sequence, min_length, max_length);

    for (const auto& tuple : result) {
        std::cout << "(" << std::get<0>(tuple) << ", " << std::get<1>(tuple) << ")" << std::endl;
    }

    return 0;
}
