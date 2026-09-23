// Bottom-k minimizer Jaccard graph.
// Each node keeps up to --top-k neighbours with Jaccard >= --min-sim.
// Edges are written in both directions so a later resolver can walk either way.

#include <algorithm>
#include <cstdint>
#include <deque>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kK = 21;
constexpr int kWindow = 11;
constexpr int kSketch = 512;
constexpr uint64_t kMask = (kK == 32) ? ~0ULL : ((1ULL << (2 * kK)) - 1ULL);

int base_code(char ch) {
    switch (ch) {
        case 'A':
        case 'a':
            return 0;
        case 'C':
        case 'c':
            return 1;
        case 'G':
        case 'g':
            return 2;
        case 'T':
        case 't':
            return 3;
        default:
            return -1;
    }
}

uint64_t reverse_complement(uint64_t code) {
    uint64_t out = 0;
    for (int i = 0; i < kK; ++i) {
        out = (out << 2) | ((code & 3ULL) ^ 3ULL);
        code >>= 2;
    }
    return out;
}

std::vector<uint64_t> sketch_sequence(const std::string& sequence) {
    std::vector<uint64_t> minimizers;
    std::deque<std::pair<uint64_t, int>> window;
    uint64_t code = 0;
    int filled = 0;
    int index = 0;
    for (char ch : sequence) {
        const int base = base_code(ch);
        if (base < 0) {
            filled = 0;
            code = 0;
            window.clear();
            continue;
        }
        code = ((code << 2) | static_cast<uint64_t>(base)) & kMask;
        ++filled;
        if (filled < kK) {
            continue;
        }
        const uint64_t canonical = std::min(code, reverse_complement(code));
        while (!window.empty() && window.back().first > canonical) {
            window.pop_back();
        }
        window.emplace_back(canonical, index);
        while (!window.empty() && window.front().second <= index - kWindow) {
            window.pop_front();
        }
        if (index >= kWindow - 1) {
            minimizers.push_back(window.front().first);
        }
        ++index;
    }
    std::sort(minimizers.begin(), minimizers.end());
    minimizers.erase(std::unique(minimizers.begin(), minimizers.end()), minimizers.end());
    if (static_cast<int>(minimizers.size()) > kSketch) {
        minimizers.resize(kSketch);
    }
    return minimizers;
}

double jaccard(const std::vector<uint64_t>& left, const std::vector<uint64_t>& right) {
    if (left.empty() && right.empty()) {
        return 0.0;
    }
    size_t i = 0;
    size_t j = 0;
    size_t inter = 0;
    while (i < left.size() && j < right.size()) {
        if (left[i] == right[j]) {
            ++inter;
            ++i;
            ++j;
        } else if (left[i] < right[j]) {
            ++i;
        } else {
            ++j;
        }
    }
    const size_t uni = left.size() + right.size() - inter;
    if (uni == 0) {
        return 0.0;
    }
    return static_cast<double>(inter) / static_cast<double>(uni);
}

struct Sequence {
    std::string name;
    std::string bases;
};

std::vector<Sequence> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot read fasta");
    }
    std::vector<Sequence> records;
    std::string line;
    Sequence current;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            if (!current.name.empty()) {
                records.push_back(std::move(current));
                current = Sequence();
            }
            const auto end = line.find_first_of(" \t", 1);
            current.name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
            continue;
        }
        current.bases.append(line);
    }
    if (!current.name.empty()) {
        records.push_back(std::move(current));
    }
    return records;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "usage: kmer_knn FASTA EDGES.tsv [--top-k 8] [--min-sim 0.15]\n";
        return 2;
    }
    const std::string fasta = argv[1];
    const std::string edges_path = argv[2];
    int top_k = 8;
    double min_sim = 0.15;
    for (int i = 3; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--top-k" && i + 1 < argc) {
            top_k = std::stoi(argv[++i]);
        } else if (arg == "--min-sim" && i + 1 < argc) {
            min_sim = std::stod(argv[++i]);
        } else {
            std::cerr << "unknown argument: " << arg << "\n";
            return 2;
        }
    }
    if (top_k < 1) {
        std::cerr << "top-k must be >= 1\n";
        return 2;
    }
    std::vector<Sequence> records;
    try {
        records = read_fasta(fasta);
    } catch (const std::exception& exc) {
        std::cerr << exc.what() << "\n";
        return 1;
    }
    const int n = static_cast<int>(records.size());
    std::vector<std::vector<uint64_t>> sketches(n);
    std::vector<int> lengths(n);
    for (int i = 0; i < n; ++i) {
        lengths[i] = static_cast<int>(records[i].bases.size());
        sketches[i] = sketch_sequence(records[i].bases);
        records[i].bases.clear();
        records[i].bases.shrink_to_fit();
    }

    std::map<std::pair<int, int>, double> directed;
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<double, int>> scored;
        scored.reserve(n > 0 ? n - 1 : 0);
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            const double sim = jaccard(sketches[i], sketches[j]);
            if (sim >= min_sim) {
                scored.emplace_back(sim, j);
            }
        }
        const int keep = std::min(top_k, static_cast<int>(scored.size()));
        if (keep > 0) {
            std::partial_sort(
                scored.begin(),
                scored.begin() + keep,
                scored.end(),
                [](const auto& a, const auto& b) {
                    if (a.first != b.first) {
                        return a.first > b.first;
                    }
                    return a.second < b.second;
                });
        }
        for (int k = 0; k < keep; ++k) {
            const int j = scored[k].second;
            const double weight = scored[k].first;
            directed[{i, j}] = std::max(directed[{i, j}], weight);
            directed[{j, i}] = std::max(directed[{j, i}], weight);
        }
    }

    std::ofstream out(edges_path);
    std::ofstream lengths_out(edges_path + ".lengths.tsv");
    if (!out || !lengths_out) {
        std::cerr << "cannot write edges\n";
        return 1;
    }
    out << "edge_id\tsource\ttarget\torientation\tweight\n";
    lengths_out << "node_id\tlength\n";
    for (int i = 0; i < n; ++i) {
        lengths_out << records[i].name << '\t' << lengths[i] << '\n';
    }
    int edge_index = 1;
    for (const auto& item : directed) {
        const int source = item.first.first;
        const int target = item.first.second;
        out << "e" << edge_index++ << '\t' << records[source].name << '\t' << records[target].name
            << "\t++\t" << item.second << '\n';
    }
    std::cerr << "nodes\t" << n << "\tedges\t" << (edge_index - 1) << "\n";
    return 0;
}
