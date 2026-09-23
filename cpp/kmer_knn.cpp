// Symmetrized k-nearest-neighbour graph from canonical 4-mer composition.
// Disjoint contigs from one genome share tetranucleotide composition even when
// they share few 21-mers. Each node keeps up to --top-k neighbours whose cosine
// is at least --min-sim. Edges are written in both directions.

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kK = 4;
constexpr int kBins = 256;

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

int reverse_complement(int code) {
    int out = 0;
    for (int i = 0; i < kK; ++i) {
        out = (out << 2) | ((code & 3) ^ 3);
        code >>= 2;
    }
    return out;
}

struct Sequence {
    std::string name;
    std::vector<double> composition;
    int length = 0;
};

Sequence composition_of(std::string name, const std::string& sequence) {
    Sequence record;
    record.name = std::move(name);
    record.length = static_cast<int>(sequence.size());
    std::vector<double> counts(kBins, 0.0);
    int code = 0;
    int filled = 0;
    double total = 0.0;
    for (char ch : sequence) {
        const int base = base_code(ch);
        if (base < 0) {
            filled = 0;
            code = 0;
            continue;
        }
        code = ((code << 2) | base) & (kBins - 1);
        ++filled;
        if (filled < kK) {
            continue;
        }
        const int canonical = std::min(code, reverse_complement(code));
        counts[canonical] += 1.0;
        total += 1.0;
    }
    if (total > 0.0) {
        for (double& value : counts) {
            value /= total;
        }
    }
    record.composition = std::move(counts);
    return record;
}

double cosine(const std::vector<double>& left, const std::vector<double>& right) {
    double dot = 0.0;
    double left_norm = 0.0;
    double right_norm = 0.0;
    for (size_t i = 0; i < left.size(); ++i) {
        dot += left[i] * right[i];
        left_norm += left[i] * left[i];
        right_norm += right[i] * right[i];
    }
    if (left_norm <= 0.0 || right_norm <= 0.0) {
        return 0.0;
    }
    return dot / std::sqrt(left_norm * right_norm);
}

std::vector<Sequence> read_fasta(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot read fasta");
    }
    std::vector<Sequence> records;
    std::string line;
    std::string name;
    std::string bases;
    auto flush = [&]() {
        if (name.empty()) {
            return;
        }
        records.push_back(composition_of(name, bases));
        name.clear();
        bases.clear();
    };
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        if (line[0] == '>') {
            flush();
            const auto end = line.find_first_of(" \t", 1);
            name = line.substr(1, end == std::string::npos ? std::string::npos : end - 1);
            continue;
        }
        bases.append(line);
    }
    flush();
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
    std::map<std::pair<int, int>, double> directed;
    for (int i = 0; i < n; ++i) {
        std::vector<std::pair<double, int>> scored;
        scored.reserve(n > 0 ? static_cast<size_t>(n - 1) : 0);
        for (int j = 0; j < n; ++j) {
            if (i == j) {
                continue;
            }
            const double sim = cosine(records[i].composition, records[j].composition);
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
            const int target = scored[k].second;
            const double weight = scored[k].first;
            directed[{i, target}] = std::max(directed[{i, target}], weight);
            directed[{target, i}] = std::max(directed[{target, i}], weight);
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
    for (const Sequence& record : records) {
        lengths_out << record.name << '\t' << record.length << '\n';
    }
    int edge_index = 1;
    for (const auto& item : directed) {
        out << "e" << edge_index++ << '\t' << records[item.first.first].name << '\t'
            << records[item.first.second].name << "\t++\t" << item.second << '\n';
    }
    std::cerr << "nodes\t" << n << "\tedges\t" << (edge_index - 1) << "\n";
    return 0;
}
