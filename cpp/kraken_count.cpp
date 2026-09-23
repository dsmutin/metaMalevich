// Aggregate Kraken2 k-mer evidence.
// Each output row is one (sequence, taxon) pair summed from the k-mer column.
// Column 3 of the Kraken line is the classifier call and is not counted here.

#include <cstdint>
#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <utility>

namespace {

bool parse_pair(const std::string& token, std::string& taxid, uint64_t& count) {
    const auto colon = token.find(':');
    if (colon == std::string::npos || colon == 0 || colon + 1 >= token.size()) {
        return false;
    }
    taxid = token.substr(0, colon);
    const std::string raw = token.substr(colon + 1);
    if (taxid.empty() || raw.empty()) {
        return false;
    }
    count = 0;
    for (char ch : taxid) {
        if (ch < '0' || ch > '9') {
            return false;
        }
    }
    for (char ch : raw) {
        if (ch < '0' || ch > '9') {
            return false;
        }
        count = count * 10 + static_cast<uint64_t>(ch - '0');
    }
    return true;
}

}  // namespace

int main(int argc, char** argv) {
    if (argc != 3 && argc != 4) {
        std::cerr << "usage: kraken_count INPUT.output COUNTS.tsv [CALLS.tsv]\n";
        return 2;
    }
    std::ifstream in(argv[1]);
    if (!in) {
        std::cerr << "kraken_count: failed to open input\n";
        return 1;
    }
    std::map<std::pair<std::string, std::string>, uint64_t> totals;
    std::map<std::string, std::pair<std::string, std::string>> calls;
    std::string line;
    int line_number = 0;
    while (std::getline(in, line)) {
        ++line_number;
        if (line.empty() || line[0] == '#') {
            continue;
        }
        std::string fields[5];
        int nfields = 0;
        std::string current;
        for (size_t i = 0; i <= line.size() && nfields < 5; ++i) {
            if (i == line.size() || line[i] == '\t') {
                fields[nfields++] = current;
                current.clear();
            } else {
                current.push_back(line[i]);
            }
        }
        if (nfields < 5 || fields[1].empty()) {
            std::cerr << "kraken_count: line " << line_number << " does not have 5 columns\n";
            return 1;
        }
        if (calls.count(fields[1])) {
            std::cerr << "kraken_count: duplicate sequence id " << fields[1] << "\n";
            return 1;
        }
        calls[fields[1]] = {fields[0], fields[2]};
        const std::string& seq = fields[1];
        const std::string& kmers = fields[4];
        std::string token;
        for (size_t i = 0; i <= kmers.size(); ++i) {
            if (i == kmers.size() || kmers[i] == ' ') {
                if (!token.empty() && token != "|:|") {
                    std::string taxid;
                    uint64_t count = 0;
                    if (parse_pair(token, taxid, count)) {
                        totals[{seq, taxid}] += count;
                    }
                }
                token.clear();
            } else {
                token.push_back(kmers[i]);
            }
        }
    }
    std::ofstream out(argv[2]);
    if (!out) {
        std::cerr << "kraken_count: failed to open output\n";
        return 1;
    }
    out << "seq_id\ttaxon_id\tcount\n";
    for (const auto& item : totals) {
        out << item.first.first << '\t' << item.first.second << '\t' << item.second << '\n';
    }
    if (argc == 4) {
        std::ofstream call_out(argv[3]);
        if (!call_out) {
            std::cerr << "kraken_count: failed to open calls output\n";
            return 1;
        }
        call_out << "seq_id\tstatus\ttaxon_id\n";
        for (const auto& item : calls) {
            call_out << item.first << '\t' << item.second.first << '\t' << item.second.second << '\n';
        }
    }
    return 0;
}
