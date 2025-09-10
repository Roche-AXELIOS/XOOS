#include "pangenome-consensus-caller-core.h"

#include <getopt.h>
#include <unistd.h>

#include <algorithm>
#include <atomic>
#include <cassert>
#include <cctype>
#include <cstring>
#include <fstream>
#include <iostream>
#include <list>
#include <mutex>
#include <sstream>
#include <string>
#include <thread>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <htslib/hts.h>
#include <htslib/sam.h>

#include <xoos/yc-decode/yc-decoder.h>

namespace xoos::pangenome_consensus_caller {

uint8_t ADJUSTED_BQ = 22;

const char SIMPLEX = '~';
const char MATCH = '*';

// the inserted character for insertion and deletion YC operation, else null
static const char INDELS[256] = {
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', 'A',  'G',  '\0',
    'A',  '\0', '\0', 'G',  'C',  'C',  '\0', '\0', '\0', '\0', '\0', '\0', 'T',  '\0', 'T',  '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};

// the two characters involved in a substitution YC operation (at 2*c and
// 2*c+1), else null
static const char SUBSTITUTIONS[512] = {
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', 'C',
    'A',  '\0', '\0', 'G',  'A',  'T',  'C',  'T',  'G',  '\0', '\0', 'T',  'A',  '\0', '\0', '\0', '\0', 'G',  'T',
    '\0', '\0', 'A',  'A',  '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', 'A',  'G',  'C',  'G',  '\0', '\0', '\0',
    '\0', 'G',  'C',  'A',  'T',  '\0', '\0', 'C',  'T',  '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0',
    '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0', '\0'};

// true for MDX=N
// TODO: if I left CIGARs in the HTSlib-native internal encoding, I could use
// its lookup tables for this
static const bool op_consumes_ref[256] = {
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, true,  false, false,
    false, false, false, false, true,  false, false, false, false, false, false, false, false, true,  true,  false,
    false, false, false, false, false, false, false, false, true,  false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false,
    false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false};

// for fixed-length SAM tag codes, the number of bytes the type occupies, else 0
// TODO: this used to be available as bam_aux_type2size in htslib, but it seems
// to have been removed...
static const uint8_t data_type_size[256] = {
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 4,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 4, 0, 0, 4, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

// complement base for ACGT- and the mismatch YC codes
static const char complement[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', '-', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'T',
    'K', 'G', 'Y', 'R', 'M', 'C', 'W', 'N', 'N', 'B', 'N', 'F', 'N', 'N', 'N', 'N', 'E', 'V', 'A', 'N', 'S', 'H',
    'N', 'D', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'};

void print_help() {
  std::cerr << "Usage:\n";
  std::cerr << "pangenome_consensus_caller [-t <threads>] [-b <excluded regions BED>] [-q <adjusted base qual>] "
               "alns.bam\n";
  std::cerr << "  -- or --\n";
  std::cerr << "pangenome_consensus_caller [-t <threads>] [-b <excluded regions BED>] [-q <adjusted base qual>] < "
               "alns.bam\n";
}

// return BED as vector of tuples of (contig ID from BAM header, pos begin, pos
// end) in sorted order
std::vector<std::tuple<int, size_t, size_t>> parse_bed(std::istream& in, bam_hdr_t* header) {
  std::vector<std::tuple<int, size_t, size_t>> parsed;

  std::string line;
  while (std::getline(in, line)) {
    std::stringstream line_strm(line);

    std::string contig;
    std::getline(line_strm, contig, '\t');
    int tid = bam_name2id(header, contig.c_str());
    if (tid >= 0) {
      std::string pos_begin;
      std::getline(line_strm, pos_begin, '\t');
      std::string pos_end;
      std::getline(line_strm, pos_end, '\t');

      parsed.emplace_back(bam_name2id(header, contig.c_str()), stoul(pos_begin), stoul(pos_end));
    }
  }

  if (!std::is_sorted(parsed.begin(), parsed.end())) {
    std::sort(parsed.begin(), parsed.end());
  }

  return parsed;
}

// return the graph character aligned to each position in a read, or '-' if the
// read alignment is an insertion. note: deletions are ignored
std::string read_aligned_chars(const std::string& sequence, const std::string& cs_string, bool rev) {
  std::string aligned_chars(sequence.size(), '\0');
  // if the read is on the reverse strand, we will traverse the sequence in
  // reverse as we go along the cs string
  size_t seq_idx = rev ? sequence.size() - 1 : 0;
  int64_t seq_incr = rev ? -1 : 1;
  for (size_t i = 0; i < cs_string.size();) {
    char op = cs_string[i++];
    // TODO: can't compile as switch-case because of nontrivial iteration on i
    if (op == ':') {
      // a match
      assert(isdigit(cs_string[i]));
      size_t j = i + 1;
      while (j < cs_string.size() && isdigit(cs_string[j])) {
        ++j;
      }
      size_t match_len = std::stoul(cs_string.substr(i, j - i));
      for (size_t k = 0; k < match_len; ++k, seq_idx += seq_incr) {
        aligned_chars[seq_idx] = sequence[seq_idx];
      }
      i = j;
    } else if (op == '*') {
      // a mismatch
      aligned_chars[seq_idx] = rev ? complement[(uint8_t)cs_string[i]] : cs_string[i];
      seq_idx += seq_incr;
      i += 2;
      // FIXME: this alternate value is needed for earlier data before the
      // Giraffe PR i += 3;
    } else if (op == '+') {
      // read is insertion
      while (isalpha(cs_string[i])) {
        aligned_chars[seq_idx] = '-';
        seq_idx += seq_incr;
        ++i;
      }
    } else if (op == '-') {
      // read is deletion
      while (isalpha(cs_string[i])) {
        ++i;
      }
    } else {
      std::cerr << "error: difference string has unrecognized operation " << op << "\n";
      exit(1);
    }
  }

  return aligned_chars;
}

// expand YC tag string to one operation per read base
std::string expand_yc(const yc_decode::YcTag& yc_tag) {
  std::string expanded(yc_tag.left_overhang, SIMPLEX);
  for (const auto& op : yc_tag.duplex_ops) {
    if (op.is_numeric) {
      for (size_t k = 0; k < op.length; ++k) {
        expanded.push_back(MATCH);
      }
    } else {
      expanded.push_back(op.code);
    }
  }
  expanded.reserve(expanded.size() + yc_tag.right_overhang);
  for (size_t i = 0; i < yc_tag.right_overhang; ++i) {
    expanded.push_back(SIMPLEX);
  }
  std::string(yc_tag.left_overhang, SIMPLEX);
  return expanded;
}

// convert an expanded YC tag into a YcTag struct
yc_decode::YcTag compress_yc(const std::string& expanded_yc_string, const yc_decode::YcTag& original_yc_tag) {
  yc_decode::YcTag yc_tag;
  yc_tag.left_overhang = original_yc_tag.left_overhang;
  yc_tag.right_overhang = original_yc_tag.right_overhang;
  yc_tag.left_overhang_is_r1 = original_yc_tag.left_overhang_is_r1;
  yc_tag.right_overhang_is_r1 = original_yc_tag.right_overhang_is_r1;

  size_t curr_match_len = 0;
  for (size_t i = yc_tag.left_overhang, n = expanded_yc_string.size() - yc_tag.right_overhang; i < n; ++i) {
    if (expanded_yc_string[i] == MATCH) {
      ++curr_match_len;
    } else {
      if (curr_match_len) {
        yc_tag.duplex_ops.emplace_back();
        auto& op = yc_tag.duplex_ops.back();
        op.is_numeric = true;
        op.length = curr_match_len;
        curr_match_len = 0;
      }
      yc_tag.duplex_ops.emplace_back();
      auto& op = yc_tag.duplex_ops.back();
      op.code = expanded_yc_string[i];
    }
  }
  if (curr_match_len) {
    yc_tag.duplex_ops.emplace_back();
    auto& op = yc_tag.duplex_ops.back();
    op.is_numeric = true;
    op.length = curr_match_len;
  }

  return yc_tag;
}

// un-RLE a CIGAR string
std::string expand_cigar(const bam1_t* read) {
  std::string expanded;
  expanded.reserve(read->core.l_qseq);

  uint32_t* cigar = bam_get_cigar(read);
  for (int i = 0, n = read->core.n_cigar; i < n; ++i) {
    size_t len = bam_cigar_oplen(cigar[i]);
    char op = bam_cigar_opchr(bam_cigar_op(cigar[i]));

    for (size_t j = 0; j < len; ++j) {
      expanded.push_back(op);
    }
  }

  return expanded;
}

// put an expanded CIGAR string into the HTSlib internal encoding
std::vector<uint32_t> encode_expanded_cigar(const std::string& expanded_cigar) {
  std::vector<uint32_t> cigar;

  char curr_op = '\0';
  uint32_t curr_op_len = 0;
  for (char op : expanded_cigar) {
    if (op == curr_op) {
      ++curr_op_len;
    } else {
      if (curr_op_len != 0) {
        cigar.push_back(bam_cigar_gen(curr_op_len, bam_cigar_table[(uint8_t)curr_op]));
      }
      curr_op = op;
      curr_op_len = 1;
    }
  }
  if (curr_op_len != 0) {
    cigar.push_back(bam_cigar_gen(curr_op_len, bam_cigar_table[(uint8_t)curr_op]));
  }
  return cigar;
}

// convert BAM sequence into string
std::string get_bam_sequence(const bam1_t* read) {
  std::string seq(read->core.l_qseq, '\0');
  auto bam_seq = bam_get_seq(read);
  for (size_t i = 0; i < seq.size(); ++i) {
    seq[i] = seq_nt16_str[bam_seqi(bam_seq, i)];
  }
  return seq;
}

// within homopolymers, reposition aligned inserts onto discordant inserts
// TODO: there are probably more general ways to do this in simple repeats that
// are not homopolymers as well
void reposition_inserts(std::string& aligned_chars, const std::string& sequence, const std::string& expanded_yc) {
  for (size_t i = 0; i < sequence.size();) {
    // identify homopolymer run
    size_t j = i + 1;
    while (j < sequence.size() && sequence[j] == sequence[i]) {
      ++j;
    }

    // locate inserts in the graph alignment and the yc string
    std::vector<size_t> yc_inserts;
    std::vector<size_t> aln_inserts;
    for (size_t k = i; k < j; ++k) {
      if (INDELS[(uint8_t)expanded_yc[k]] != '\0') {
        yc_inserts.push_back(k);
      }
      if (aligned_chars[k] == '-') {
        aln_inserts.push_back(k);
      }
    }

    // reposition them on top of each other
    for (size_t k = 0, n = std::min(yc_inserts.size(), aln_inserts.size()); k < n; ++k) {
      std::swap(aligned_chars[yc_inserts[k]], aligned_chars[aln_inserts[k]]);
    }

    i = j;
  }
}

// get the index of the interval that contains the contig and position
// assumes that intervals are sorted and non-overlapping
size_t locate_in_intervals(int contig_id, size_t pos, const std::vector<std::tuple<int, size_t, size_t>>& intervals) {
  size_t lo = 0;
  size_t hi = intervals.size();
  while (lo < hi) {
    size_t mid = (lo + hi) / 2;
    const auto& mid_interval = intervals[mid];
    if (std::get<0>(mid_interval) == contig_id) {
      if (std::get<1>(mid_interval) > pos) {
        hi = mid;
      } else if (std::get<2>(mid_interval) <= pos) {
        lo = mid + 1;
      } else {
        // we found the interval that contains the position
        return mid;
      }
    } else if (std::get<0>(mid_interval) > contig_id) {
      hi = mid;
    } else {
      lo = mid + 1;
    }
  }
  // no interval contains the position
  return lo;
}

// make a new BAM record with the modified fields
bam1_t* reencode_read(const bam1_t* read,
                      const std::string& mod_seq,
                      const std::vector<uint8_t>& mod_qual,
                      const std::string& mod_yc,
                      const std::string& mod_cigar,
                      size_t position_shift,
                      const yc_decode::YcTag& original_yc_tag) {
  yc_decode::YcTag compressed_yc_tag = compress_yc(mod_yc, original_yc_tag);
  if (bam_is_rev(read)) {
    // YC strings are reported relative to the forward strand of the read, not
    // the ref
    compressed_yc_tag.ReverseComplement();
  }
  auto compressed_yc_string = compressed_yc_tag.ToString();
  auto encoded_cigar = encode_expanded_cigar(mod_cigar);

  const auto& core = read->core;
  bam1_t* mod_read = bam_init1();
  bam_set1(mod_read,
           strlen(bam_get_qname(read)),
           bam_get_qname(read),
           core.flag,
           core.tid,
           core.pos + position_shift,
           core.qual,
           encoded_cigar.size(),
           encoded_cigar.data(),
           core.mtid,
           core.mpos,
           core.isize,
           mod_seq.size(),
           mod_seq.c_str(),
           (const char*)mod_qual.data(),
           // the aux will only shrink, so this should be enough memory
           bam_get_l_aux(read));
  mod_read->id = read->id;

  // update the bin assignment (modified from cram/cram_samtools.h)
  // TODO: there has to be a less opaque interface somewhere to do this...
  mod_read->core.bin = hts_reg2bin(mod_read->core.pos, bam_endpos(mod_read) - 1, 14, 5);

  // add the updated YC string as a new tag
  bam_aux_append(mod_read, "YC", 'Z', compressed_yc_string.size() + 1, (const uint8_t*)compressed_yc_string.c_str());

  // copy over the remaining tags
  uint8_t* aux_cursor = bam_get_aux(read);
  uint8_t* aux_end = read->data + read->l_data;
  while (aux_cursor != aux_end) {
    const char* tag = (const char*)aux_cursor;
    uint8_t tag_type = aux_cursor[2];
    size_t value_len = 0;
    switch (tag_type) {
      case 'A':
      case 'c':
      case 'C':
      case 's':
      case 'S':
      case 'i':
      case 'I':
      case 'f':
        // fixed length types
        value_len = data_type_size[tag_type];
        break;
      case 'Z':
      case 'H':
        // string (length includes null terminator)
        value_len = strlen(reinterpret_cast<char*>(aux_cursor + 3)) + 1;
        break;
      case 'B':
        // array (length include 1-byte type and 4-byte length)
        value_len = 5 + data_type_size[aux_cursor[3]] * bam_auxB_len(aux_cursor);
        break;
      default:
        std::cerr << "error: invalid SAM tag\n";
        exit(1);
        break;
    }

    if (!(aux_cursor[0] == 'Y' && aux_cursor[1] == 'C') && !(aux_cursor[0] == 'G' && aux_cursor[1] == 'R')) {
      // this is some other tag besides the YC (which was already written) and
      // the graph alignment (which may no longer describe the read)

      // copy it over
      bam_aux_append(mod_read, tag, tag_type, value_len, aux_cursor + 3);
    }

    // move to the next aux value
    aux_cursor += 3 + value_len;
  }

  return mod_read;
}

// main algorithm function in which we resolve discordancies using the
// pangenome-space alignment
bam1_t* clip_discordant(const bam1_t* read, const std::vector<std::tuple<int, size_t, size_t>>& exclusion_bed) {
  // extract tags
  uint8_t* YC_aux = bam_aux_get(read, "YC");
  uint8_t* GR_aux = bam_aux_get(read, "GR");
  if (YC_aux == nullptr || GR_aux == nullptr || GR_aux[1] == '*') {
    // need both tags to be able to execute the clipping algorithm, just copy
    // and return
    return bam_copy1(bam_init1(), read);
  }
  std::string YC_str = bam_aux2Z(YC_aux);
  std::string GR_str = bam_aux2Z(GR_aux);

  // translate GR tag into bases aligned at each read position
  std::string seq = get_bam_sequence(read);
  const uint8_t* qual = bam_get_qual(read);
  std::string aligned_chars = read_aligned_chars(seq, GR_str, bam_is_rev(read));

  // figure out the bounds of any softclips in the graph alignment
  size_t graph_nonclip_begin = 0;
  size_t graph_nonclip_end = aligned_chars.size();
  while (graph_nonclip_begin < aligned_chars.size() && aligned_chars[graph_nonclip_begin] == '-') {
    ++graph_nonclip_begin;
  }
  while (graph_nonclip_end != 0 && aligned_chars[graph_nonclip_end - 1] == '-') {
    --graph_nonclip_end;
  }

  std::string expanded_cigar = expand_cigar(read);

  yc_decode::YcTag parsed_yc = yc_decode::DeserializeYcTag(YC_str);
  if (bam_is_rev(read)) {
    parsed_yc.ReverseComplement();
  }
  std::string expanded_yc = expand_yc(parsed_yc);

  // try to locate inserts over discordant bases
  reposition_inserts(aligned_chars, seq, expanded_yc);

  // cursors that maintain our position in the BED intervals
  size_t ref_cursor = read->core.pos;
  size_t bed_cursor = locate_in_intervals(read->core.tid, ref_cursor, exclusion_bed);

  // the fields that we will edit and modify
  std::string mod_seq;
  std::vector<uint8_t> mod_qual;
  std::string mod_yc;
  std::string mod_cigar;
  mod_seq.reserve(seq.size());
  mod_qual.reserve(seq.size());
  mod_yc.reserve(seq.size());
  mod_cigar.reserve(expanded_cigar.size());

  // trackers that we may need if the alignment adjusts the reference position
  bool have_consumed_ref = false;
  size_t first_removed_run = 0;

  // FIXME: the outer iteration should really be over the graph alignment, not
  // the cigar, otherwise we won't correct reads whose alignment has been lost
  // due to being off-reference

  size_t read_idx = 0;
  for (size_t i = 0; i < expanded_cigar.size(); ++i) {
    bool do_adjudication = true;
    if ((bed_cursor < exclusion_bed.size() && std::get<0>(exclusion_bed[bed_cursor]) == read->core.tid &&
         ref_cursor >= std::get<1>(exclusion_bed[bed_cursor]) && ref_cursor < std::get<2>(exclusion_bed[bed_cursor])) ||
        read_idx < graph_nonclip_begin || read_idx >= graph_nonclip_end) {
      // we are in a BED region that we are excluding from pangenome-based
      // discordancy adjudication or this is a softclip in the graph alignment,
      // in which case we take the graph alignment to be uninformative
      do_adjudication = false;
    }

    char cigar_op = expanded_cigar[i];
    if (cigar_op == 'D' || cigar_op == 'N') {
      // operation does not correspond to a base in the read
      mod_cigar.push_back(cigar_op);
    } else {
      // operation corresponds to a base in the read (and also quals and YC
      // string)
      uint8_t bq = qual[read_idx];
      char base = seq[read_idx];
      char yc_op = expanded_yc[read_idx];
      bool include_base = true;
      if (do_adjudication && yc_op != MATCH && yc_op != SIMPLEX) {
        if (seq[read_idx] == aligned_chars[read_idx]) {
          // the discordant base matches the pangenome, bump the base qual up to
          // simplex level
          bq = ADJUSTED_BQ;
        } else if (aligned_chars[read_idx] == '-' && INDELS[(uint8_t)yc_op] != '\0') {
          // the pangenome supports the shorter of the duplex reads, so we will
          // clip this base out
          include_base = false;
        } else {
          // is this a substitution operation?
          char secondary = SUBSTITUTIONS[((uint16_t)yc_op) << 1];
          if (secondary != '\0') {
            if (secondary == base) {
              // we actually grabbed the primary base, so switch to the other
              // one
              secondary = SUBSTITUTIONS[(((uint16_t)yc_op) << 1) | 1];
            }
            if (aligned_chars[read_idx] == base) {
              // the pangenome supports the current base call, bump up its base
              // quality
              bq = ADJUSTED_BQ;
            } else if (aligned_chars[read_idx] == secondary) {
              // the pangenome supports the discordant base call, switch to it
              // and bump up the quality
              base = secondary;
              bq = ADJUSTED_BQ;
            }
          }
        }
      }

      if (include_base) {
        // preserve the info corresponding to this base in the output
        mod_cigar.push_back(cigar_op);
        mod_seq.push_back(base);
        mod_qual.push_back(bq);
        mod_yc.push_back(yc_op);
        have_consumed_ref = have_consumed_ref || op_consumes_ref[(uint8_t)cigar_op];
      } else if (op_consumes_ref[(uint8_t)expanded_cigar[i]]) {
        // we are deleting a base aligned with an operation that consumed
        // reference sequence, so we need to replace it with a deletion to
        // maintain the integrity of the remaining alignment
        mod_cigar.push_back('D');
        if (!have_consumed_ref) {
          // this will require adjusting the position of the alignment
          ++first_removed_run;
        }
      }
      ++read_idx;
    }

    // check if we've moved past this BED interval
    if (bed_cursor < exclusion_bed.size() && op_consumes_ref[(uint8_t)cigar_op]) {
      ++ref_cursor;
      if (ref_cursor >= std::get<2>(exclusion_bed[bed_cursor])) {
        ++bed_cursor;
      }
    }
  }

  // figure out the interval of the modified cigar that is inside any softclips
  size_t nonclip_begin = 0;
  size_t nonclip_end = mod_cigar.size();
  while (nonclip_begin < mod_cigar.size() && (mod_cigar[nonclip_begin] == 'S' || mod_cigar[nonclip_begin] == 'H')) {
    ++nonclip_begin;
  }
  while (nonclip_end != 0 && (mod_cigar[nonclip_end - 1] == 'S' || mod_cigar[nonclip_end - 1] == 'H')) {
    --nonclip_end;
  }

  if (nonclip_begin < nonclip_end && (mod_cigar[nonclip_begin] == 'D' || mod_cigar[nonclip_begin] == 'N' ||
                                      mod_cigar[nonclip_end - 1] == 'D' || mod_cigar[nonclip_end - 1] == 'N')) {
    // non-read-consuming operations have been left "hanging" at the end of the
    // alignment, which some tools see as invalid, so we have to remove them

    size_t inner_begin = nonclip_begin;
    size_t inner_end = nonclip_end;
    while (inner_begin < inner_end && (mod_cigar[inner_begin] == 'D' || mod_cigar[inner_begin] == 'N')) {
      ++inner_begin;
      ++first_removed_run;
    }
    while (inner_end > inner_begin && (mod_cigar[inner_end - 1] == 'D' || mod_cigar[inner_end - 1] == 'N')) {
      --inner_end;
    }
    mod_cigar.erase(inner_end, nonclip_end - inner_end);
    mod_cigar.erase(nonclip_begin, inner_begin - nonclip_begin);
  }

  return reencode_read(read, mod_seq, mod_qual, mod_yc, mod_cigar, first_removed_run, parsed_yc);
}

int VgDuplexClipperMain(int argc, char** argv) {
  size_t num_threads = 1;

  std::string bed_file;

  while (true) {
    static struct option options[] = {{"threads", required_argument, NULL, 't'},
                                      {"bed", required_argument, NULL, 'b'},
                                      {"adjusted-qual", required_argument, NULL, 'q'},
                                      {"help", no_argument, NULL, 'h'},
                                      {NULL, 0, NULL, 0}};
    int o = getopt_long(argc, argv, "t:b:q:h", options, NULL);

    if (o == -1) {
      // end of options
      break;
    }
    switch (o) {
      case 't':
        num_threads = std::stoul(std::string(optarg));
        break;
      case 'b':
        bed_file = optarg;
        break;
      case 'q':
        ADJUSTED_BQ = std::stoul(std::string(optarg));
        break;
      case 'h':
        print_help();
        return 0;
      default:
        print_help();
        return 1;
    }
  }

  if (argc - optind > 1) {
    std::cerr << "error: unused positional arguments\n";
    print_help();
    return 1;
  }

  // open BAM from positional argument or from
  samFile* bam_in = NULL;
  const char* bam_file = NULL;
  if (argc - optind == 0) {
    bam_file = "-";
  } else {
    bam_file = argv[optind];
  }
  bam_in = sam_open(bam_file, "rb");
  if (!bam_in) {
    std::cerr << "error: could not open BAM file: " << bam_file << '\n';
    return 1;
  }

  bam_hdr_t* header = sam_hdr_read(bam_in);
  if (!header) {
    std::cerr << "error: could not read header in BAM file: " << bam_file << '\n';
    sam_close(bam_in);
    return 1;
  }

  samFile* bam_out = sam_open("-", "wb");
  if (!bam_out) {
    std::cerr << "error: could not open stdout as a BAM file\n";
    bam_hdr_destroy(header);
    sam_close(bam_in);
    return 1;
  }

  // read in the BED of excluded regions (if any)
  std::vector<std::tuple<int, size_t, size_t>> exclusion_bed;
  if (!bed_file.empty()) {
    std::ifstream bed_in(bed_file);
    if (!bed_in) {
      std::cerr << "error: could not open BED file " << bed_file << '\n';
      bam_hdr_destroy(header);
      sam_close(bam_in);
      sam_close(bam_out);
      return 1;
    }
    exclusion_bed = std::move(parse_bed(bed_in, header));
  }

  if (sam_hdr_write(bam_out, header) < 0) {
    std::cerr << "error: could not write BAM header to stdout\n";
    bam_hdr_destroy(header);
    sam_close(bam_in);
    sam_close(bam_out);
    return 1;
  }

  // the size of one batch of reads sent to a worker thread
  static const size_t batch_size = (1 << 13);  // ~8000
  // the number of batches we will aim to keep in the queue
  size_t num_queue_batches = 16 * num_threads;
  // the largest number of batches we will possibly keep in the queue if we
  // increase the size
  size_t max_num_queue_batches = 64 * num_queue_batches;

  // have we finished reading in the BAM?
  std::atomic_bool finished_reading(false);
  // queue of batches
  std::list<std::vector<bam1_t*>> batches;
  // lock for intereacting with the queue
  std::mutex queue_lock;
  // lock for writing to stdout
  std::mutex write_lock;

  // task in which we do the clipping and write the output, either in a loop or
  // for a single iteration
  auto follower_task = [&](bool in_loop) {
    while (true) {
      // try to get a batch
      bool was_empty = false;
      std::vector<bam1_t*> batch;
      queue_lock.lock();
      if (batches.empty()) {
        queue_lock.unlock();
        was_empty = true;
      } else {
        batch = std::move(batches.front());
        batches.pop_front();
        queue_lock.unlock();
      }

      if (!was_empty) {
        // we got a batch, do the clipping
        std::vector<bam1_t*> clipped;
        clipped.reserve(batch.size());
        for (auto read : batch) {
          clipped.push_back(clip_discordant(read, exclusion_bed));
        }

        // write the output
        write_lock.lock();
        for (auto mod_read : clipped) {
          if (sam_write1(bam_out, header, mod_read) < 0) {
            std::cerr << "error: failed to write BAM records to stdout\n";
            exit(1);
          }
        }
        write_lock.unlock();

        // clean up
        for (size_t i = 0; i < batch.size(); ++i) {
          bam_destroy1(batch[i]);
          bam_destroy1(clipped[i]);
        }

        if (!in_loop) {
          // we've finished with the one iteration we want to do
          break;
        }
      } else if (!in_loop || finished_reading.load()) {
        // the queue is empty and we don't need to wait for any more
        break;
      } else {
        // the queue is empty so lets sleep to wait for it to fill up again (1
        // millisec)
        usleep(1000);
      }
    }
  };
  // task in which we do file reading, and when the queue is full also clipping
  auto leader_task = [&]() {
    bool more_left = true;  // redundant on finished_reading, but avoids
                            // contention on it with the followers
    while (more_left) {
      // TODO: is it possible to do this with an atomic int instead and still
      // have it track to the true size of the queue?
      queue_lock.lock();
      size_t queue_size = batches.size();
      queue_lock.unlock();

      if (queue_size == num_queue_batches && num_queue_batches < max_num_queue_batches) {
        // the queue is filled up, we can do some processing (which should
        // reduce the queue's size to below the max limit)

        follower_task(false);

        queue_lock.lock();
        queue_size = batches.size();
        queue_lock.unlock();
        if (queue_size < num_queue_batches / 2) {
          // too much of the queue emptied while we were doing work, increase
          // the size to avoid this happening again
          num_queue_batches *= 2;
        }
      }

      // read in a batch
      std::vector<bam1_t*> batch;
      batch.reserve(batch_size);
      for (size_t i = 0; i < batch_size; ++i) {
        bam1_t* read = bam_init1();
        int code = sam_read1(bam_in, header, read);
        if (code >= 0) {
          batch.push_back(read);
        } else {
          more_left = false;
          finished_reading.store(true);
        }
      }

      // add it to the queue
      queue_lock.lock();
      batches.emplace_back(std::move(batch));
      queue_lock.unlock();
    }

    // enter the follower task in a loop to help clear out the remaining queue
    follower_task(true);
  };

  // one leader thread
  std::vector<std::thread> workers;
  workers.emplace_back(leader_task);
  // the rest as followers
  while (workers.size() < num_threads) {
    workers.emplace_back(follower_task, true);
  }
  // barrier sync
  for (auto& worker : workers) {
    worker.join();
  }

  // cleanup
  sam_close(bam_out);
  sam_close(bam_in);
  bam_hdr_destroy(header);

  return 0;
}

}  // namespace xoos::pangenome_consensus_caller
