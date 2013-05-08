// Copyright (c) 2005 The Walter and Eliza Hall Institute
// 
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
// 
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
// 
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
// IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR
// ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
// CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

#include "exportpred.hh"

#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <fstream>
#include <iterator>

#include <GHMM/string_funcs.hh>

#include <getopt.h>
#include <string.h>

#define BLOCK_SIZE 1024

template<typename OutputIterator>
static inline void addSeq(OutputIterator &out,
                          const std::string &name,
                          const std::vector<char *> &seq,
                          int seq_len) {
  std::string s;
  s.reserve(seq_len);

  for (int i = 0; i < seq_len; i += BLOCK_SIZE) {
    s.append(seq[i / BLOCK_SIZE], std::min(BLOCK_SIZE, seq_len - i));
  }
  assert(s.size() ==(unsigned int)seq_len);
  *out++ = std::make_pair(name, s);
}

static inline void addChar(std::vector<char *> &seq,
                           int &seq_len,
                           int &seq_size,
                           char c) {
  if (!isspace(c)) {
    while (seq_len >= seq_size) {
      seq.push_back(new char[BLOCK_SIZE]);
      seq_size += BLOCK_SIZE;
    }
    seq[seq_len / BLOCK_SIZE][seq_len % BLOCK_SIZE] = c;
    seq_len++;
  }
}

template<typename OutputIterator>
void readFasta(std::istream &in,
               OutputIterator out,
               std::string (*name_xform)(const std::string &)) {
  std::string name = "";
  int c;
  std::vector<char *> seq;

  std::string l;
  int seq_len = 0;
  int seq_size = 0;
  int state = 0;

  while ((c = in.get()) != EOF) {
    switch (state) {
    case 0: {
      // looking for name
      if (c == '>') {
        state = 1;
        name = "";
      }
      break;
    }
    case 1: {
      // name
      if (c == '\n') {
        state = 2;
        if (name_xform) name = name_xform(name);
      } else {
        name.push_back(c);
      }
      break;
    }
    case 2: {
      // sequence, first char of column
      if (c == '>') {
        if (seq_len) addSeq(out, name, seq, seq_len);
        seq_len = 0;
        state = 1;
        name = "";
      } else if (c != '\n') {
        addChar(seq, seq_len, seq_size, c);
        state = 3;
      }
      break;
    }
    case 3: {
      // sequence, not first char of column
      if (c == '\n') {
        state = 2;
      } else {
        addChar(seq, seq_len, seq_size, c);
      }
      break;
    }
    }
  }

  if (seq_len) addSeq(out, name, seq, seq_len);

  for (std::vector<char *>::iterator i = seq.begin(); i != seq.end(); ++i) {
    delete [] *i;
  }
}



static int a_spacer_raw_distrib[] = {
  12, 13, 13, 14, 14, 14, 14, 15, 15, 15, 15, 16, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 27, 27, 27, 28, 30, 31, 31, 31, 33
};



GHMM::Model::Ptr makePEXELmodel(bool signalp_model, const std::string exclude_5) {
  GHMM::UTIL::Alphabet::Ptr alphabet = new GHMM::UTIL::Alphabet();
  alphabet->addCharTokenRange('A','Z');
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::EMISSION::Base::Ptr background, met;

  GHMM::LENGTH::Discrete::Ptr a_leader_length, a_spacer_length;

  GHMM::LENGTH::Geometric::Ptr a_tail_length, c_tail_length;

  GHMM::LENGTH::Discrete::Ptr a_hydrophobic_length;

  a_spacer_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                            a_spacer_raw_distrib,
                                                            ENDOF(a_spacer_raw_distrib),
                                                            1,
                                                            60));

  MATH::DPDF::Ptr background_distrib = ep->parse("\n\
             A:  78883 C:  71359 D: 260979 E: 288230\n\
             F: 175488 G: 114068 H:  97688 I: 373389\n\
             K: 473828 L: 304967 M:  88773 N: 581084\n\
             P:  80295 Q: 111860 R: 106760 S: 256676\n\
             T: 164816 V: 154088 W:  19966 Y: 230362");
  background = new GHMM::EMISSION::Stateless(background_distrib);

  met = new GHMM::EMISSION::Stateless(ep->parse("M: 1"));

  std::vector<MATH::DPDF::Ptr> RLE;
  std::ostringstream pos5;
#define ADD5(ch, count) do { if (exclude_5.find(ch) == std::string::npos) { pos5 << ch << ": " << count << " "; } } while(0)
  ADD5('C', 1);
  ADD5('D', 14);
  ADD5('E', 100);
  ADD5('F', 1);
  ADD5('G', 2);
  ADD5('H', 1);
  ADD5('I', 1);
  ADD5('K', 2);
  ADD5('L', 1);
  ADD5('M', 1);
  ADD5('N', 2);
  ADD5('P', 1);
  ADD5('Q', 22);
  ADD5('R', 1);
  ADD5('S', 6);
  ADD5('T', 4);
  ADD5('V', 1);
  ADD5('W', 1);
  ADD5('Y', 2);
  std::string pos5str = pos5.str();

  ep->parseMotif(std::back_inserter(RLE),
      "A: 6 C: 10 D: 7 E: 7 F: 16 G: 12 H: 8 I: 18 K: 28 L: 11 M: 6 N: 20 P: 6 Q: 8 R: 9 S: 38 T: 10 V: 8 W: 6 Y: 11",
      "R: 145",
      "A: 3 C: 7 D: 1 E: 2 F: 5 G: 1 H: 4 I: 26 K: 11 L: 15 M: 4 N: 26 P: 1 Q: 4 R: 4 S: 27 T: 10 V: 7 W: 2 Y: 5",
      "L: 1",
      "A: 33 C: 9 D: 1 E: 2 F: 3 G: 5 H: 2 I: 3 K: 2 L: 4 M: 1 N: 7 P: 1 Q: 1 R: 1 S: 48 T: 19 V: 8 W: 1 Y: 14",
      // "C: 1 D: 14 E: 100 F: 1 G: 2 H: 1 I: 1 K: 2 L: 1 M: 1 N: 2 P: 1 Q: 22 R: 1 S: 6 T: 4 V: 1 W: 1 Y: 2",
      pos5str.c_str(),
      "A: 7 C: 9 D: 6 E: 7 F: 10 G: 10 H: 13 I: 10 K: 15 L: 17 M: 7 N: 21 P: 13 Q: 8 R: 10 S: 16 T: 21 V: 21 W: 5 Y: 19",
      NULL);

  std::vector<MATH::DPDF::Ptr> RLE_late;
  ep->parseMotif(std::back_inserter(RLE_late),
      "A: 6 C: 10 D: 7 E: 7 F: 16 G: 12 H: 8 I: 18 K: 28 L: 11 M: 6 N: 20 P: 6 Q: 8 R: 9 S: 38 T: 10 V: 8 W: 6 Y: 11",
      "R: 145",
      "A: 3 C: 7 D: 1 E: 2 F: 5 G: 1 H: 4 I: 26 K: 11 L: 15 M: 4 N: 26 P: 1 Q: 4 R: 4 S: 27 T: 10 V: 7 W: 2 Y: 5",
      "L: 1",
      "A: 33 C: 9 D: 1 E: 2 F: 3 G: 5 H: 2 I: 3 K: 2 L: 4 M: 1 N: 7 P: 1 Q: 1 R: 1 S: 48 T: 19 V: 8 W: 1 Y: 14",
      "A: 33 C: 9 D: 1 E: 2 F: 3 G: 5 H: 2 I: 3 K: 2 L: 4 M: 1 N: 7 P: 1 Q: 1 R: 1 S: 48 T: 19 V: 8 W: 1 Y: 14",
      // "C: 1 D: 14 E: 100 F: 1 G: 2 H: 1 I: 1 K: 2 L: 1 M: 1 N: 2 P: 1 Q: 22 R: 1 S: 6 T: 4 V: 1 W: 1 Y: 2",
      pos5.str().c_str(),
      "A: 7 C: 9 D: 6 E: 7 F: 10 G: 10 H: 13 I: 10 K: 15 L: 17 M: 7 N: 21 P: 13 Q: 8 R: 10 S: 16 T: 21 V: 21 W: 5 Y: 19",
      NULL);

  a_tail_length = new GHMM::LENGTH::Geometric(364);

  c_tail_length = new GHMM::LENGTH::Geometric(755);

  GHMM::StateBase::Ptr a_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr a_spacer = GHMM::UTIL::makeState(a_spacer_length, background);
  GHMM::StateBase::Ptr a_RLE = GHMM::UTIL::makeMotifState(RLE);
  GHMM::StateBase::Ptr a_RLE_late = GHMM::UTIL::makeMotifState(RLE_late);
  GHMM::StateBase::Ptr a_tail = GHMM::UTIL::makeState(a_tail_length, background);

  GHMM::StateBase::Ptr c_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr c_tail = GHMM::UTIL::makeState(c_tail_length, background);

  GHMM::ModelBuilder mb;

  std::pair<std::string, std::string> makeSignalPModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);
  std::pair<std::string, std::string> makeSSModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);

  std::pair<std::string, std::string> ss_states;
  if (signalp_model) {
    ss_states = makeSignalPModel(mb, alphabet);
  } else {
    ss_states = makeSSModel(mb, alphabet);
  }

  mb.addState("a-met",         a_met);
  mb.addState("a-spacer",      a_spacer);
  mb.addState("a-RLE",         a_RLE);
  mb.addState("a-RLE-late",    a_RLE_late);
  mb.addState("a-tail",        a_tail);

  mb.addState("c-met",         c_met);
  mb.addState("c-tail",        c_tail);

  mb.addStateTransition(GHMM::Model::BEGIN, "a-met", 400);

  mb.addStateTransition(GHMM::Model::BEGIN, "c-met", 5009);

  mb.addStateTransition("a-met",            ss_states.first,  1);
  mb.addStateTransition(ss_states.second,   "a-spacer",       1);
  mb.addStateTransition("a-spacer",         "a-RLE",          100);
  mb.addStateTransition("a-spacer",         "a-RLE-late",     5);
  mb.addStateTransition("a-RLE",            "a-tail",         1);
  mb.addStateTransition("a-RLE-late",       "a-tail",         1);
  mb.addStateTransition("a-tail",           GHMM::Model::END, 1);

  mb.addState("d-tail", c_tail);
  mb.addStateTransition(ss_states.second,   "d-tail",         0.01);
  mb.addStateTransition("d-tail",           GHMM::Model::END, 1);

  mb.addStateTransition("c-met",            "c-tail",         1);
  mb.addStateTransition("c-tail",           GHMM::Model::END, 1);

  return mb.make();
}

static const struct option options[] = {
//{ "signalp",           no_argument,                0,            's' },
  { "input",             required_argument,          0,            'i' },
  { "exclude-5",         required_argument,          0,            '5' },
  { "output",            required_argument,          0,            'o' },
  { "RLE-threshold",     required_argument,          0,            'R' },
  { 0,                   0,                          0,            0   }
};

bool parseContainsState(const GHMM::Model::Ptr &model, int state, GHMM::Traceback::Ptr tbp) {
  while (tbp != NULL) {
    if (tbp->state == state) return true;
    tbp = tbp->prev;
  }
  return false;
}

std::string genParse(const std::string &sequence, const GHMM::Model::Ptr &model, GHMM::Traceback::Ptr tbp) {
  std::ostringstream out;
  std::vector<std::string> parse;
  std::string::size_type pos = sequence.size();
  while (tbp != NULL) {
    std::ostringstream o;
    pos -= tbp->length;
    o << "[" << model->stateName(tbp->state) << ":" << sequence.substr(pos, tbp->length) << "]";
    parse.push_back(o.str());
    tbp = tbp->prev;
  }
  std::copy(parse.rbegin(), parse.rend(), std::ostream_iterator<std::string>(out, ""));
  return out.str();
}

void usage(const char *progname) {
  std::cout << "Usage: " << progname << " [arguments]" << std::endl;

  std::cout << "\
\n\
--input=file            -i file        read sequences from file (-:stdin)\n\
--output=file           -o file        write results to file (-:stdout)\n\
--RLE-threshold=float   -R float       RLE threshold for positive prediction\n\
                                       (default: 0.0)\n\
--exclude-5             -5 str         exclude given residues in position 5\n\
\n\
";
//--signalp               -s             use SignalP-based HMM for SS prediction
//                                       (default: off)
}

int main(int argc, char **argv) {
  double RLE_threshold = 0.0;

  std::list<std::pair<std::string, std::string> > seq_list;
  std::string output = "-";

  int ch;
  bool signalp_model = false;
  std::string exclude_5;

  while ((ch = getopt_long(argc, argv, "i:o:R:K:hkr" /* "s" */, options, NULL)) != -1) {
    switch (ch) {
//     case 's': {
//       signalp_model = true;
//       break;
//     }
    case 'i': {
      if (!strcmp(optarg, "-")) {
        std::ifstream in(optarg);
        int old_sz = seq_list.size();
        readFasta(std::cin, std::inserter(seq_list, seq_list.begin()), NULL);
        std::cerr << seq_list.size() - old_sz << " sequences read from stdin" << std::endl;
      } else {
        std::ifstream in(optarg);
        int old_sz = seq_list.size();
        readFasta(in, std::inserter(seq_list, seq_list.begin()), NULL);
        std::cerr << seq_list.size() - old_sz << " sequences read from " << optarg << std::endl;
      }
      break;
    }
    case '5': {
      exclude_5 += optarg;
      break;
    }
    case 'o': {
      output = optarg;
      break;
    }
    case 'R': {
      RLE_threshold = strtod(optarg, NULL);
      break;
    }
    case 'h':
    case '?': {
      usage(argv[0]);
      exit(0);
    }
    }
  }

  GHMM::Model::Ptr model = makePEXELmodel(signalp_model, exclude_5);

  std::list<std::pair<std::string, std::string> >::iterator i, e;

  std::vector<std::pair<double, std::string> > rle_out;

  for (i = seq_list.begin(), e = seq_list.end(); i != e; ++i) {
    std::string &name((*i).first);
    std::string &sequence((*i).second);
    int *seq_raw = new int[sequence.size()];

    for (int j = 0; j < (int)sequence.size(); j++) {
      if (isalpha(sequence[j])) {
        seq_raw[j] = toupper(sequence[j]) - 'A';
      } else {
        seq_raw[j] = 'X' - 'A';
      }
    }

    GHMM::Parse::Ptr parse = new GHMM::Parse();
    parse->parse(model, seq_raw, seq_raw + sequence.size());
    double alpha_rle, alpha_bkg;
    alpha_rle = parse->alpha(model->stateNumber("a-tail"), 0);
    alpha_bkg = parse->alpha(model->stateNumber("c-tail"), 0);
#if 0
    std::cerr << name
              << " alpha_rle:" << alpha_rle
              << " alpha_bkg:" << alpha_bkg
              << " alpha_ssonly=" << parse->alpha(model->stateNumber("d-tail"), 0) << std::endl;
#endif

    if (alpha_rle - alpha_bkg > RLE_threshold) {
      std::ostringstream out;
      const double mu = 15.55600978;
      const double lambda = 0.08604353;

      double x =  alpha_rle - alpha_bkg;
      double px = exp(-exp(-lambda * -x * mu));
      out << name << "\t"
          << (parseContainsState(model,
                                 model->stateNumber("a-RLE-late"),
                                 parse->psi(model->stateNumber("a-tail"), 0)) ?
              "RLE-late" : "RLE") << "\t"
          << alpha_rle - alpha_bkg << "\t"
          << px << "\t"
          << genParse(sequence, model, parse->psi(model->stateNumber("a-tail"), 0));
      rle_out.push_back(std::make_pair(alpha_rle - alpha_bkg, out.str()));
    }
  }

  std::sort(rle_out.begin(), rle_out.end());

  if (output == "-") {
    for (std::vector<std::pair<double, std::string> >::reverse_iterator i = rle_out.rbegin(); i != rle_out.rend(); ++i) {
      std::cout << (*i).second << std::endl;
    }
  } else {
    std::ofstream out;
    out.open(output.c_str());

    for (std::vector<std::pair<double, std::string> >::reverse_iterator i = rle_out.rbegin(); i != rle_out.rend(); ++i) {
      out << (*i).second << std::endl;
    }
  }
}
