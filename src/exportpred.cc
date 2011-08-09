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
#include "predict_pexel.hh"

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



#if VERSION == 1

static int a_spacer_raw_distrib[] = {
  1, 2, 3, 6, 8, 8, 9, 10, 10, 11, 11, 12, 13, 13, 13, 13, 14, 14,
  14, 14, 15, 15, 15, 16, 16, 17, 17, 17, 17, 17, 17, 18, 18, 18, 18,
  18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
  19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21, 21, 21, 21, 21, 21, 21,
  21, 22, 22, 22, 22, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 24, 24,
  24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25, 25, 25, 26, 26,
  26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 27, 27,
  27, 28, 29, 29, 29, 29, 29, 30, 31, 31, 31, 31, 31, 32, 32, 33, 35,
  37, 41, 43, 48, 53, 61, 89, 110, 129, 166, 182, 202, 209, 324, 495
};

#endif

#if VERSION == 2

static int a_spacer_raw_distrib[] = {
  9, 12, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 16,
  16, 16, 17, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18, 18,
  18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19,
  19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20,
  20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 20, 21,
  21, 21, 21, 21, 21, 22, 22, 22, 23, 23, 23, 23, 23, 23, 23, 23, 23,
  24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 24, 25, 25, 25, 25, 25, 25,
  25, 25, 25, 25, 25, 25, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26, 26,
  26, 26, 26, 26, 26, 26, 26, 26, 27, 27, 27, 29, 31, 31, 31, 31, 33,
  33
};

#endif

GHMM::Model::Ptr makePEXELmodel() {
  GHMM::UTIL::Alphabet::Ptr alphabet = new GHMM::UTIL::Alphabet();
  alphabet->addCharTokenRange('A','Z');
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::EMISSION::Base::Ptr background, met;

  std::vector<MATH::DPDF::Ptr> RLE(7, MATH::DPDF::Ptr());

  GHMM::LENGTH::Discrete::Ptr a_leader_length, a_spacer_length;

  GHMM::LENGTH::Geometric::Ptr a_tail_length, c_tail_length;

#if VERSION == 1
  GHMM::LENGTH::Uniform::Ptr a_hydrophobic_length;
#endif

#if VERSION == 2
  GHMM::LENGTH::Discrete::Ptr a_hydrophobic_length;
#endif

  a_spacer_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                            a_spacer_raw_distrib,
                                                            ENDOF(a_spacer_raw_distrib),
                                                            1,
                                                            60));

  background = new GHMM::EMISSION::Stateless(ep->parse("\n\
             A:  78883 C:  71359 D: 260979 E: 288230\n\
             F: 175488 G: 114068 H:  97688 I: 373389\n\
             K: 473828 L: 304967 M:  88773 N: 581084\n\
             P:  80295 Q: 111860 R: 106760 S: 256676\n\
             T: 164816 V: 154088 W:  19966 Y: 230362"));

  met = new GHMM::EMISSION::Stateless(ep->parse("M: 1"));

#if VERSION == 1
  RLE[0] = ep->parse("A: 1 C: 6 D: 2 E: 2 F: 11 G: 7 H: 3 I: 15 K: 23 L: 6 M: 1 N: 15 P: 1 Q: 3 R: 6 S: 38 T: 5 V: 3 W: 1 Y: 6");
  RLE[1] = ep->parse("K: 3 R: 152");
  // RLE[2] = ep->parse("I: 28 K: 11 L: 15 N: 32 S: 26 T: 9");
  RLE[2] = ep->parse("A: 2 C: 6 E: 1 F: 4 H: 3 I: 28 K: 11 L: 15 M: 2 N: 32 Q: 3 R: 3 S: 26 T: 9 V: 6 W: 1 Y: 3");
  // RLE[3] = ep->parse("L: 151");
  RLE[3] = ep->parse("F: 1 I: 2 L: 151 N: 1");
  // RLE[4] = ep->parse("A: 35 S: 48 T: 18 Y: 15 C: 9");
  RLE[4] = ep->parse("A: 35 C: 9 E: 2 F: 2 G: 5 H: 1 I: 2 K: 1 L: 3 N: 7 S: 48 T: 18 V: 7 Y: 15");
  // RLE[5] = ep->parse("D: 11 E: 109 Q: 21");
  RLE[5] = ep->parse("C: 2 D: 11 E: 109 G: 2 H: 1 K: 1 Q: 21 S: 4 T: 3 Y: 1");
  RLE[6] = ep->parse("A: 1 C: 4 D: 1 E: 6 F: 5 G: 6 H: 7 I: 5 K: 10 L: 14 M: 3 N: 15 P: 9 Q: 3 R: 5 S: 10 T: 17 V: 18 Y: 16");
#endif

#if VERSION == 2
  RLE[0] = ep->parse("M: 1 P: 1 W: 1 A: 2 D: 2 E: 2 H: 3 Q: 3 V: 4 T: 5 C: 6 Y: 6 G: 7 L: 7 R: 7 F: 13 I: 15 N: 15 K: 26 S: 40");
  RLE[1] = ep->parse("K: 7 R: 159");
  RLE[2] = ep->parse("W: 1 A: 2 E: 2 M: 2 H: 3 R: 3 Q: 4 Y: 4 F: 5 V: 6 C: 7 T: 9 K: 12 L: 16 S: 28 I: 29 N: 33");
  RLE[3] = ep->parse("N: 1 F: 2 I: 2 L: 161");
  RLE[4] = ep->parse("E: 1 H: 1 F: 2 I: 2 K: 3 L: 3 G: 4 N: 8 V: 8 C: 10 Y: 17 T: 18 A: 37 S: 52");
  RLE[5] = ep->parse("H: 1 K: 1 Y: 1 C: 2 G: 3 T: 3 S: 5 D: 15 Q: 21 E: 114");
  RLE[6] = ep->parse("D: 1 A: 3 M: 3 Q: 3 C: 4 F: 5 R: 5 G: 6 I: 6 E: 7 H: 8 P: 9 K: 10 S: 11 L: 15 N: 17 T: 17 Y: 17 V: 19");
#endif

  a_tail_length = new GHMM::LENGTH::Geometric(364);

  c_tail_length = new GHMM::LENGTH::Geometric(755);

  GHMM::StateBase::Ptr a_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr a_spacer = GHMM::UTIL::makeState(a_spacer_length, background);
  GHMM::StateBase::Ptr a_RLE = GHMM::UTIL::makeMotifState(RLE);
  GHMM::StateBase::Ptr a_tail = GHMM::UTIL::makeState(a_tail_length, background);

  GHMM::StateBase::Ptr c_met = GHMM::UTIL::makeState(NULL, met);
  GHMM::StateBase::Ptr c_tail = GHMM::UTIL::makeState(c_tail_length, background);

  GHMM::ModelBuilder mb;

  std::pair<std::string, std::string> makeSignalPModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);
  std::pair<std::string, std::string> makeSSModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet);

#ifdef SIGNALP_MODEL
  std::pair<std::string, std::string> ss_states = makeSignalPModel(mb, alphabet);
#else
  std::pair<std::string, std::string> ss_states = makeSSModel(mb, alphabet);
#endif

#ifdef RLE_PATTERN
  mb.addState("a-met",         a_met);
  mb.addState("a-spacer",      a_spacer);
  mb.addState("a-RLE",         a_RLE);
  mb.addState("a-tail",        a_tail);
#endif

  mb.addState("c-met",         c_met);
  mb.addState("c-tail",        c_tail);

#ifdef RLE_PATTERN
  mb.addStateTransition(GHMM::Model::BEGIN, "a-met", 400);
#endif

  mb.addStateTransition(GHMM::Model::BEGIN, "c-met", 5009);

#ifdef RLE_PATTERN
  mb.addStateTransition("a-met",            ss_states.first,  1);
  mb.addStateTransition(ss_states.second,   "a-spacer",       1);
  mb.addStateTransition("a-spacer",         "a-RLE",          1);
  mb.addStateTransition("a-RLE",            "a-tail",         1);
  mb.addStateTransition("a-tail",           GHMM::Model::END, 1);

  mb.addState("d-tail", c_tail);
  mb.addStateTransition(ss_states.second,   "d-tail",         0.01);
  mb.addStateTransition("d-tail",           GHMM::Model::END, 1);
#endif

  mb.addStateTransition("c-met",            "c-tail",         1);
  mb.addStateTransition("c-tail",           GHMM::Model::END, 1);

  return mb.make();
}

static const struct option options[] = {
  { "input",             required_argument,          0,            'i' },
  { "output",            required_argument,          0,            'o' },
  { "RLE-threshold",     required_argument,          0,            'R' },
  { "no-RLD",            no_argument,                0,            'r' },
  { 0,                   0,                          0,            0   }
};

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
                                       (default: 4.3)\n\
--no-RLE                -r             turn off RLE prediction\n\
\n\
";
}

int main(int argc, char **argv) {
  double RLE_threshold = 4.3;
  bool do_RLE = true;

  std::list<std::pair<std::string, std::string> > seq_list;
  std::string output = "-";

  int ch;

  while ((ch = getopt_long(argc, argv, "i:o:R:K:hkr", options, NULL)) != -1) {
    switch (ch) {
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
    case 'o': {
      output = optarg;
      break;
    }
    case 'R': {
      RLE_threshold = strtod(optarg, NULL);
      break;
    }
    case 'r': {
      do_RLE = false;
      break;
    }
    case 'h':
    case '?': {
      usage(argv[0]);
      exit(0);
    }
    }
  }

  GHMM::Model::Ptr model = makePEXELmodel();

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
      out << name << "\t"
          << "RLE" << "\t"
          << alpha_rle - alpha_bkg << "\t"
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
