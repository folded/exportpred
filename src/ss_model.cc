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

static int leader_raw_distrib[] = {
  1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 4, 4, 4, 4, 4, 5, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 8, 8, 8, 8, 9, 9, 10, 10, 11, 12, 12, 12, 12, 13, 13, 13, 14, 15, 15, 15, 16, 16, 17, 19, 19, 20, 20, 23, 24, 25, 25, 26, 26, 28, 28, 28, 28, 29, 29, 29, 32, 35, 37, 37, 37, 39, 39, 40, 42, 43, 43, 46, 47, 48, 48, 49, 49, 49, 49, 49, 49, 49, 50, 50, 50, 50, 50, 51, 51, 51, 51, 51, 51, 51, 51, 51, 51, 52, 52, 52, 52, 52, 52, 53, 53, 54, 58
};

static int hydrophobic_raw_distrib[] = {
  10, 11, 11, 11, 11, 11, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 13, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 14, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 15, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 16, 17, 17, 17, 18, 18, 18, 18, 18, 18, 18, 18, 18, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 19, 20, 20, 20, 20, 20, 20, 22, 22, 23, 24, 24, 24, 24
};

std::pair<std::string, std::string> makeSSModel(GHMM::ModelBuilder &mb, GHMM::UTIL::Alphabet::Ptr &alphabet) {
  GHMM::UTIL::EmissionDistributionParser::Ptr ep = new GHMM::UTIL::EmissionDistributionParser(alphabet);

  GHMM::LENGTH::Discrete::Ptr leader_length;
  GHMM::EMISSION::Base::Ptr ss_leader, hydrophobic_distrib;

  GHMM::LENGTH::Discrete::Ptr hydrophobic_length;

  leader_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(2.0),
                                                            leader_raw_distrib,
                                                            ENDOF(leader_raw_distrib),
                                                            1,
                                                            80));

  hydrophobic_length = new GHMM::LENGTH::Discrete(MATH::smooth(MATH::GaussianKernel(1.0),
                                                                 hydrophobic_raw_distrib,
                                                                 ENDOF(hydrophobic_raw_distrib),
                                                                 10,
                                                                 25));

  ss_leader = new GHMM::EMISSION::Stateless(ep->parse("\n\
             W:  20 P:  44 A:  52 Q:  73 H:  74 G:  90\n\
             D: 103 V: 109 M: 111 C: 116 T: 127 E: 156\n\
             R: 191 L: 239 F: 248 Y: 274 I: 285 S: 366\n\
             N: 516 K: 639"));

  hydrophobic_distrib = new GHMM::EMISSION::Stateless(ep->parse("\n\
             D:   4 E:   6 Q:   6 R:   6 H:  11 K:  24\n\
             P:  24 W:  24 M:  27 A:  37 N:  40 G:  70\n\
             T:  87 C: 106 S: 113 Y: 132 V: 199 F: 396\n\
             I: 508 L: 559"));

  GHMM::StateBase::Ptr leader = GHMM::UTIL::makeState(leader_length, ss_leader);
  GHMM::StateBase::Ptr hydrophobic = GHMM::UTIL::makeState(hydrophobic_length, hydrophobic_distrib);

  mb.addState("a-leader",      leader);
  mb.addState("a-hydrophobic", hydrophobic);

  mb.addStateTransition("a-leader", "a-hydrophobic",  1);

  return std::make_pair(std::string("a-leader"), std::string("a-hydrophobic"));

}

