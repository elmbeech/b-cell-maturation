////////
// date: 2023-03-20
// license: BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)
// author: bue
//
// description:
//   qunatify hamming distance (affinity) from to protein sequences of variable size.
////////

#include <cstdio>
#include <vector>

#include "./custom.h"
#include "./hamming.h"

double alignment( Vector_Variable antigenSequence, Vector_Variable antibodySequence ) {

    // strip and print input
    std::vector<double> antigen_copy {};
    std::vector<double> antibody_copy {};

    for (double element : antigenSequence.value) {
        if (element != PAD) {
            antigen_copy.push_back(element);
        }
    }
    for (double element : antibodySequence.value) {
        if (element != PAD) {
            antibody_copy.push_back(element);
        }
    }

    printf("***antigen_sequence: ");
    for (double element : antigen_copy) {
        printf("{%g}", element);
    }
    printf("***.\n");
    printf("***antibody_sequence: ");
    for (double element : antibody_copy) {
        printf("{%g}", element);
    }
    printf("***.\n");


    // find smaller slide sequence and pad the longer one
    std::vector<double> padded_sequence {};
    std::vector<double> slide_sequence {};

    int i_antigen = antigen_copy.size();
    int i_antibody = antibody_copy.size();

    if (i_antigen <= i_antibody) {
        slide_sequence = antigen_copy;
        for (double element : antigen_copy) {
            padded_sequence.push_back(PAD);
        }
        for (double element : antibody_copy) {
            padded_sequence.push_back(element);
        }
        for (double element : antigen_copy) {
            padded_sequence.push_back(PAD);
        }
    }
    else {
        slide_sequence = antigen_copy;
        for (double element : antibody_copy) {
            padded_sequence.push_back(PAD);
        }
        for (double element : antigen_copy) {
            padded_sequence.push_back(element);
        }
        for (double element : antibody_copy) {
            padded_sequence.push_back(PAD);
        }
    }

    // print padded sequence
    printf("***padded_sequence: ");
    for (double element : padded_sequence) {
        printf("{%d}", (int) element);
    }
    printf("***.\n");
    printf("***side_sequence: ");
    for (double element : slide_sequence) {
        printf("{%d}", (int) element);
    }
    printf("***.\n");


    // get hamming distance.
    double i_hammdist_max = 0.0;
    int i_slide = slide_sequence.size();
    int i_template = padded_sequence.size() - i_slide;

    for (size_t i=0; i <= i_template; i++) {
        double i_hammdist = 0.0;
        printf("***sequence intersection: ");
        for (size_t j=i; j < (i + i_slide); j++) {
            if (padded_sequence[j] == slide_sequence[j-i]) {
                i_hammdist = i_hammdist + 1.0;
            }
            // print slide sequences
            printf("{%d}", (int) padded_sequence[j]);
        }
        printf("*** hamming distance: %g.\n", i_hammdist);
        if (i_hammdist_max < i_hammdist) {
            i_hammdist_max = i_hammdist;
        }
    }

    // calcualte hammingdistance score.
    if (i_slide < aminoComplete) {
        printf("Warning : aminoComplete {%d} is greater than the smaller squence {%d}, hamming score can never reach 1!\n", (int) aminoComplete, i_slide);
    }

    double r_hammscore =  i_hammdist_max / aminoComplete;
    if (r_hammscore > 1) {
        r_hammscore = 1;
    }

    // output
    return(r_hammscore);
}
