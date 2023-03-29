////////
// title: hamming.cpp
//
// language: C++
// date: 2023-03
// license: BSD 3-Clause
//
// description:
//   custome functions used in custom.cpp
////////

#include <cstdio>
#include <vector>
#include "./custom.h"
#include "./hamming.h"


// get vector_variable, check thereby if the variable name exist
std::vector<double> get_vector_variable( Cell* pCell, std::string name ) {
    int index = pCell->custom_data.find_vector_variable_index(name);
    if (index < 0 || index >= pCell->custom_data.vector_variables.size())
        throw std::invalid_argument("The cell has no vector with name `"+ name +"`");
    return pCell->custom_data.vector_variables[index].value;
}


// specify equal probabilities to choose from ALPHABET with choose_event
size_t alphabetLength = sizeof(ALPHABET) / sizeof(ALPHABET[0]);
std::vector<double> probabilities(alphabetLength, 1.0 / alphabetLength);


// generate antibody/antigen sequence given actual length of coding sequence, and length of vector
std::vector<double> generateSequence( int lenSequence ) {
    int lenPad = LEN_VECTOR_SEQUENCE - lenSequence;
    std::vector<double> sequence;
    for (size_t i = 0; i < lenSequence; ++i) {
        sequence.push_back(ALPHABET[choose_event(probabilities)]);
    }
    for (size_t i = 0; i < lenPad; ++i) {
        sequence.push_back(PAD);
    }
    printSequence(sequence);
    return sequence;
}


// Mutate antibody/antigen sequence given sequence and number of mutations (doesn't affect padding)
// This function can mutate the same letter twice, which can be fixed by changing mutateProbabilities after each mutation
// May change if necessary
void mutateSequence( std::vector<double>& sequence, int mutations ) {
    // Find length of actual sequence without padding
    size_t sequenceLength = 0;
    for (const auto& value : sequence) {
        if (value != PAD) {
            ++sequenceLength;
        } else {
            break;
        }
    }

    // Create equal probability vector to choose from the sequence
    std::vector<double> mutateProbabilities( sequenceLength, 1.0 / sequenceLength );

    // Mutate the sequence
    for (size_t i = 0; i < mutations; ++i) {
        sequence[choose_event(mutateProbabilities)] = ALPHABET[choose_event(probabilities)];
    }
}


// Print antibody/antigen sequence
void printSequence( std::vector<double>& sequence ) {
    printf("sequence: ");
    for (const auto& element : sequence) {
        if ((int)element == PAD) {
            printf(" %d", (int)element);
        } else {
            printf(" %c", (int)element);
        }
    }
    printf("\n");
}


// qunatify hamming distance from to antigen antibody sequences of variable size.
double alignment( Vector_Variable antigenSequence, Vector_Variable antibodySequence ) {

    // strip input
    std::vector<double> antigenCode {};
    std::vector<double> antibodyCode {};

    for (double element : antigenSequence.value) {
        if (element != PAD) {
            antigenCode.push_back(element);
        }
    }
    for (double element : antibodySequence.value) {
        if (element != PAD) {
            antibodyCode.push_back(element);
        }
    }

    // print seqeunces
    printSequence(antigenCode);
    printSequence(antibodyCode);

    // find smaller slide sequence and pad the longer one
    std::vector<double> paddedSequence {};
    std::vector<double>* topadSequence {};
    std::vector<double> slideSequence {};

    int i_antigen = antigenCode.size();
    int i_antibody = antibodyCode.size();

    if (i_antigen <= i_antibody) {
        slideSequence = antigenCode;
        topadSequence = &antibodyCode;
    }
    else {
        slideSequence = antibodyCode;
        topadSequence = &antigenCode;
    }

    for (double element : slideSequence) {
        paddedSequence.push_back(PAD);
    }
    for (double element : *topadSequence) {
        paddedSequence.push_back(element);
    }
    for (double element : slideSequence) {
        paddedSequence.push_back(PAD);
    }

    // print seqeunces
    printSequence(paddedSequence);
    printSequence(slideSequence);

    // get hamming distance.
    double i_hammdist_max = 0.0;
    int i_slide = slideSequence.size();
    int i_padded = paddedSequence.size() - i_slide;

    for (size_t i=0; i <= i_padded; i++) {
        double i_hammdist = 0.0;
        printf("***sequence intersection: ");
        for (size_t j=i; j < (i + i_slide); j++) {
            if (paddedSequence[j] == slideSequence[j-i]) {
                i_hammdist = i_hammdist + 1.0;
            }
            // print slide sequences
            if ((int) paddedSequence[j] == 0) {
                printf(" %d", (int) paddedSequence[j]);
            } else {
                printf(" %c", (int) paddedSequence[j]);
            }
        }
        printf("*** hamming distance: %g.\n", i_hammdist);
        if (i_hammdist_max < i_hammdist) {
            i_hammdist_max = i_hammdist;
        }
    }

    // calcualte hammingdistance score.
    if (i_slide < LEN_AMINOCOMPLETE) {
        printf("Warning : LEN_AMINOCOMPLETE {%d} is greater than the smaller squence {%d}, hamming score can never reach 1!\n", (int) LEN_AMINOCOMPLETE, i_slide);
    }

    double r_hammscore =  i_hammdist_max / LEN_AMINOCOMPLETE;
    if (r_hammscore > 1) {
        r_hammscore = 1;
    }

    // output
    return(r_hammscore);
}
