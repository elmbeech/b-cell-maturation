#include <cstdio>
#include <vector>

//b,j,o,u,x,z
double amminoComplete = 5;   // 100% affinity
double alphabet[] {'a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y'};
std::vector<double> antigenSequence {'c','a','d','d','c','e'};
std::vector<double> antibodySequence {'a','d','d','c','c','c','d','e','f','a'};

int main() {
    int i_antigen = antigenSequence.size();
    int i_antibody = antibodySequence.size();
    std::vector<double> paddedSequence {};
    std::vector<double> slideSequence {};

    // find smaller sequence and pad
    if (i_antigen <= i_antibody) {
        slideSequence = antigenSequence;
        for (double element : antigenSequence) {
            paddedSequence.push_back(NULL);
        }
        for (double element : antibodySequence) {
            paddedSequence.push_back(element);
        }
        for (double element : antigenSequence) {
            paddedSequence.push_back(NULL);
        }
    }
    else {
        slideSequence = antibodySequence;
        for (double element : antibodySequence) {
            paddedSequence.push_back(NULL);
        }
        for (double element : antigenSequence) {
            paddedSequence.push_back(element);
        }
        for (double element : antibodySequence) {
            paddedSequence.push_back(NULL);
        }
    }

    // get hamming distance.
    double i_hamm_max = 0;
    int i_slide = slideSequence.size();
    for (size_t i=i_slide; i < paddedSequence.size(); i++) {
        double i_hamm = 0;
        for (size_t j=i; j < i + i_slide; j++) {
            if (paddedSequence[i+j] == slideSequence[j]) {
                i_hamm = i_hamm + 1;
            }
        }
        if (i_hamm_max < i_hamm) {
            i_hamm_max = i_hamm;
        }
    }

    // calcualte hammingdistance score.
    double r_hamm =  i_hamm_max / amminoComplete;
    //if (r_hamm > 1);
    //    r_hamm = 1;
    //return(r_hamm);
    printf("hamming score: %g = %g / %g\n", r_hamm, i_hamm_max, amminoComplete);
}
