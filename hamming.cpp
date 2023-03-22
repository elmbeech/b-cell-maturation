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

    printf("***antigenSequence: ");
    for (double element : antigenSequence) {                         
        printf("{%g}", element);                                            
    }                                                                       
    printf("***.\n");

    printf("***antibodySequence: ");
    for (double element : antibodySequence) {                         
        printf("{%g}", element);                                            
    }                                                                       
    printf("***.\n");

    // find smaller slide sequence and pad the longer one
    if (i_antigen <= i_antibody) {
        slideSequence = antigenSequence;
        for (double element : antigenSequence) {
            paddedSequence.push_back(0);
        }
        for (double element : antibodySequence) {
            paddedSequence.push_back(element);
        }
        for (double element : antigenSequence) {
            paddedSequence.push_back(0);
        }
    }
    else {
        slideSequence = antigenSequence;
        for (double element : antibodySequence) {
            paddedSequence.push_back(0);
        }
        for (double element : antigenSequence) {
            paddedSequence.push_back(element);
        }
        for (double element : antibodySequence) {
            paddedSequence.push_back(0);
        }
    }

    // print padded sequence
    printf("***paddedSequence: ");
    for (double element : paddedSequence) {                         
        printf("{%d}", (int) element);                                            
    }                                                                       
    printf("***.\n");
    printf("***sideSequence: ");
    for (double element : slideSequence) {                         
        printf("{%d}", (int) element);                                            
    }                                                                       
    printf("***.\n\n");


    // get hamming distance.
    double i_hammdist_max = 0.0;
    int i_slide = slideSequence.size();
    int i_template = paddedSequence.size() - i_slide;

    for (size_t i=0; i <= i_template; i++) {
        double i_hammdist = 0.0;
        printf("***sequence intersection: ");
        for (size_t j=i; j < (i + i_slide); j++) {
            if (paddedSequence[j] == slideSequence[j-i]) {
                i_hammdist = i_hammdist + 1.0;
            }
            // print slide sequences
            printf("{%d}", (int) paddedSequence[j]);
        }
        printf("*** hamming distance: %g.\n\n", i_hammdist);
        if (i_hammdist_max < i_hammdist) {
            i_hammdist_max = i_hammdist;
        }
    }

    // calcualte hammingdistance score.
    double r_hammscore =  i_hammdist_max / amminoComplete;
    //if (r_hammscore > 1);
    //    r_hammscore = 1;
    //return(r_hammscore);
    printf("final hamming score: %g = %g / %g\n", r_hammscore, i_hammdist_max, amminoComplete);
}
