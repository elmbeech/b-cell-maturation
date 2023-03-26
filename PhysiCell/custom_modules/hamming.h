#include "../core/PhysiCell.h"

using namespace PhysiCell;

//letters not in the human amino acid alphabet: b,j,o,u,x,z
static const double ALPHABET[] {'a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y'};  // humman amino acide alphabet
static const double PAD {0};

// number of matching antigen antibody amino sequences that account for 100% affinity.
static double aminoComplete = 10;

double alignment( Vector_Variable antigenSequence, Vector_Variable antibodySequence );