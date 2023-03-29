#include "../core/PhysiCell.h"

using namespace PhysiCell;

double alignment( Vector_Variable antigenSequence, Vector_Variable antibodySequence );
std::vector<double> generateSequence( int lenSequence );
std::vector<double> get_vector_variable( Cell* pCell, std::string name );
void mutateSequence( std::vector<double>& sequence, int mutations );
void printSequence( std::vector<double>& sequence );
