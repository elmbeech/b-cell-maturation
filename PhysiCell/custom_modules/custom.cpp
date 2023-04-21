/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2021, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/
#include <cstdio>
#include <stdlib.h>
#include <string>
#include <time.h>
#include <typeinfo>
#include <vector>

#include "./custom.h"


//Cell_Definition* invader;  //TODO
Cell_Definition* tfhelper_cell;
Cell_Definition* bnaive_cell;
Cell_Definition* bfollicular_cell;
//Cell_Definition* bmemory_cell;  //TODO
Cell_Definition* bplasma_cell;
//Cell_Definition* antibody;  //TODO


// custom constantes and variables

//letters not in the human amino acid alphabet: b,j,o,u,x,z
//static const double ALPHABET[] {'a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y'};  // humman amino acide alphabet
static const double ALPHABET[] {'a','t','c','g'};  // oligo nucleotide alphabet
static const double PAD = 0;

// LEN_VECTOR_SEQUENCE >= max(LEN_ANTIGEN_SEQUENCE, LEN_ANTIBODY_SEQUENCE) >= min(LEN_ANTIGEN_SEQUENCE, LEN_ANTIBODY_SEQUENCE) >= LEN_AMINOCOMPLETE
static const int LEN_VECTOR_SEQUENCE = 24;
static const int LEN_ANTIBODY_SEQUENCE = 24;
static const int LEN_ANTIGEN_SEQUENCE = 24;
static const int LEN_AMINOCOMPLETE = 16;  // number of matching antigen antibody amino sequences that account for 100% affinity.
static const double MUTATION_PER_SEQUENCE = 3.0;  // number antibody sequence mutations per follicular B cell division.
static const double MUTATION_CHANCE = 0.5;  // about every daugther cell mutates.
static const std::vector<double> EMPTY_VECTOR (LEN_VECTOR_SEQUENCE, PAD);  // generate empty antigen antybody vector.

// gereate sequence manual
//static const std::vector<double> EMPTY_VECTOR {PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD,PAD};
//static const std::vector<double> AGENSEQ_VECTOR {'a','t','c','g','a','a','t','t','c','c','g','g','a','t','c','g'};


// specify equal probabilities to choose from ALPHABET with choose_event
size_t alphabetLength = sizeof(ALPHABET) / sizeof(ALPHABET[0]);
std::vector<double> probabilities(alphabetLength, 1.0 / alphabetLength);


void create_cell_types( void )
{
    // set the random seed
    SeedRandom( parameters.ints("random_seed") );

    /*
       Put any modifications to default cell definition here if you
       want to have "inherited" by other cell types.

       This is a good place to set default functions.
    */

    initialize_default_cell_definition();
    cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );

    cell_defaults.functions.volume_update_function = standard_volume_update_function;
    cell_defaults.functions.update_velocity = standard_update_cell_velocity;

    cell_defaults.functions.update_migration_bias = NULL;
    cell_defaults.functions.update_phenotype = NULL; // update_cell_and_death_parameters_O2_based;
    cell_defaults.functions.custom_cell_rule = NULL;
    cell_defaults.functions.contact_function = NULL;

    cell_defaults.functions.add_cell_basement_membrane_interactions = NULL;
    cell_defaults.functions.calculate_distance_to_membrane = NULL;

    /*
       This parses the cell definitions in the XML config file.
    */

    initialize_cell_definitions_from_pugixml();

    /*
       This builds the map of cell definitions and summarizes the setup.
    */

    build_cell_definitions_maps();

    /*
       This intializes cell signal and response dictionaries
    */

    setup_signal_behavior_dictionaries();

    /*
       Put any modifications to individual cell definitions here.

       This is a good place to set custom functions.
    */

    cell_defaults.functions.update_phenotype = phenotype_function;
    cell_defaults.functions.custom_cell_rule = custom_function;
    cell_defaults.functions.contact_function = contact_function;

    //create_invader_type();   // TODO
    create_tfhelper_cell_type();
    create_bnaive_cell_type();
    create_bfollicular_cell_type();
    //create_bmemory_cell_type();  // TODO
    create_bplasma_cell_type();
    //create_antibody_type();  // TODO

    /*
       This builds the map of cell definitions and summarizes the setup.
    */

    display_cell_definitions( std::cout );

    return;
}

void setup_microenvironment( void )
{
    // set domain parameters

    // put any custom code to set non-homogeneous initial conditions or
    // extra Dirichlet nodes here.
    srand(time(NULL));

    // initialize BioFVM

    initialize_microenvironment();

    return;
}

void setup_tissue( void )
{
    double Xmin = microenvironment.mesh.bounding_box[0];
    double Ymin = microenvironment.mesh.bounding_box[1];
    double Zmin = microenvironment.mesh.bounding_box[2];

    double Xmax = microenvironment.mesh.bounding_box[3];
    double Ymax = microenvironment.mesh.bounding_box[4];
    double Zmax = microenvironment.mesh.bounding_box[5];

    if( default_microenvironment_options.simulate_2D == true )
    {
        Zmin = 0.0;
        Zmax = 0.0;
    }

    double Xrange = Xmax - Xmin;
    double Yrange = Ymax - Ymin;
    double Zrange = Zmax - Zmin;

    // create some of each type of cell

    Cell* pC;

    for( int k=0; k < cell_definitions_by_index.size() ; k++ )
    {
        Cell_Definition* pCD = cell_definitions_by_index[k];
        std::cout << "Placing cells of type " << pCD->name << " ... " << std::endl;
        for( int n = 0 ; n < parameters.ints("number_of_cells") ; n++ )
        {
            std::vector<double> position = {0,0,0};
            position[0] = Xmin + UniformRandom()*Xrange;
            position[1] = Ymin + UniformRandom()*Yrange;
            position[2] = Zmin + UniformRandom()*Zrange;

            pC = create_cell( *pCD );
            pC->assign_position( position );
        }
    }
    std::cout << std::endl;

    // load cells from your CSV file (if enabled)
    load_cells_from_pugixml();

    return;
}

std::vector<std::string> my_coloring_function( Cell* pCell )
{ return paint_by_number_cell_coloring(pCell); }

void phenotype_function( Cell* pCell, Phenotype& phenotype, double dt )
{ return; }

void custom_function( Cell* pCell, Phenotype& phenotype , double dt )
{ return; }

void contact_function( Cell* pMe, Phenotype& phenoMe , Cell* pOther, Phenotype& phenoOther , double dt )
{ return; }


// get vector_variable, check thereby if the variable name exist
std::vector<double> get_vector_variable( Cell* pCell, std::string name ) {
    int index = pCell->custom_data.find_vector_variable_index(name);
    if (index < 0 || index >= pCell->custom_data.vector_variables.size())
        throw std::invalid_argument("The cell has no vector with name `"+ name +"`");
    return pCell->custom_data.vector_variables[index].value;
}


// follicular T helper cell
void create_tfhelper_cell_type(void)  {
    tfhelper_cell = find_cell_definition("Tf_helper");

    // custom vector variable antigen
    std::vector<double> antigenSequence = generateSequence(LEN_ANTIGEN_SEQUENCE);
    tfhelper_cell->custom_data.add_vector_variable("antigenSequence", antigenSequence);
    // custom vector variable antibody
    std::vector<double> antibodySequence = EMPTY_VECTOR;
    tfhelper_cell->custom_data.add_vector_variable("antibodySequence", antibodySequence);
    // custom vector variable anker
    std::vector<double> coordinateAnchor = {PAD, PAD, PAD};
    tfhelper_cell->custom_data.add_vector_variable("coordinateAnchor", coordinateAnchor);
    // custom variables
    tfhelper_cell->custom_data.add_variable("mutate", -1.0);
    tfhelper_cell->custom_data.add_variable("hamming_fract", -1.0);
    tfhelper_cell->custom_data.add_variable("pressure_fract", -1.0);
    tfhelper_cell->custom_data.add_variable("apoptosis", -1.0);
    tfhelper_cell->custom_data.add_variable("apoptosis_fract", -1.0);
    tfhelper_cell->custom_data.add_variable("b_anchor", 0.0);

    // update phenotype
    tfhelper_cell->functions.update_phenotype = tfhelper_cell_phenotype;
    tfhelper_cell->functions.custom_cell_rule = tfhelper_cell_custom;
}

// executes evert biological timestep!
void tfhelper_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt) {
    std::vector<double> foreignAntigen = get_vector_variable(pCell, "antigenSequence");
    //int antigenIndex = pCell->custom_data.find_vector_variable_index("antigenSequence");
    //Vector_Variable antigenSequence = pCell->custom_data.vector_variables[antigenIndex];

    int coordinateAnchorIndex = pCell->custom_data.find_vector_variable_index("coordinateAnchor");

    int numTouching = pCell->state.neighbors.size();
    for (int i = 0; i < numTouching; i++) {
        Cell* neighbor = pCell->state.neighbors[i];
        if (neighbor->type_name == "B_naive") {
            // manipulate Tf helper cell
            //pCell->attach_cell(neighbor);
            //pCell->is_movable = false;
            //pCell->functions.update_phenotype = NULL;
            //neighbor->is_movable = false;

            // Anchor the cell
            //  Here we just get the position for the anchor and set b_anchor to 1
            pCell->phenotype.motility.is_motile = false;
            pCell->custom_data["b_anchor"] = 1.0;
            pCell->custom_data.vector_variables[coordinateAnchorIndex].value = pCell->position;

            //transfer antigen to B cell
            int antigenIndex = neighbor->custom_data.find_vector_variable_index("antigenSequence");
            neighbor->custom_data.vector_variables[antigenIndex].value = foreignAntigen;
            // break out of the for loop
            break;
        }
    }

    // increase migration bias in higher quorum factor
    double q = get_single_signal( pCell, "Quorum_factor");
    double b0 = get_single_base_behavior( pCell, "migration bias");
    double bM = 1;
    double b = b0 + (bM-b0)*linear_response_function( q , 0 , 1 );
    set_single_behavior( pCell , "migration bias" , b );

    // reduce migration speed in higher quorum factor
    //double s0 = get_single_base_behavior( pCell, "migration speed");
    //double sM = 0.1*s0;
    //double s = s0 + (sM-s0)*decreasing_linear_response_function( q , 0.5, 0.75 );
    //set_single_behavior( pCell , "migration speed" , s );
}

// executes every mechanical timestep!
void tfhelper_cell_custom(Cell* pCell, Phenotype& phenotype , double dt) {
    int coordinateAnchorIndex = pCell->custom_data.find_vector_variable_index("coordinateAnchor");
    pCell->velocity -= pCell->custom_data["b_anchor"] * (pCell->position - pCell->custom_data.vector_variables[coordinateAnchorIndex].value);
}


// naive B cell
void create_bnaive_cell_type(void)  {
    // from xml
    bnaive_cell = find_cell_definition("B_naive");

    // custom vector variable antigen
    std::vector<double> antigenSequence = EMPTY_VECTOR;
    bnaive_cell->custom_data.add_vector_variable("antigenSequence", antigenSequence);
    // custom vector variable antibody
    std::vector<double> antibodySequence = EMPTY_VECTOR;
    bnaive_cell->custom_data.add_vector_variable("antibodySequence", antibodySequence);
    // custom vector variable anker
    std::vector<double> coordinateAnchor = { PAD, PAD, PAD };
    bnaive_cell->custom_data.add_vector_variable("coordinateAnchor", coordinateAnchor);
    // custom variables
    bnaive_cell->custom_data.add_variable("mutate", -1.0);
    bnaive_cell->custom_data.add_variable("hamming_fract", -1.0);
    bnaive_cell->custom_data.add_variable("pressure_fract", -1.0);
    bnaive_cell->custom_data.add_variable("apoptosis", -1.0);
    bnaive_cell->custom_data.add_variable("apoptosis_fract", -1.0);

    // update phenotype
    bnaive_cell->functions.update_phenotype = bnaive_cell_phenotype;
}

void bnaive_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt) {

    // antigen sequence
    int antigenIndex = pCell->custom_data.find_vector_variable_index("antigenSequence");
    Vector_Variable antigenSequence = pCell->custom_data.vector_variables[antigenIndex];

    // shodul I transform to follicular B cell?
    if (antigenSequence.value != EMPTY_VECTOR) {
        // antibody sequence
        int antibodyIndex = pCell->custom_data.find_vector_variable_index("antibodySequence");
        pCell->custom_data.vector_variables[antibodyIndex].value = generateSequence(LEN_ANTIBODY_SEQUENCE);

        // phenotype
        set_single_behavior(pCell, "transform to B_follicular", 9e9);

        // print
        printf("\nCell ID %d: Yay, got antigen, transform to B_follicular cell!\n", pCell->ID);
        printSequence(antigenSequence.value, "Antigen: ");
        Vector_Variable antibodySequence = pCell->custom_data.vector_variables[antibodyIndex];
        printSequence(antibodySequence.value, "Antibody: ");
    }
}


// follicular B cell
void create_bfollicular_cell_type(void)  {
    // from xml
    bfollicular_cell = find_cell_definition("B_follicular");

    // custem vector variable antigen
    std::vector<double> antigenSequence = EMPTY_VECTOR;
    bfollicular_cell->custom_data.add_vector_variable("antigenSequence", antigenSequence);
    // custom vector variable antibody
    std::vector<double> antibodySequence = EMPTY_VECTOR;
    bfollicular_cell->custom_data.add_vector_variable("antibodySequence", antibodySequence);
    // custom vector variable anker
    std::vector<double> coordinateAnchor = { PAD, PAD, PAD };
    bfollicular_cell->custom_data.add_vector_variable("coordinateAnchor", coordinateAnchor);
    // custom variables
    bfollicular_cell->custom_data.add_variable("mutate", -1.0);
    bfollicular_cell->custom_data.add_variable("hamming_fract", -1.0);
    bfollicular_cell->custom_data.add_variable("pressure_fract", -1.0);
    bfollicular_cell->custom_data.add_variable("apoptostic", -1.0);
    bfollicular_cell->custom_data.add_variable("apoptosis_fract", -1.0);

    // update phenotype
    bfollicular_cell->functions.update_phenotype = bfollicular_cell_phenotype;
}

void bfollicular_cell_phenotype(Cell* pCell, Phenotype& phenotype , double dt) {
    // if cell not in apoptosis mode
    if (get_single_signal(pCell, "dead") < 0.5) {

        // load antigen and anibody sequence
        int antigenIndex = pCell->custom_data.find_vector_variable_index("antigenSequence");
        Vector_Variable antigenSequence = pCell->custom_data.vector_variables[antigenIndex];

        int antibodyIndex = pCell->custom_data.find_vector_variable_index("antibodySequence");
        Vector_Variable antibodySequence = pCell->custom_data.vector_variables[antibodyIndex];

        // check for cell cycle time and if divided less than 1[h] ago and possibly mutate
        double generationTime = phenotype.cycle.data.elapsed_time_in_phase;
        if (generationTime < 60.0) {  // elapsed_time_in_phase is in [min]
            if (pCell->custom_data["mutate"]  < -0.5) {
                printf("Cell ID %d: just divided, %g[min] ago.\n", pCell->ID, generationTime);
                std::vector<double> flip {1 - MUTATION_CHANCE, MUTATION_CHANCE};
                int choice = choose_event(flip);
                if (choice == 0) { pCell->custom_data["mutate"]  = 0.0; }
                else { pCell->custom_data["mutate"] = MUTATION_PER_SEQUENCE; }
            }
            printf("Cell ID %d: mutation to go: %d\n", pCell->ID, (int)pCell->custom_data["mutate"]);
            if (pCell->custom_data["mutate"] > 0.5) {
                printSequence(antigenSequence.value, "Antigen: ");
                printSequence(antibodySequence.value, "Antibody: ");
                mutateSequence(antibodySequence.value);
                pCell->custom_data.vector_variables[antibodyIndex].value = antibodySequence.value;
                printSequence(antibodySequence.value, "Mutated : ");
                pCell->custom_data["mutate"] = pCell->custom_data["mutate"] - 1.0;
            }
            else {
                printSequence(antigenSequence.value, "Antigen: ");
                printSequence(antibodySequence.value, "Antibody: ");
                printf("No mutation!\n");
            }
        }
        else {
            pCell->custom_data["mutate"] = -1.0;
        }

        // get alignment signal
        double fracHamming = alignment(antigenSequence, antibodySequence, false);
        //int hammingIndex = pCell->custom_data.find_variable_index("hamming");
        //pCell->custom_data.variables[hammingIndex].value = hammscore;
        pCell->custom_data["hamming_fract"] = fracHamming;

        // should I transform to a plasma or a memory B cell?
        if (fracHamming > 0.9) {
            printf("Yay, high hamming score, transform to B_plasma cell!\n");
            set_single_behavior(pCell, "transform to B_plasma", 9e9);
        }
        else {
            // get pressure signal
            double pressure = get_single_signal(pCell, "pressure");
            // get pressure respons
            double s0Pressure {0.0};  // min pressure
            double s1Pressure {10.0};  // max pressure evaluated by measurement
            double fracPressure = linear_response_function(pressure, s0Pressure, s1Pressure);  // value between 0 and 1
            pCell->custom_data["pressure_fract"] = fracPressure;

            // get and output apoptosis state value
            double apoptosis = get_single_behavior(pCell, "apoptosis");
            pCell->custom_data["apoptosis"] = apoptosis;

            // set and output apoptosis rate
            // default min rate value for apoptosis is 5.31667e-05 [1/min]
            // this means the probabiliry to die in a 60 min time setp is 60[min] * 5.31667e-05[1/min] = 0.003190002 (~ 3 per mille)
            double s0Apoptosis = get_single_base_behavior(pCell, "apoptosis" );  // min pressure
            double s1Apoptosis = 0.98 / 60;  // we specify, at max 98% of the cells should enter apoptosis within 60 [1/min], this specifies the steepness
            //double rApoptosis = s0Apoptosis + (s1Apoptosis - s0Apoptosis) * rPressure;
            //double rApoptosis = s0Apoptosis + (s1Apoptosis - s0Apoptosis) * (rPressure - fracHamming) / 2;
            double fracApoptosis = s0Apoptosis + (s1Apoptosis - s0Apoptosis) * fracPressure * (1 - fracHamming);
            set_single_behavior(pCell, "apoptosis" , fracApoptosis);
            pCell->custom_data["apoptosis_fract"] = fracApoptosis;
        }
    }
}


// plasma B cell
void create_bplasma_cell_type( void )  {
    // from xml
    bplasma_cell = find_cell_definition("B_plasma");

    // custom vector variable antigen
    std::vector<double> antigenSequence = EMPTY_VECTOR;
    bplasma_cell->custom_data.add_vector_variable("antigenSequence", antigenSequence);
    // custom vector variable antibody
    std::vector<double> antibodySequence = EMPTY_VECTOR;
    bplasma_cell->custom_data.add_vector_variable("antibodySequence", antibodySequence);
    // custom vector variable anker
    std::vector<double> coordinateAnchor = { PAD, PAD, PAD };
    bplasma_cell->custom_data.add_vector_variable("coordinateAnchor", coordinateAnchor);
    // custom variables
    bplasma_cell->custom_data.add_variable("mutate", -1.0);
    bplasma_cell->custom_data.add_variable("hamming_fract", -1.0);
    bplasma_cell->custom_data.add_variable("pressure_fract", -1.0);
    bplasma_cell->custom_data.add_variable("apoptosis", -1.0);
    bplasma_cell->custom_data.add_variable("apoptosis_fract", -1.0);

    // update phenotype
    bplasma_cell->functions.update_phenotype = NULL;
    //bplasma_cell->functions.update_phenotype = bplasma_cell_phenotype;  // TODO
}


// print antibody/antigen sequence
void printSequence( std::vector<double>& sequence, std::string prefix = "sequence: " ) {
    for (char element : prefix) printf("%c", element);
    for (double element : sequence) {
        if (element == PAD) {
            printf("{%d}", (int)element);
        } else {
            printf("{%c}", (int)element);
        }
    }
    printf("\n");
}


// generate antibody/antigen sequence given actual length of coding sequence, and length of vector
std::vector<double> generateSequence(int lenSequence) {
    int lenPad = LEN_VECTOR_SEQUENCE - lenSequence;
    std::vector<double> sequence;
    for (size_t i = 0; i < lenSequence; ++i) {
        sequence.push_back(ALPHABET[choose_event(probabilities)]);
    }
    for (size_t i = 0; i < lenPad; ++i) {
        sequence.push_back(PAD);
    }
    //printSequence(sequence, "Generated sequence: ");
    return sequence;
}


// Mutate antibody/antigen sequence at one random non pad character position.
void mutateSequence(std::vector<double>& sequence) {
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
    std::vector<double> mutateProbabilities(sequenceLength, 1.0 / sequenceLength);
    // Mutate the sequence
    sequence[choose_event(mutateProbabilities)] = ALPHABET[choose_event(probabilities)];
}


// qunatify hamming distance from to antigen antibody sequences of variable size.
double alignment(Vector_Variable antigenSequence, Vector_Variable antibodySequence , bool verbose = false) {

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

    // print stripped seqeunces
    if (verbose) {
        printSequence(antigenCode, "Antigen code");
        printSequence(antibodyCode,  "Antibody code");
    }

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
    if (verbose) {
        printSequence(paddedSequence, "Padded sequence: ");
        printSequence(slideSequence, "Slide sequence: ");
    }

    // get hamming distance.
    double i_hammdist_max = 0.0;
    int i_slide = slideSequence.size();
    int i_padded = paddedSequence.size() - i_slide;

    for (size_t i=0; i <= i_padded; i++) {
        if (verbose) printf("Sequence intersection: ***");
        double i_hammdist = 0.0;
        for (size_t j=i; j < (i + i_slide); j++) {
            if (paddedSequence[j] == slideSequence[j-i]) {
                i_hammdist = i_hammdist + 1.0;
            }
            if (verbose) {
                if ((int) paddedSequence[j] == PAD) printf("{%d}", (int) paddedSequence[j]);
                else printf("{%c}", (int) paddedSequence[j]);
            }
        }
        if (i_hammdist_max < i_hammdist) {
            i_hammdist_max = i_hammdist;
        }
        if (verbose) printf("*** hamming distance: %g.\n", i_hammdist);
    }

    // calcualte hamming distance score.
    if (i_slide < LEN_AMINOCOMPLETE) {
        printf("Warning : LEN_AMINOCOMPLETE {%d} is greater than the smaller squence {%d}, hamming score can never reach 1!\n", (int) LEN_AMINOCOMPLETE, i_slide);
    }

    double r_hammscore = i_hammdist_max / LEN_AMINOCOMPLETE;
    if (r_hammscore > 1) {
        r_hammscore = 1;
    }

    // output
    if (verbose) printf("Sequence intersection: ***");
    //printf("hamming distance score: %g.\n", r_hammscore);
    return(r_hammscore);
}
