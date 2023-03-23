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
#include <typeinfo>
#include <vector>
#include "./custom.h"

Cell_Definition* naive_bcell;

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

        create_naive_bcell_type();
	
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



//letters not in the human amino acid alphabet: b,j,o,u,x,z
static const double alphabet[] {'a','c','d','e','f','g','h','i','k','l','m','n','p','q','r','s','t','v','w','y'};  // humman amino acide alphabet
static const double pad {0};

// number of matching antigen antibody amino sequences that account for 100% affinity.
static double aminoComplete = 10;

//temporary (equal) probabilities to choose from alphabet with choose_event
size_t alphabetLength = sizeof(alphabet) / sizeof(alphabet[0]);
std::vector<double> probabilities(alphabetLength, 1.0 / alphabetLength);

// Generate antibody/antigen sequence given actual length of code, and length of padding
std::vector<double> generateSequence( int sequenceLength, int padLength ) {
    std::vector<double> sequence;
    for (size_t i = 0; i < sequenceLength; ++i) {
        sequence.push_back(alphabet[choose_event(probabilities)]);
    }
    for (size_t i = 0; i < padLength; ++i) {
        sequence.push_back(pad);
    }
    return sequence;
}

// Mutate antibody/antigen sequence given sequence and number of mutations (doesn't affect padding)
// This function can mutate the same letter twice, which can be fixed by changing mutateProbabilities after each mutation
// May change if necessary
void mutateSequence( std::vector<double>& sequence, int mutations ) {
    // Find length of actual sequence without padding
    size_t sequenceLength = 0;
    for (const auto& value : sequence) {
        if (value != pad) {
            ++sequenceLength;
        } else {
            break;
        }
    }

    // Create equal probability vector to choose from the sequence
    std::vector<double> mutateProbabilities(sequenceLength, 1.0 / sequenceLength);

    // Mutate the sequence
    for (size_t i = 0; i < mutations; ++i) {
        sequence[choose_event(mutateProbabilities)] = alphabet[choose_event(probabilities)];
    }
}

// Print antibody/antigen sequence
void printSequence( std::vector<double>& sequence ) {
    printf("sequence: ");
    for (const auto& value : sequence) {
        if ((int)value==0) {
            printf(" %d", (int)value);
        } else {
            printf(" %c", (int)value);
        }
    }
    printf("\n");
}

void create_naive_bcell_type( void )  {

	naive_bcell = find_cell_definition( "B_naive" );

        // antigen variable
	//std::vector<double> antigenSequence = {pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad,pad};
        std::vector<double> antigenSequence {'c','a','d','d','c','e','n','k','l','l','c',pad,pad,pad,pad,pad};
	naive_bcell->custom_data.add_vector_variable( "antigenSequence", antigenSequence );
        long antigenLength = antigenSequence.size();
        printf("number of antigenSequence ELEMENTS: %ld\n", antigenLength);

        // antibody variable
        std::vector<double> antibodySequence {'a','d','d','c','c','c','d','e','f','a','l','l','c','c','d','a'};
	naive_bcell->custom_data.add_vector_variable( "antibodySequence", antibodySequence );
        long antibodyLength = antibodySequence.size();
        printf("number of antibodySequence ELEMENTS: %ld\n", antibodyLength);

        // update phenotype
	naive_bcell->functions.update_phenotype = naive_bcell_phenotype; 
}

void naive_bcell_phenotype( Cell* pCell, Phenotype& phenotype , double dt ) {

        // antigen sequence
	int antigenIndex = pCell->custom_data.find_vector_variable_index("antigenSequence");
	Vector_Variable antigenSequence = pCell->custom_data.vector_variables[antigenIndex];
        printf("number of antigenSequence elements: %ld\n", antigenSequence.value.size());

        // antibody sequence
	int antibodyIndex = pCell->custom_data.find_vector_variable_index("antibodySequence");
	Vector_Variable antibodySequence = pCell->custom_data.vector_variables[antibodyIndex];
        printf("number of antibodySequence elements: %ld\n", antibodySequence.value.size());

        // get alignment signal
        double hammscore = alignment ( antigenSequence,  antibodySequence );
        printf("alignment hamming score: %g\n", hammscore);

        // get pressure signal
        double pressure = get_single_signal( pCell , "pressure");
        // min pressure will be 0 [?]
        // max pressure - I have no idea. pragmatically set to 10. [?]
        double s0Pressure {0.0};
        double s1Pressure {10000.0};
        // get the response functions
        double rPressure = linear_response_function( pressure, s0Pressure, s1Pressure);
        printf("pressure min: {%g}\tmax: {%g}\tdetected: {%g}\tresponse_fraction:{%g} \n", s0Pressure, s1Pressure, pressure, rPressure);

        // apoptosis response
        // get min rate value for apoptosis
        // the default rate is 5.31667e-05
        // this means the probability to die in a 6 min time step is 6[min] * 5.31667e-05[1/min] = 0.0003190002
        // this means the probabiliry to die in a 60 min time setp is 6[min] * 5.31667e-05[1/min] = 0.003190002 (~ 3 per mille)
        double s0Apoptosis = get_single_base_behavior( pCell, "apoptosis" );
        // get max rate value for apoptosis
        // let's say, at max pressure 98% of the cells should enter apoptosis within 60 [1/min]
        // 60[min] * rM[1/min] = 1
        double s1Apoptosis = 0.98 / 60;

        // bue 20130322: getting the complete formula adjusted will need some analysis. 
	// rule of pressure and alignment score to steer apoptosis rate
	//double rApoptosis = s0Apoptosis + (s1Apoptosis - s0Apoptosis) * (rPressure + (1 - hammscore)) / 2;
	double rApoptosis = s0Apoptosis + (s1Apoptosis - s0Apoptosis) * rPressure;
        set_single_behavior( pCell, "apoptosis" , rApoptosis );
        printf("apoptosis min: {%g}\tmax: {%g}\tset: {%g}\n", s0Apoptosis, s1Apoptosis, rApoptosis);
}

double alignment( Vector_Variable antigenSequence, Vector_Variable antibodySequence ) {

    // strip and print input
    std::vector<double> antigen_sequence {};
    std::vector<double> antibody_sequence {};

    for (double element : antigenSequence.value) {
        if (element != pad) {
            antigen_sequence.push_back(element);
        }
    }
    for (double element : antibodySequence.value) {
        if (element != pad) {
            antibody_sequence.push_back(element);
        }
    }

    printf("***antigen_sequence: ");
    for (double element : antigen_sequence) {
        printf("{%g}", element);
    }
    printf("***.\n");
    printf("***antibody_sequence: ");
    for (double element : antibody_sequence) {
        printf("{%g}", element);
    }
    printf("***.\n");


    // find smaller slide sequence and pad the longer one
    std::vector<double> padded_sequence {};
    std::vector<double> slide_sequence {};

    int i_antigen = antigen_sequence.size();
    int i_antibody = antibody_sequence.size();

    if (i_antigen <= i_antibody) {
        slide_sequence = antigen_sequence;
        for (double element : antigen_sequence) {
            padded_sequence.push_back(pad);
        }
        for (double element : antibody_sequence) {
            padded_sequence.push_back(element);
        }
        for (double element : antigen_sequence) {
            padded_sequence.push_back(pad);
        }
    }
    else {
        slide_sequence = antigen_sequence;
        for (double element : antibody_sequence) {
            padded_sequence.push_back(pad);
        }
        for (double element : antigen_sequence) {
            padded_sequence.push_back(element);
        }
        for (double element : antibody_sequence) {
            padded_sequence.push_back(pad);
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
