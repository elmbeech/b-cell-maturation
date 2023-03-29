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
#include <stdlib.h>
#include <time.h>
#include <typeinfo>
#include <vector>

#include "./custom.h"
#include "./hamming.h"

//Cell_Definition* invader;  //TODO
Cell_Definition* tfhelper_cell;
Cell_Definition* bnaive_cell;
Cell_Definition* bfollicular_cell;
Cell_Definition* bplasma_cell;
//Cell_Definition* bmemory_cell;
//Cell_Definition* antibody;  //TODO


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
        create_bplasma_cell_type();
        //create_bmemory_cell_type();  // BUE: only exist through transition
	//create_antibody_type();  // BUE: only exist through transition

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


// follicular T helper cell
void create_tfhelper_cell_type( void )  {
	tfhelper_cell = find_cell_definition( "Tf_helper" );

	std::vector<double> antigenSequence = generateSequence( LEN_ANTIGEN_SEQUENCE );
	tfhelper_cell->custom_data.add_vector_variable( "antigenSequence", antigenSequence );

	std::vector<double> antibodySequence = EMPTY_VECTOR;
	tfhelper_cell->custom_data.add_vector_variable( "antibodySequence", antibodySequence);

	tfhelper_cell->functions.update_phenotype = tfhelper_cell_phenotype;
}

void tfhelper_cell_phenotype( Cell* pCell, Phenotype& phenotype , double dt ) {
	std::vector<double> foreignAntigen = get_vector_variable(pCell, "antigenSequence");

	int numTouching = pCell->state.neighbors.size();
	for (int i = 0; i < numTouching; i++) {
		Cell* neighbor = pCell->state.neighbors[i];
		if (neighbor->type_name == B_NAIVE_NAME) {
			pCell->is_movable = false;
			pCell->attach_cell(neighbor);
			pCell->functions.update_phenotype = NULL;

			//Transfer antigen to B cell
			int antigenIndex = neighbor->custom_data.find_vector_variable_index("antigenSequence");
			neighbor->custom_data.vector_variables[antigenIndex].value = foreignAntigen;

			// set_single_behavior(neighbor, "transform to B_follicular", 1e9); //FIXME
                        // BUE 20230328: I tink the naive B cell should initate its transfromation, not the Tf helper cell
			return;
		}
	}

	// increase migration bias in higher quorum factor
	double q = get_single_signal( pCell, "Quorum_factor");
	double b0 = get_single_base_behavior( pCell, "migration bias");
	double bM = 1;
	double b = b0 + (bM-b0)*linear_response_function( q , 0 , 1 );
	set_single_behavior( pCell , "migration bias" , b );

	// reduce migration speed in higher quorum factor
	// double s0 = get_single_base_behavior( pCell, "migration speed");
	// double sM = 0.1*s0;
	// double s = s0 + (sM-s0)*decreasing_linear_response_function( q , 0.5, 0.75 );
	// set_single_behavior( pCell , "migration speed" , s );
}


// naive B cell
void create_bnaive_cell_type( void )  {
        // from xml
	bnaive_cell = find_cell_definition( "B_naive" );

        // antigen variablae
	std::vector<double> antigenSequence = EMPTY_VECTOR;
	bnaive_cell->custom_data.add_vector_variable( "antigenSequence", antigenSequence );

        // antibody variable
	std::vector<double> antibodySequence = EMPTY_VECTOR;
	bnaive_cell->custom_data.add_vector_variable( "antibodySequence", antibodySequence );

	// update phenotype
	bnaive_cell->functions.update_phenotype = bnaive_cell_phenotype;
}

void bnaive_cell_phenotype( Cell* pCell, Phenotype& phenotype , double dt ) {

        // antigen sequence
	int antigenIndex = pCell->custom_data.find_vector_variable_index("antigenSequence");
	Vector_Variable antigenSequence = pCell->custom_data.vector_variables[antigenIndex];


	// shodul I transform to follicular B cell?
	if ( antigenSequence.value != EMPTY_VECTOR ) {
	    printf("\nYay, got antigen, transform to B_follicular cell!\n");

            // antibody sequence
	    int antibodyIndex = pCell->custom_data.find_vector_variable_index("antibodySequence");
	    Vector_Variable antibodySequence = pCell->custom_data.vector_variables[antibodyIndex];
	    antibodySequence.value = generateSequence( LEN_ANTIBODY_SEQUENCE );

            // phenotype
	    set_single_behavior(pCell, "transform to B_follicular", 9e9);

            printSequence(antigenSequence.value);
            printSequence(antibodySequence.value);
        }
}


// follicular B cell
void create_bfollicular_cell_type( void )  {
        // from xml
	bfollicular_cell = find_cell_definition( "B_follicular" );

        // antigen variablae
	std::vector<double> antigenSequence = EMPTY_VECTOR;
	bfollicular_cell->custom_data.add_vector_variable( "antigenSequence", antigenSequence );

        // antibody variable
	std::vector<double> antibodySequence = EMPTY_VECTOR;
	bfollicular_cell->custom_data.add_vector_variable( "antibodySequence", antibodySequence );

	// update phenotype
	bfollicular_cell->functions.update_phenotype = bfollicular_cell_phenotype;
}


void bfollicular_cell_phenotype( Cell* pCell, Phenotype& phenotype , double dt ) {

	printf("I am a B_follicular cell!\n");
        // load antigen and anibody sequence
	int antigenIndex = pCell->custom_data.find_vector_variable_index("antigenSequence");
	Vector_Variable antigenSequence = pCell->custom_data.vector_variables[antigenIndex];
        printSequence( antigenSequence.value);

	int antibodyIndex = pCell->custom_data.find_vector_variable_index("antibodySequence");
	Vector_Variable antibodySequence = pCell->custom_data.vector_variables[antibodyIndex];
        printSequence( antibodySequence.value);

	// mutate inherited antibody sequnece by cell division!!!
        mutateSequence( antibodySequence.value, MUTATION );

        // get alignment signal
        double hammscore = alignment ( antigenSequence,  antibodySequence );
        printf("alignment hamming score: %g\n", hammscore);

	// shodul I transform to a plasma or a memory B cell?
	if ( hammscore > 0.9) {
	    set_single_behavior(pCell, "transform to B_plasma", 9e9);
        }
	else{
	    printf("Alignment not specific enough, stay follicular.\n");
	}

        // get pressure signal
        double pressure = get_single_signal( pCell , "pressure");
        // min pressure will be 0 [?]
        // max pressure - I have no idea. pragmatically set to 10. [?]
        double s0Pressure {0.0};
        double s1Pressure {10.0};
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

// plasma B cell
void create_bplasma_cell_type( void )  {
        // from xml
	bplasma_cell = find_cell_definition( "B_plasma" );

        // antigen variablae
	std::vector<double> antigenSequence = EMPTY_VECTOR;
	bplasma_cell->custom_data.add_vector_variable( "antigenSequence", antigenSequence );

        // antibody variable
	std::vector<double> antibodySequence = EMPTY_VECTOR;
	bplasma_cell->custom_data.add_vector_variable( "antibodySequence", antibodySequence );

	// update phenotype
	//bplasma_cell->functions.update_phenotype = bplasma_cell_phenotype;
}

