// MyoSimCpp.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <string>

//#include "gsl_math.h"

// Includes
#include "muscle.h"

int main(int argc, char* argv[])
{
    /**
    Main function
    + the entry point for MyoSimCpp
    **/

    // Variables
    std::string model_file_string;
    std::string options_file_string;
    std::string protocol_file_string;
    std::string results_file_string;
    std::string system_id;

    muscle* p_muscle;                       // Pointer to a muscle object

    printf("no_of_args: %i\n\n", argc);

    // Set inputs
    model_file_string = argv[1];
    options_file_string = argv[2];
    protocol_file_string = argv[3];
    results_file_string = argv[4];

    // Initialize
    p_muscle = new muscle(model_file_string, options_file_string);

    p_muscle->implement_protocol(protocol_file_string, results_file_string);

    // Tidy up
    delete p_muscle;

    // Close
    printf("Closing CMyoSim\n");

    return(1);

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
