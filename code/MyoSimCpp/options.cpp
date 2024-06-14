/**
/* @file		options.cpp
/* @brief		Source file for an options object
/* @author		Ken Campbell
*/

#include <cstdio>
#include <iostream>
#include <filesystem>
#include <string>

#include "global_definitions.h"
#include "JSON_functions.h"
#include "options.h"

#include "gsl_math.h"

#include "rapidjson\document.h"
#include "rapidjson\filereadstream.h"


// Constructor
options::options(std::string set_options_file_string)
{
    // Initialise

    // Code
    std::cout << "In options constructor\n";

    options_file_string = set_options_file_string;

    // Set some variables
    sprintf_s(rate_file_string, _MAX_PATH, "");

    // Set defining variables
    initialise_options_from_JSON_file_string(options_file_string);
}

// Destructor
options::~options(void)
{
    // Code
    std::cout << "In options\n";

    // Tidy up
}

void options::initialise_options_from_JSON_file_string(std::string JSON_options_file_string)
{
    //! Code initialises a model from a JSON formatted file


    errno_t file_error;
    FILE* fp;
    char readBuffer[65536];

    // Code
    file_error = fopen_s(&fp, JSON_options_file_string.c_str(), "rb");
    if (file_error != 0)
    {
        std::cout << "Error opening JSON model file: " << JSON_options_file_string;
        exit(1);
    }

    rapidjson::FileReadStream is(fp, readBuffer, sizeof(readBuffer));

    rapidjson::Document doc;
    doc.ParseStream(is);

    fclose(fp);

    // Now work through the JSON file
    JSON_functions::check_JSON_member_object(doc, "options");
    const rapidjson::Value& opt = doc["options"];

    JSON_functions::check_JSON_member_number(opt, "max_rate");
    max_rate = opt["max_rate"].GetDouble();

    JSON_functions::check_JSON_member_number(opt, "bin_x_min");
    bin_x_min = opt["bin_x_min"].GetDouble();

    JSON_functions::check_JSON_member_number(opt, "bin_x_max");
    bin_x_max = opt["bin_x_max"].GetDouble();

    JSON_functions::check_JSON_member_number(opt, "bin_x_width");
    bin_x_width = opt["bin_x_width"].GetDouble();

    if (JSON_functions::check_JSON_member_exists(opt, "movement_enhancement_width"))
    {
        movement_enhancement_width = opt["movement_enhancement_width"].GetDouble();
    }
    else
    {
        movement_enhancement_width = GSL_NAN;
    }

    // Myofibrils
    if (JSON_functions::check_JSON_member_exists(opt, "myofibrils"))
    {
        const rapidjson::Value& myof = opt["myofibrils"];

        JSON_functions::check_JSON_member_number(myof, "force_tolerance");
        myofibril_force_tolerance = myof["force_tolerance"].GetDouble();

        JSON_functions::check_JSON_member_number(myof, "max_iterations");
        myofibril_max_iterations = myof["max_iterations"].GetInt();

        JSON_functions::check_JSON_member_number(myof, "max_delta_hs_length");
        myofibril_max_delta_hs_length = myof["max_delta_hs_length"].GetDouble();
    }
    else
    {
        myofibril_force_tolerance = GSL_NAN;
        myofibril_max_iterations = 1;
        myofibril_max_delta_hs_length = GSL_NAN;
    }
    
    // Rate files
    if (JSON_functions::check_JSON_member_exists(opt, "rate_files"))
    {
        JSON_functions::check_JSON_member_object(opt, "rate_files");
        const rapidjson::Value& rf = opt["rate_files"];

        sprintf_s(rate_file_string, _MAX_PATH, rf["file"].GetString());
    }
}
