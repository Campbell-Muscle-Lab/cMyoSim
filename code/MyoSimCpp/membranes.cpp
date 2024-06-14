/**
* @file		m_state.cpp
* @brief	Source file for the m_state class
* @author	Ken Campbell
*/

#include <cstdio>

#include "membranes.h"
#include "half_sarcomere.h"

// Constructor
membranes::membranes(half_sarcomere* set_p_parent_hs)
{
	// Initialise
	p_parent_hs = set_p_parent_hs;

	memb_Ca_cytosol = 0;
}

// Destructor
membranes::~membranes(void)
{
	// Tidy up
}