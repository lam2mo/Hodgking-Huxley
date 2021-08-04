/*
    FILE NAME: HH_Model_class.hpp
  DESCRIPTION: Part of Hodgking-Huxley model, used in HH-model.cpp

       AUTHOR: Daniel Mejia Raigosa
       E-MAIL: danielmejia55@gmail.com
       GITHUB: http://github.com/Daniel-M/Hodgking-Huxley
         DATE: 17 November 2014
      VERSION: 1.0
    
    LICENSE
    
    This file is part of "Hodgking-Huxley".

    "HH_Model_ODES.hpp" is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    "HH_Model_ODES.hpp" is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with "HH_Model_ODES.hpp".  If not, see <http://www.gnu.org/licenses/>.
    
    
*/

#include "HH_Includes.hpp"


/*
 * Gating functions
 */

real alpha_n(real V);
real beta_n(real V);
real alpha_m(real V);
real beta_m(real V);
real alpha_h(real V);
real beta_h(real V);


/* ==========================================================================
/     Defining the ODE system for the Hodgking Huxley Model
/ ==========================================================================*/
//\file HH_Model_ODES.hpp 
/*!
 * \brief Defines the ODE system proposed by Hodgking and Huxley as a class
 *
 * The default constructor receives a std::vector<real> that contains the parameters of the model in the order
 *   * Membrane capacitance 
 *   * Induced current on axon, 0 means no external current 
 *   * Na conductances
 *   * Na Nernst Potential
 *   * K Conductance
 *   * K Nernst Potential
 *   * Leakage conductance (Due to a Cl current)
 *   * Leakage Nernst potential (Due to a Cl current)
 *
 * So parameter[0] contains membrane capacitance and so on
 */
  
class hh_model
{
	private:
	
		std::vector<real> parameters;

	public:

		hh_model( std::vector<real> params ) : parameters(params) { }

        void operator()( const std::vector<real>&, std::vector<real>&, const real );
};


void hhSolver(int POINTS, real tf, hh_model hhmodel,std::vector<real> y,std::string out_datafile);

