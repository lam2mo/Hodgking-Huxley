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

/* ==========================================================================
/     Defining the ODE system for the Hodgking Huxley Model
/ ==========================================================================*/
//\file HH_Model_ODES.hpp 

#include "HH_Includes.hpp"
#include "HH_Model_class.hpp"

using namespace boost::numeric::odeint;

/* ===============================================================================
/ Defining some auxiliary "gating functions" for the Hodgking Huxley Model
/ These are taken from book "Mathematical aspects of the Hodgking-Huxley
/ Neural theory"
/ From Janne Cronin
/ Chapter 2.3 - The work of Hodgking-Huxley, pp 50-52
/ ===============================================================================*/

/* Jim Sochacki
   The gating functions from HH 1952 paper */

/*! \def alpha_n
 * \brief The alpha gating function for the n gate
 *
 */
real alpha_n(real V)
{
     return (0.01*(10-V))/(exp((10-V)/10.0)-1);
}

/*! \def beta_n
 * \brief The beta gating function for the n gate
 *
 */
real beta_n(real V)
{
     return 0.125*exp(-V/80.0);
}

/*! \def alpha_m
 * \brief The alpha gating function for the m gate
 *
 */
real alpha_m(real V)
{
     return (0.1*(25.0-V))/(1-exp((25.0-V)/10.0));
}

/*! \def beta_m
 * \brief The beta gating function for the m gate
 *
 */
real beta_m(real V)
{
     return 4.0*exp(-V/18.0);
}

/*! \def alpha_h
 * \brief The alpha gating function for the h gate
 *
 */
real alpha_h(real V)
{
     return 0.07*exp(-V/20.0);
}

/*! \def beta_h
 * \brief The beta gating function for the h gate
 *
 */
real beta_h(real V)
{
     return 1.0/(1.0+exp((30.0-V)/10.0));
}

void hh_model::operator()( const std::vector<real> &y , std::vector<real> &f , const real /* t */ )
{
    f[0]=(1/parameters[0])*(parameters[1]-parameters[2]*pow(y[1],3)*y[2]*(y[0]-parameters[3])-parameters[4]*pow(y[3],4)*(y[0]-parameters[5])-parameters[6]*(y[0]-parameters[7]));
    f[1]=alpha_m(y[0])*(1-y[1])-beta_m(y[0])*y[1];
    f[2]=alpha_h(y[0])*(1-y[2])-beta_h(y[0])*y[2];
    f[3]=alpha_n(y[0])*(1-y[3])-beta_n(y[0])*y[3];
}

/*! \brief Structure that estores each set of points during the integration step.
 *
 * This structre stores the time value and the variable set in one file, this is
 * the libboost odeint observer so the vectors x and t contain each time step value for time,
 * voltaje, and gating probabilities h, m, and n.
 *
 * */ 

struct push_back_state_and_time
{
    std::vector< std::vector<real> >& m_states;
    std::vector< real >& m_times;

    push_back_state_and_time( std::vector< std::vector<real> > &states , std::vector< real > &times )
	: m_states( states ) , m_times( times ) { }

    void operator()( const std::vector<real> &x , real t )
    {
    	m_states.push_back( x );
    	m_times.push_back( t );
    }
};


/*! \brief The ODE Solver for the HH Model
 *
 * This function performs the numerical solition of the hh_model object taking the values stored on y 
 * and placing them on the file out_datafile
 * \param POINTS Number of points desired to make the numerical integral.
 * \param tf	Final time (in ms) of the integration.
 * \param hhmodel	hh_model object that contains the ODES and gating functions.
 * \param y			a vector that stores the solution (and initial conditions)
 * \param out_datafile	The filename where the solution is going to be placed
 */

void hhSolver(int POINTS, real tf, hh_model hhmodel,std::vector<real> y,std::string out_datafile)
{
	std::vector<std::vector<real>> y_vec;
	std::vector<real> times;

	std::ofstream hh_action_potential(out_datafile.c_str());
	std::ofstream hh_data((out_datafile+"-all").c_str());

	//hh_model hhmodel(parametros);

	real h_step=tf/POINTS;

    // Original:
    //   size_t steps = integrate(hhmodel,y,0.0,tf,h_step,push_back_state_and_time( y_vec , times ));

    // Full list of steppers:
    //   https://www.boost.org/doc/libs/1_66_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html
    //
    typedef std::vector<real> state_t;
    //euler<state_t> stepper;
    runge_kutta4<state_t> stepper;
    //runge_kutta_dopri5<state_t> stepper;
	size_t steps = integrate_adaptive(stepper,hhmodel,y,real(0.0),tf,h_step,push_back_state_and_time( y_vec , times ));
  
	/* output */
	for( size_t i=0; i<=steps; i++ )
	{
		//std::cout << times[i] << '\t' << y_vec[i][0] << std::endl; 
		hh_action_potential << times[i] << '\t' << y_vec[i][0]  << std::endl;
		
		hh_data << times[i] << "\t";
		
		for(int j=0;j<4;j++)
		{
			hh_data	<< y_vec[i][j] << "\t";
		}
	
		hh_data << std::endl;
	}
}
