/*
    LICENSE

    This file is an implementation of "Hodgkin-Huxley" using Boost ODEint.

    "hh-boost.cpp" is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    "hh-boost.cpp" is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with "hh-boost.cpp".  If not, see <http://www.gnu.org/licenses/>.
*/

#include <cstdlib>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include <getopt.h>

#include <boost/numeric/odeint.hpp>
using namespace boost::numeric::odeint;

#ifdef USE_CADNA
    #include <cadna.h>
    typedef double_st real_t;
#else
    typedef double real_t;
#endif


// ==========================================================================
//          DATA STRUCTURES
// ==========================================================================

/**
 * @brief Primary Hodgkin-Huxley neuron model
 *
 * The default constructor receives a std::vector<real_t> that contains the
 * parameters of the model in the following order:
 *   * Membrane capacitance
 *   * Induced current on axon, 0 means no external current
 *   * Na conductances
 *   * Na Nernst Potential
 *   * K Conductance
 *   * K Nernst Potential
 *   * Leakage conductance (Due to a Cl current)
 *   * Leakage Nernst potential (Due to a Cl current)
 *
 */

class hh_model
{
	private:
		std::vector<real_t> parameters;

	public:

		hh_model(std::vector<real_t> params) : parameters(params) { }
		void operator()(const std::vector<real_t> &y, std::vector<real_t> &f, const real_t /* t */);
};

/**
 * @brief Neuron state vector
 */
typedef std::vector<real_t> state_t;

/**
 * @brief Structure that stores each set of points during the integration step.
 *
 * This structure stores the time value and the variable set in one file, this is
 * the libboost odeint observer so the vectors x and t contain each time step value for time,
 * voltage, and gating probabilities h, m, and n.
 */
struct push_back_state_and_time
{
    std::vector< std::vector<real_t> >& m_states;
    std::vector< real_t >& m_times;

    push_back_state_and_time(std::vector< std::vector<real_t> > &states, std::vector<real_t> &times)
	: m_states(states), m_times(times) { }

    void operator()( const std::vector<real_t> &x, real_t t)
    {
    	m_states.push_back(x);
    	m_times.push_back(t);
    }
};


// ==========================================================================
//          ODE SOLVER
// ==========================================================================

/* ===============================================================================
 * These are taken from book "Nonlinear dynamics in physiology and medicine"
 * From Anne Beuter et al.
 * Chapter 4 - Excitable cells, pp 104
 ===============================================================================*/

real_t alpha_n(real_t V)
{
    return (0.01*(V+50.0))/(1.0-exp(-(V+50.0)/10.0));
}

real_t beta_n(real_t V)
{
    return 0.125*exp(-(V+60.0)/80.0);
}

real_t alpha_m(real_t V)
{
    return (0.1*(V+35.0))/(1-exp(-(V+35.0)/10.0));
}

real_t beta_m(real_t V)
{
    return 4.0*exp(-(V+60.0)/18.0);
}

real_t alpha_h(real_t V)
{
    return 0.07*exp(-(V+60.0)/20.0);
}

real_t beta_h(real_t V)
{
    return 1.0/(1.0+exp(-(V+30.0)/10.0));
}

/**
 * @brief Main ODE equations
 */
void hh_model::operator()(const std::vector<real_t> &y, std::vector<real_t> &f, const real_t /*t*/)
{
    f[0]=(1/parameters[0])*(parameters[1]-parameters[2]*pow(y[1],3)*y[2]*(y[0]-parameters[3])-parameters[4]*pow(y[3],4)*(y[0]-parameters[5])-parameters[6]*(y[0]-parameters[7]));
    f[1]=alpha_m(y[0])*(1-y[1])-beta_m(y[0])*y[1];
    f[2]=alpha_h(y[0])*(1-y[2])-beta_h(y[0])*y[2];
    f[3]=alpha_n(y[0])*(1-y[3])-beta_n(y[0])*y[3];
}

/**
 * @brief ODE solution driver
 *
 * This function performs the numerical solution of the hh_model object taking
 * the initial values stored in y and writing output to the file out_datafile.
 *
 * @param points        Number of points desired to make the numerical integral.
 * @param tf	        Final time (in ms) of the integration.
 * @param hhmodel	    hh_model object that contains the ODES and gating functions.
 * @param y			    Vector for the solution (and initial conditions)
 * @param out_datafile	The filename where the solution is going to be placed
 * @param all_output    Include full output as well as potential
 */
void hhSolver(int points, real_t tf, hh_model hhmodel, std::vector<real_t> y,
        std::string out_datafile, char all_output)
{
	std::vector<std::vector<real_t>> y_vec;
	std::vector<real_t> times;

	std::ofstream hh_action_potential(out_datafile.c_str());
	std::ofstream hh_data((out_datafile+"-all").c_str());

	real_t h_step=tf/points;

    // Original:
    //   size_t steps = integrate(hhmodel,y,0.0,tf,h_step,push_back_state_and_time( y_vec , times ));

    // Full list of steppers:
    //   https://www.boost.org/doc/libs/1_66_0/libs/numeric/odeint/doc/html/boost_numeric_odeint/getting_started/overview.html
    //
    euler<state_t> stepper;
    //runge_kutta4<state_t> stepper;
    //runge_kutta_dopri5<state_t> stepper;
    size_t steps = integrate_adaptive(stepper, hhmodel, y, real_t(0.0), tf, h_step,
            push_back_state_and_time(y_vec, times));

	// output
	for( size_t i=0; i<=steps; i++ )
	{
		//std::cout << times[i] << '\t' << y_vec[i][0] << std::endl;
		hh_action_potential << times[i] << '\t' << y_vec[i][0]  << std::endl;

        if (all_output == 'y') {
            hh_data << times[i] << "\t";
            for(int j=0;j<4;j++) {
                hh_data	<< y_vec[i][j] << "\t";
            }
            hh_data << std::endl;
        }
	}
}


// ==========================================================================
//          UTILITY FUNCTIONS
// ==========================================================================

void getInfoFromFile(std::string filename, std::vector<real_t> buffer)
{
    std::ifstream data;
    char sBuffer[20];
    int i=0;

    data.open (filename.c_str(), std::ifstream::in);
    while(data.eof()!=true)
    {
		data.getline(sBuffer,20);
		buffer.push_back(atof(sBuffer));
		i++;
    }
    data.close();
}

// ==========================================================================
//          ENTRY POINT
// ==========================================================================
int main(int argc, char **argv)
{
    // program options
    char all_output='n', automatic='y';
    std::string paramfile = "parameters.txt";
    std::string iniconfile = "init_cond.txt";
    std::string basefile = "results";
    std::string ext = "txt";
    std::string out_datafile = basefile + "." + ext;

    // simulation state and options
    std::vector<real_t> y;
    std::vector<real_t> parameters;
    int points = 500;
    real_t tf  = 25.0;   // 25 ms for time integration

    // parse command line options
    int c;
    extern char *optarg;
    while ((c = getopt(argc, argv, "hvgAt:p:k:i:b:e:")) != -1)
        switch (c)  {

        case 'h':
            std::cout << "\nUsage: \n" << argv[0] << " [options] \n" << std::endl;
            std::cout <<"\t-h  : This help message"<< std::endl;
            std::cout <<"\t-a  : Enable all output (off by default)" << std::endl;
            std::cout <<"\t-A  : Automatic parameters and initial conditions (using built-in defaults)"<< std::endl;
            std::cout <<"\t-t# : Integrate from zero to # in milliseconds (Default= " << tf << ")" << std::endl;
            std::cout <<"\t-p# : Set number of points used (Default " << points << ")" << std::endl;
            std::cout <<"\t-i@ : Initial conditions file [-i filename.ext] (Default \"" << iniconfile << "\")" << std::endl;
            std::cout <<"\t-k@ : Parameters file [-kfilename.ext] (Default \"" << paramfile << "\")" << std::endl;
            std::cout <<"\t-b@ : Set basename for output files [-b basename] (Default \"" << basefile << "\")" << std::endl;
            std::cout <<"\t-e@ : Set extension for basename output files [-e ext = basename.ext] (Default \"" << ext << "\")" << std::endl;
            std::cout << "\nNOTE: If no options are given, the -A flag is used by default\n" << std::endl;
            std::cout << std::endl;
            exit(1);
            break;

        case 'a': all_output='y';                           break;
        case 'A':                            automatic='y'; break;
        case 't': tf=atof(optarg);                          break;
        case 'p': points = atoi(optarg);                    break;
        case 'i': iniconfile.assign(optarg); automatic='n'; break;
        case 'k': paramfile.assign(optarg);  automatic='n'; break;
        case 'b': basefile.assign(optarg);                  break;
        case 'e': ext.assign(optarg);                       break;

        case '?':
            std::cout << "Unknown option " << std::endl;
            std::cout << "Aborting..." << std::endl;
        default:
            abort();
      }

#ifdef USE_CADNA
    cadna_init(-1, 0, 51, 4);
#endif

    // set initial conditions and parameters
    if (automatic=='y') {

        // initial conditions AUTOMATED

        y.push_back(-60);   // Initial transmembrane potential, assuming resting potential (-60 mV)
        y.push_back(0);     // Initial state for Gating function m
        y.push_back(0);     // Initial state for Gating function h
        y.push_back(0);     // Initial state for Gating function n

        // # Membrane capacitance in Farads
        // Cm = 0.001;
        //
        // # Ion conductances in mMho
        // gNa = 120;
        // gK  = 36;
        // gL  = 0.3;
        //
        // # Ion equilibrium potentials in mVolts
        // vNa = -115;
        // vL  = 10.613;
        // vK  = 12;

        parameters.push_back(0.01);   // Membrane capacitance
        parameters.push_back(0);      // induced current on axon, 0 means no external current
        parameters.push_back(1.2);    // Na conductances
        parameters.push_back(55.17);  // Na Nernst Potential
        parameters.push_back(0.36);   // K Conductance
        parameters.push_back(-72.14); // K Nernst Potential
        parameters.push_back(0.003);  // Leakage conductance (Due to a Cl current)
        parameters.push_back(-49.42); // Leakage Nernst potential (Due to a Cl current)

    } else {

        // initial conditions FROM FILE
        getInfoFromFile(iniconfile, y);
        getInfoFromFile(paramfile, parameters);
    }

    // print information
    std::cout << "Initial Conditions" << std::endl;
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y[2] << std::endl;
    std::cout << y[3] << std::endl;
    std::cout << "Parameters" << std::endl;
    std::cout << parameters[0] << std::endl;
    std::cout << parameters[1] << std::endl;
    std::cout << parameters[2] << std::endl;
    std::cout << parameters[3] << std::endl;
    std::cout << parameters[4] << std::endl;
    std::cout << parameters[5] << std::endl;
    std::cout << parameters[6] << std::endl;
    std::cout << parameters[7] << std::endl;

    // iterative solution
    hh_model hhmodel(parameters);
    hhSolver(points, tf, hhmodel, y, out_datafile, all_output);

#ifdef USE_CADNA
    cadna_end();
#endif

    return EXIT_SUCCESS;
}


