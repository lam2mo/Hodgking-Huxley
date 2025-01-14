/*
    FILE NAME: HH-model.cpp
  DESCRIPTION: This code uses GSL to calculate the solution to the Hodgking Huxley Model

		[!] The system must be normalized to get the feeling of the biophysical behavior
	   	[!] See the txt files to more info about the model
      
       AUTHOR: Daniel Mejia Raigosa
       E-MAIL: danielmejia55@gmail.com
       GITHUB: http://github.com/Daniel-M/Hodgking-Huxley
         DATE: 18 November 2014
      VERSION: 3.0
    
    LICENSE
    
    This file is part of "Hodgking-Huxley".

    "HH-model.cpp" is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    "HH-model.cpp" is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with "HH-model.cpp".  If not, see <http://www.gnu.org/licenses/>.
    
    
*/

#include "HH_Includes.hpp"
#include "HH_Model_class.hpp"
//#include "mglGraphics.hpp"


// ==========================================================================
// 			FILE I/O HELPERS
// ==========================================================================

void getInfoFromFile(std::string filename, real buffer[])
{
    std::ifstream data;
    char sBuffer[20];
    int i=0;
    

    data.open (filename.c_str(), std::ifstream::in);
  
    while(data.eof()!=true)
    {
		data.getline(sBuffer,20);
		buffer[i]=atof(sBuffer);
		i++;
    }
   
       data.close();
}

void getInfoFromFile(std::string filename, std::vector<real> buffer)
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
// 			MAIN CODE
// ==========================================================================

int main(int argc, char **argv){

    cadna_init(-1, 0, 51, 4);

    int POINTS=20000;
    char verbosity='n',gengraphs='g', automatic='y';
    
    std::string paramfile = "parameters.txt";
    std::string iniconfile = "init_cond.txt";
    std::string basefile = "results-HH-model";
    std::string ext = "txt";
    
    std::string out_datafile = basefile + "." + ext;
    
    // ==========================================================================

    real h_step;  /// h_step This will be the h_step for integrate, will be redefined on line 333
	std::vector<real> parametros; /// parametros[8] Here one stores the ODES parameters for solvig, i.e membrane capacitance, conductances and so on
    
    real t=0.0;       /// t Starting time, this variable is a "time buffer"
    real tf=25.0;	/// tf 25 ms for time integration


// ==========================================================================
//   Console commands to run the programm
// ==========================================================================
    
    int c;
	extern int optind;
	extern char *optarg;
	
	while ((c = getopt(argc, argv, "hvgAt:p:k:i:b:e:")) != -1)
	  switch (c)  {
	  case 'h':
	    std::cout << "\nUsage: \n" << argv[0] << " [options] \n" << std::endl;
	    std::cout <<"\t-h  : This help message"<< std::endl;
	    std::cout <<"\t-v  : Set verbosity level,(No verbosity by default)"<< std::endl;
	    std::cout <<"\t-g  : Don't generate PNG graphics,(Graphics by default)"<< std::endl;
        std::cout <<"\t-A  : Automatic parameters and initial conditions (using built-in defaults)"<< std::endl;
        std::cout <<"\t-t# : Integrate from zero to # (In milliseconds),(Default= " << tf << ")" << std::endl;
	    std::cout <<"\t-p# : Set number of points used (Default " << POINTS << ")" << std::endl;
	    std::cout <<"\t-i@ : Initial conditions file [-i filename.ext] (Default \"" << iniconfile << "\")" << std::endl;
	    std::cout <<"\t-k@ : Parameters file [-kfilename.ext] (Default \"" << paramfile << "\")" << std::endl;
	    std::cout <<"\t-b@ : Set basename for output files [-b basename] (Default \"" << basefile << "\")" << std::endl;
	    std::cout <<"\t-e@ : Set extension for basename output files [-e ext = basename.ext] (Default \"" << ext << "\")" << std::endl;

	    std::cout << "\nNOTE: If no options are given, the -A flag is used by default\n" << std::endl;

	    std::cout << std::endl;

	    exit(1);
	    break;

	  case 'v':  verbosity='v';
	    break;
	    
	  case 'g':  gengraphs='n';
	    break;
        
	  case 'A':  automatic='y'; // Automatic mode
        break;

	  case 't': tf=atof(optarg);
	    break;

	  case 'p': POINTS = atoi(optarg);
	    break;
	   
	  case 'i': iniconfile.assign(optarg);
		automatic='n'; //Disable automode 
	    break;
	    
	  case 'k': paramfile.assign(optarg);
		automatic='n'; //Disable automode 

	    break;
	    
	  case 'b': basefile.assign(optarg);
	    break;
	    
	  case 'e': ext.assign(optarg);
	    break;

	  case '?':
           std::cout << "Unknown option " << std::endl;
	       std::cout << "Aborting..." << std::endl;
             
         default:
            abort ();
	    
	  }
	  
	// Here the solution is stored
	// definde here in order to have a container for the initial conditions
	std::vector<real> y; 

	if (automatic=='y')
    {
               
// ==========================================================================
// Initial conditions AUTOMATED
// ==========================================================================

			y.push_back(-60);   // Initial transmembrane potential, assuming resting potential (-60 mV) *when automated* [-A flag]
		    y.push_back(0);     // Initial state for Gating function m *when automated* [-A flag]
		    y.push_back(0);  // Initial state for Gating function h *when automated* [-A flag]
		    y.push_back(0);  // Initial state for Gating function n *when automated* [-A flag]
 
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
    
			parametros.push_back(0.01);   // Membrane capacitance  *when automated* [-A flag]
		    parametros.push_back(0);    // induced current on axon, 0 means no external current *when automated* [-A flag]
		    parametros.push_back(1.2);   // Na conductances *when automated* [-A flag]
		    parametros.push_back(55.17);  // Na Nernst Potential *when automated* [-A flag]
		    parametros.push_back(0.36);   // K Conductance *when automated* [-A flag]
		    parametros.push_back(-72.14); // K Nernst Potential *when automated* [-A flag]
		    parametros.push_back(0.003);  // Leakage conductance (Due to a Cl current) *when automated* [-A flag]
		    parametros.push_back(-49.42); // Leakage Nernst potential (Due to a Cl current) *when automated* [-A flag]
 
	        
    }
    else
    {
// ==========================================================================
// Initial conditions FROM FILE
// ==========================================================================        
        
            getInfoFromFile(iniconfile,y);
            getInfoFromFile(paramfile,parametros);
    }
        
        
    std::cout << "Initial Conditions" << std::endl;
    
    std::cout << y[0] << std::endl;
    std::cout << y[1] << std::endl;
    std::cout << y[2] << std::endl;
    std::cout << y[3] << std::endl;
    
    std::cout << "Parameters" << std::endl;
    
    std::cout << parametros[0] << std::endl;
    std::cout << parametros[1] << std::endl;
    std::cout << parametros[2] << std::endl;
    std::cout << parametros[3] << std::endl;      
    std::cout << parametros[4] << std::endl;
    std::cout << parametros[5] << std::endl;
    std::cout << parametros[6] << std::endl;
    std::cout << parametros[7] << std::endl;
		    
	    
// ==========================================================================  
// ==========================================================================
//  Now the iterative process of solution using ODEINT by Libboost
// ==========================================================================
// ==========================================================================

	hh_model hhmodel(parametros);

    hhSolver(POINTS, tf, hhmodel, y, out_datafile);
 
	if(gengraphs!='n')
	{
		//mglDrawFunction(parametros,POINTS,out_datafile);
	}

    cadna_end();

	exit(0);
}//End of Main Code


