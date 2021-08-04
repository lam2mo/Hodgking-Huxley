/*
    FILE NAME: HH_Includes.hpp
  DESCRIPTION: This code uses GSL to calculate the solution to the Hodgking Huxley Model

        [!] The system must be normalized to get the feeling of the biophysical behavior
        [!] See the txt files to more info about the model
      
       AUTHOR: Daniel Mejia Raigosa
       E-MAIL: danielmejia55@gmail.com
       GITHUB: http://github.com/Daniel-M/Hodgking-Huxley
         DATE: 13 Sept 2012
      VERSION: 1.0
    
    LICENSE
    
    This file is part of "Hodgking-Huxley".

    "HH_Includes.hpp" is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    "HH_Includes.hpp" is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with "HH_Includes.hpp".  If not, see <http://www.gnu.org/licenses/>.
    
    
*/

// Program Headers

#ifndef HH_INCLUDES
#define HH_INCLUDES

#include <iostream>
#include <string>
#include <iomanip>
#include <fstream>
#include <cmath>

// C old headers
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

// ODEINT Headers

#include <vector>
#include <boost/numeric/odeint.hpp>

#include <cadna.h>
typedef double_st real;

#endif

