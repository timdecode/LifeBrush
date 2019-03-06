//============================================================
// Problem.h
//
// A Problem reads in a problem description from the command
// line and computes the corresponding solution (i.e., a
// tangent direction field on the specified mesh with the
// specified singularities).  Input is given by a regular
// ASCII text file containing lines of the following format
// in any order:
//
//    in [path to input mesh file]
//    out [path to output data]
//    vertex [0-based vertex ID] [target holonomy]
//    generator [0-based generator ID] [target holonomy]
//    angle [initial field angle]
//
// Terms in square brackets [] need to be specified by
// the user.  An example input file can be found in
// test/problem.txt
//

#ifndef PROBLEM_H
#define PROBLEM_H

#include <vector>
#include <iosfwd>
#include <string>

namespace tcods
{
    typedef std::pair<int,double> Singularity;
    typedef std::pair<int,double> Generator;
    
    class Problem
    {
    public:
        Problem( void );               // default constructor
        Problem( std::istream& in );   // construct a problem from a valid istream
        void read( std::istream& in ); // load a problem from a valid istream
        void solve( void ) const;      // solve the problem
        
    protected:
        std::string inputPath;
        std::string outputPath;
        
        double fieldAngle;
        
        std::vector<Singularity> singularities;
        
        std::vector<Generator> generators;
    };
}

#endif

