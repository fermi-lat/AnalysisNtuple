/**
 *
 * @class UBinterpolate
 *
 * @brief interpolation class for the unbiased energy algorithm of Carmelo
 *
 * This class provides the interpolation values for Carmelo's implementation
 * of the unbiased energy algorithm.
 * The constructor reads in the interpolation table, the function
 * interpolate() returns an interpolated value for a given logE-zDir pair.
 *
 * @author Michael Kuss
 *
 * $Header: /nfs/slac/g/glast/ground/cvs/AnalysisNtuple/src/UBinterpolate.cxx,v 1.2 2012/01/17 16:52:00 cohen Exp $
 */

#include "UBinterpolate.h"

#include "facilities/Util.h"

#include <fstream>
#include <sstream>
#include <cstdlib>

UBinterpolate::UBinterpolate(std::string calibFileName) : m_calibFileName(calibFileName) {
    facilities::Util::expandEnvVar(&m_calibFileName);
    std::ifstream fin(m_calibFileName.c_str());
    if ( !fin ) {
        std::cerr << "Couldn't open " << m_calibFileName << std::endl;
        std::exit(1);
    }
    std::string s;
    int counter = 0;
    while ( !fin.eof() ) {
        std::getline(fin, s);
        if ( s.size() == 0 )
            continue;
        std::stringstream line(s);
        while ( !line.eof() ) {
            float f;
            line >> f;
            switch ( counter ) {
            case 0:
                m_x.push_back(f);
                break;
            case 1:
                m_y.push_back(f);
                break;
            default:
                m_z.push_back(f);
                break;
            }
        }
        ++counter;
    }
}

const float UBinterpolate::interpolate(const float logE, const float zDir) const {
    static const int OUTSIDE = 0;
    // checking if we are outside of the calibrated area
    if ( logE < getX(0) )
        return OUTSIDE;
    if ( logE > getX(m_x.size()-1) )
        return OUTSIDE;
    if ( zDir > getY(0) )
        return OUTSIDE;
    if ( zDir < getY(m_y.size()-1) )
        return OUTSIDE;

    // checking where we are
    unsigned int i = 0;
    unsigned int j = 0;
    for ( i=1; i<m_x.size(); ++i )
        if ( logE < getX(i) )
            break;
    for ( j=1; j<m_y.size(); ++j )
        if ( zDir > getY(j) )
            break;

    const float x   = logE;
    const float y   = zDir;
    const float x0  = getX(i-1);
    const float x1  = getX(i);
    const float y0  = getY(j-1);
    const float y1  = getY(j);
    const float z00 = getZ(i-1,j-1);
    const float z01 = getZ(i-1,j);
    const float z10 = getZ(i  ,j-1);
    const float z11 = getZ(i  ,j);

    // interpolating first in x
    const float z0 = ( z00 * ( x1 - x ) - z10 * ( x0 - x ) ) / ( x1 - x0 );
    const float z1 = ( z01 * ( x1 - x ) - z11 * ( x0 - x ) ) / ( x1 - x0 );

    // interpolating z0 and z1 in y
    const float z = ( z0 * ( y1 - y ) - z1 * ( y0 - y ) ) / ( y1 - y0 );

    /*
    std::cout << "we are " << i << ' ' << j << std::endl;
    std::cout << x0 << ' ' << y0 << ' ' << z00 << std::endl;
    std::cout << x1 << ' ' << y0 << ' ' << z10 << std::endl;
    std::cout << x0 << ' ' << y1 << ' ' << z01 << std::endl;
    std::cout << x1 << ' ' << y1 << ' ' << z11 << std::endl;
    std::cout << z0 << ' ' << z1 << ' ' << z << std::endl;
    */

    return z;
}
