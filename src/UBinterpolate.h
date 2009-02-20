/**
 * @class UBinterpolate
 *
 * @brief interpolation class for the unbiased energy algorithm of Carmelo
 *
 * @author Michael Kuss
 *
 * $Header$
 */

#ifndef __UBINTERPOLATE_H__
#define __UBINTERPOLATE_H__

#include <iostream>
#include <string>
#include <vector>

class UBinterpolate {
 public:
    UBinterpolate(std::string calibFileName);
    ~UBinterpolate() {}

    /// returns an interpolated value for zDir and LogE
    const float interpolate(const float logE, const float zDir) const;

    /// print the interpolation table (for debugging)
    void print() const {
        std::cout << "calib file name: " << m_calibFileName << std::endl;
        std::cout << "member vector sizes: " << m_x.size() << ' ' << m_y.size() << ' ' << m_z.size() << std::endl;
        std::cout.precision(7);
        for ( unsigned int j=0; j<m_y.size(); ++j )
            for ( unsigned int i=0; i<m_x.size(); ++i )
                std::cout << getX(i) << ' ' << getY(j) << ' ' << getZ(i,j) << std::endl;
    }

 private:
    std::string m_calibFileName;
    /// vector of LogE
    std::vector<float> m_x;
    /// vector of ZDir
    std::vector<float> m_y;
    /// matrix of values
    std::vector<float> m_z;

    const float getX(const int i) const { return m_x[i]; }
    const float getY(const int i) const { return m_y[i]; }
    /// returns a value from the interpolation table
    const float getZ(const int i, const int j) const { return m_z[i+m_x.size()*j]; }
};

#endif
