 /*
 * Created on Jul 12, 2016
 * Author: Wei Chen
 */

#pragma once
#ifndef FREE_ENERGY_H
#define	FREE_ENERGY_H

template<typename T>
class GrandPotential : public cppoptlib::Problem<T> {
  public:
    using typename cppoptlib::Problem<T>::TVector;
    // this is just the objective (NOT optional)
    T value(const TVector &density);
    void gradient(const TVector &density, TVector &grad);
}

#endif
