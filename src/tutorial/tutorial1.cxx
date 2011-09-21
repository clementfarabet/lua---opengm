#include <marray/marray_hdf5.hxx>

#include <opengm/discretespace.hxx>
#include <opengm/explicitfactor.hxx>
#include <opengm/graphicalmodel.hxx>
#include <opengm/serialization.hxx>
#include <opengm/adder.hxx>
#include <opengm/minimizer.hxx>
#include <opengm/decomposer.hxx>

int main() 
{
    typedef double Energy;
    typedef opengm::DiscreteSpace Space;
    typedef opengm::ExplicitFactor<Energy> Factor;
    typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;

    size_t gridSize = 10;
    Energy alpha = 10000;
    std::string filename = "gm.h5";

    // discrete state space and graphical model
    std::vector<size_t> numbersOfStates(gridSize*gridSize, 2);
    Space space(numbersOfStates.begin(), numbersOfStates.end());
    GraphicalModel gm;

    // single site factors
    for(size_t j=0; j<space.dimension(); ++j) {
        Factor factor(space, &j, &j+1);
        factor(0) = static_cast<Energy>(rand());
        factor(1) = static_cast<Energy>(RAND_MAX) - factor(0);
        gm.addFactor(factor);
    }

    // Potts factors
    for(size_t y0=0; y0<gridSize; ++y0) {
        for(size_t x0=0; x0<gridSize; ++x0) {
            if(x0 != gridSize-1) {
                size_t x1 = x0+1;
                size_t y1 = y0;
                size_t vi[] = {x0 + gridSize*y0, x1 + gridSize*y1};
                Factor factor(space, vi, vi+2);
                factor(0,0) = 0.0;
                factor(0,1) = alpha;
                factor(1,0) = alpha;
                factor(1,1) = 0.0;
                gm.addFactor(factor);
            }
            if(y0 != gridSize-1) {
                size_t x1 = x0;
                size_t y1 = y0+1;
                size_t vi[] = {x0 + gridSize*y0, x1 + gridSize*y1};
                Factor factor(space, vi, vi+2);
                factor(0,0) = 0.0;
                factor(0,1) = alpha;
                factor(1,0) = alpha;
                factor(1,1) = 0.0;
                gm.addFactor(factor);
            }
        }
    }

    // serialize graphical model
    marray::Vector<double> serialization;
    opengm::serialize(gm, serialization);

    // save graphical model
    hid_t file = marray::hdf5::createFile(filename);
    marray::hdf5::save(file, "graphical-model", serialization);
    marray::hdf5::closeFile(file);

    return 0;
}
