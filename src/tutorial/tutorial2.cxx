#include <iostream>
#include <marray/marray_hdf5.hxx>

int main() 
{
    std::string filename = "result.h5";
    try {
        hid_t file = marray::hdf5::openFile(filename);
        marray::Vector<double> runtimes;
        marray::hdf5::load(file, "runtimes", runtimes);
        marray::Vector<double> values;
        marray::hdf5::load(file, "values", values);
        marray::Vector<size_t> state;
        marray::hdf5::load(file, "state", state);
        marray::hdf5::closeFile(file);

        std::cout << "E=" << values[values.size()-1]
            << " at iteration " << values.size()
            << ", after t=" << runtimes[runtimes.size()-1] / 1000
            << "s" << std::endl
            << "Approximate optimal state:" << std::endl;

        for(size_t y=0; y<10; ++y) {
            for(size_t x=0; x<10; ++x) {
                std::cout << state[x + y*10];
            }
            std::cout << std::endl;
        }

    }
    catch(std::runtime_error& e) {
        std::cout << "error: " << e.what() << std::endl;
    }
    return 0;
}
