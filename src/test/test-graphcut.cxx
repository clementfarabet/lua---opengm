#include <arraypow2/arraypow2.hxx>
#include <opengm/inference/graphcut.hxx>

class GraphCutTest {
public:
    typedef float Energy;
    typedef opengm::DiscreteSpace Space;
    typedef opengm::ExplicitFactor<Energy, ArrayPow2<Energy> > Factor;
    typedef opengm::GraphCut<Factor> GraphCut;
    typedef GraphCut::gm_type GraphicalModel;

    GraphCutTest() : gridSize_(15) {
        // build the discrete space
        std::vector<size_t> numbersOfStates(gridSize_*gridSize_, 2);
        space_ = opengm::DiscreteSpace(numbersOfStates.begin(), numbersOfStates.end());

        // single site factors
        for(size_t j=0; j<space_.dimension(); ++j) {
            size_t vi[] = {j};
            Factor factor(space_, vi, vi+1);
            factor(0) = static_cast<Energy>(rand());
            factor(1) = static_cast<Energy>(RAND_MAX) - factor(0);
            gm_.addFactor(factor);
        }

        // Potts factors
        float penalty = static_cast<Energy>(RAND_MAX) / 5.0f;
        for(size_t y0=0; y0<gridSize_; ++y0) {
            for(size_t x0=0; x0<gridSize_; ++x0) {
                if(x0 != gridSize_-1) {
                    size_t x1 = x0+1;
                    size_t y1 = y0;
                    size_t q0 = x0 + gridSize_*y0;
                    size_t q1 = x1 + gridSize_*y1;
                    size_t vi[] = {q0, q1};
                    Factor factor(space_, vi, vi+2);
                    factor(0,0) = 0.0f;
                    factor(0,1) = penalty;
                    factor(1,0) = penalty;
                    factor(1,1) = 0.0f;
                    gm_.addFactor(factor);
                }
                if(y0 != gridSize_-1) {
                    size_t x1 = x0;
                    size_t y1 = y0+1;
                    size_t q0 = x0 + gridSize_*y0;
                    size_t q1 = x1 + gridSize_*y1;
                    size_t vi[] = {q0, q1};
                    Factor factor(space_, vi, vi+2);
                    factor(0,0) = 0.0f;
                    factor(0,1) = penalty;
                    factor(1,0) = penalty;
                    factor(1,1) = 0.0f;
                    gm_.addFactor(factor);
                }
            }
        }
    }

    template<class StateIterator>
    void printState(StateIterator it) {
        float energy = gm_.evaluate(it);
        std::cout << "E = " << energy << std::endl;
        for(size_t y=0; y<gridSize_; ++y) {
            for(size_t x=0; x<gridSize_; ++x) {
                std::cout << *it << ' ';
                ++it;
            }
            std::cout << std::endl;
        }
    }

    void run() {
        GraphCut gc(gm_);

        std::vector<size_t> state(gm_.space().dimension());

        // cut by push relabel max flow
        gc.setMaxFlowAlgorithm(GraphCut::PUSH_RELABEL);
        gc.infer();
        gc.arg(state);
        // printState(state.begin());

        // cut by Edmonds Karp max flow
        gc.setMaxFlowAlgorithm(GraphCut::EDMONDS_KARP);
        gc.infer();
        gc.arg(state);
        // printState(state.begin());

        // cut by Kolmogorov max flow
        gc.setMaxFlowAlgorithm(GraphCut::KOLMOGOROV);
        gc.infer();
        gc.arg(state);
        // printState(state.begin());
    }

private:
    size_t gridSize_;
    Space space_;
    GraphicalModel gm_;
};

int main() {
    std::cout << "GraphCut Test... ";
       { GraphCutTest t; t.run(); }	
    std::cout << "passed." << std::endl;

    return 0;
}
