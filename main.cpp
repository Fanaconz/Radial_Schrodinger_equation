#include "boundSolver.h"

using namespace std;

int main(int argc, char** argv) {
    
    std::function<double(double)> v = \
        [](double r) -> double {return -10/r;};
    
    boundSolver bs(v, 1, 10.0);
    bs.writeWF("energ.txt");
    bs.printE();
    return 0;
}
