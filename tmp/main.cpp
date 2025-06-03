#include <iostream>
#include <armadillo>
#include "utils.h"
#include "peltL2.h"

int main()
{
    try
    {
        const arma::mat signal = loadSignal("/Users/charles/Documents/GitHub/rupturesRcpp/tmp/data/signal.csv");
        std::cout << "Matrix dimensions: " << signal.n_rows << " x " << signal.n_cols << std::endl;
        const auto bkps = peltL2(signal, 30000, 1, 1);

        for (const int& elem : bkps)
        {
            std::cout << elem << " ";
        }
        std::cout << std::endl;
        std::cout << "N BKPS: " << bkps.size() << std::endl;

    }
    catch (const std::exception& e)
    {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}