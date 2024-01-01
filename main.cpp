#include <eigen3/Eigen/Dense>
#include <iostream>

int main(int argc, char const *argv[])
{
    int N = 1;
    if (argc > 1)
    {
        N = std::stoi(argv[1]);
    }
    
    Eigen::VectorXd cluster = Eigen::VectorXd::Zero(N);
    
    return 0;
}
