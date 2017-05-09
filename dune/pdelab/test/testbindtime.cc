#ifdef HAVE_CONFIG_H
#include "config.h"
#endif
#include <iostream>
#include <dune/common/classname.hh>
#include <dune/pdelab/function/bindtime.hh>

int main()
{
    auto f = [](double t, int i)
    {
        return i*t;
    };
    auto f2 = [](int i, double t)
    {
        return i*t;
    };

    auto f_t = Dune::PDELab::bindTime(f,Dune::Indices::_1);
    auto f2_t = Dune::PDELab::bindTime(f2,Dune::Indices::_2);

    f_t.setTime(1.5);
    f2_t.setTime(0.75);

    std::cout << f(1.5,2) << std::endl;
    std::cout << f2(4,0.75) << std::endl;
    std::cout << f_t(2) << std::endl;
    std::cout << f2_t(4) << std::endl;
}
