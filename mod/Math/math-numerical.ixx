module;

#include<cmath>

export module math:Numerical;

import Config;

import :Variable;
import :Calculus;

namespace{

    constexpr int stop = 15;
    constexpr auto signal = []<jf::types::integral_or_floating_point_c Argtype>( const Argtype& arg)
                { 
                    return (arg<0? -1 : 1);
                };
} //  namespace


export namespace jf::math{

    template<jf::types::integral_or_floating_point_c Number, class FuncType>
    constexpr auto bissec(Number interval1, Number interval2, FuncType&& func, double precision = 0.000005)
    {
        
        if(func(interval1) * func(interval2) > 0){ throw("can´t calculate because: func(interval1) * func(interval2) > 0. "); }
        
        double result = 0.0;
        for(int count{}; count < stop; ++count)
        {
            result = (interval1 + interval2) / 2;
            if(std::abs(func(result)) < precision){ return result; }

            if(signal(func(result)) == signal(func(interval1))){ interval1 = result; }
            else{interval2 = result; } 
        }

        return result;

    }// bissec
     
    template<jf::types::integral_or_floating_point_c Number, class FuncType>
    constexpr auto false_position(Number interval1, Number interval2, FuncType&& func, double precision = 0.000005)
    {
       auto fun1 = func(interval1);
       auto fun2 = func(interval2);
        if(fun1 * fun2 > 0){ throw( "can´t calculate because: func(interval1) * func(interval2) > 0. "); }
        
        
        double result = 0.0;
        for(int count{}; count < stop; ++count)
        {
            result = (interval1*fun2 - interval2*fun1) / (fun2 - fun1);
            if(std::abs(func(result)) < precision){ return result; }

            if(signal(func(result)) == signal(func(interval1))){ interval1 = result; }
            else{interval2 = result; } 

            fun1 = func(interval1);
            fun2 = func(interval2);       
        }

        return result;

    }// false_position
    
    /*
     * take a function of type f(x) = 0 and  a x0 the first extimative 
     * @return the first repeate x decimal precision*/
    template<jf::types::integral_or_floating_point_c Number, class FuncType>
    constexpr auto fixed_point(FuncType&& func, Number x0, int order = 4)
    {
     
        double result = 0.0;
        double x1 = func(x0);

        double trunc1 = std::trunc(x1 * std::pow(10, order)) / std::pow(10, order);
        double trunc2 = 0.0;

        for(int count{}; count < stop; ++count)
        {
            result = func(x1);
            trunc2 = std::trunc(result * std::pow(10, order)) / std::pow(10, order);
            if(trunc1 == trunc2){ return result; }

            x1 = result;
            trunc1 = trunc2;

        }

        return result;

    }// fixed_point

    template<jf::types::integral_or_floating_point_c Number, class FuncType>
    constexpr auto newton_method(FuncType&& func, Number x0, double funcDelimiter = 0.000001, int order = 4)
    {

        double result = 0.0;
        double x1 = x0 - func(x0) / jf::math::five_point_stencil<1>(func, x0);

        double trunc1 = std::trunc(x1 * std::pow(10, order)) / std::pow(10, order);
        double trunc2 = 0.0;

        for(int count{}; count < stop; ++count)
        {
            result = x1 - func(x1) / jf::math::five_point_stencil<1>(func, x1);
            trunc2 = std::trunc(result * std::pow(10, order)) / std::pow(10, order);
            if(trunc1 == trunc2 || std::abs(func(result)) < funcDelimiter){ return result; }

            x1 = result;
            trunc1 = trunc2;

        }

        return result;
    } // newton_method

    template<jf::types::integral_or_floating_point_c Number, class FuncType>
    constexpr auto secant_method(FuncType&& func, Number x0, Number x1, double funcDelimiter = 0.000001)
    {

        auto fun0 = func(x0);
        auto fun1 = func(x1);
        double result = 0.0;


        for(int count{}; count < stop; ++count)
        {
            result = (x0 * fun1 - x1*fun0) / (fun1 - fun0);
            if(std::abs(func(result)) < funcDelimiter) return result;
            
            x0 = x1;
            x1 = result;
            fun0 = func(x0);
            fun1 = func(x1);

        }

        return result;
    } // secant_method
} // namespace jf::math


