module;

#include<cmath>

export module math:Numerical;

import Config;

import :Variable;
import :Calculus;

namespace types = jf::types;

namespace{

    constexpr int stop = 15;
    constexpr auto signal = []<types::integral_or_floating_point_c Argtype>( const Argtype& arg)
                { 
                    return (arg<0? -1 : 1);
                };
} //  namespace


export namespace jf::math{

    template<types::integral_or_floating_point_c Number, class FuncType>
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
     
    template<types::integral_or_floating_point_c Number, class FuncType>
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
    template<types::integral_or_floating_point_c Number, class FuncType>
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

    template<types::integral_or_floating_point_c Number, class FuncType>
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

    template<types::integral_or_floating_point_c Number, class FuncType>
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
      
     
      
    //// equations to make zero the terms under pricipal diagonal
    /// under a00:
    /// L1 = L1 - a10/a00 * L0 ==> L0 = line0 
    /// L2 = L2 - a20/a00 * L0
    /// under a11:
    /// L2 = L2 - a21/a11 * L1
    auto gauss_method( auto&& f1, auto&& f2, auto&& f3, auto&& point){
        // terms multipling x, y and z
        auto [x0, y0, z0] = types::process_func_arg_pairs(f1, std::tuple{1, 0, 0}, f1, std::tuple{0, 1, 0}, f1, std::tuple{0, 0, 1});

        auto [x1, y1, z1] = types::process_func_arg_pairs(f2, std::tuple{1, 0, 0}, f2, std::tuple{0, 1, 0}, f2, std::tuple{0, 0, 1});

        auto [x2, y2, z2] = types::process_func_arg_pairs(f3, std::tuple{1, 0, 0}, f3, std::tuple{0, 1, 0}, f3, std::tuple{0, 0, 1});
        
        // making a10 == 0
        auto a10_a00 = x1 / x0;
        // x1 = x1 - a10_a00 * x0; // x1 == 0
        y1 = y1 - a10_a00 * y0;
        z1 = z1 - a10_a00 * z0;
        std::get<1>(point) = std::get<1>(point) - a10_a00 * std::get<0>(point);
        // making a20 == 0
        auto a20_a00 = x2 / x0;
        //x2 = x2 - a20_a00 * x0; // x2 == 0
        y2 = y2 - a20_a00 * y0;
        z2 = z2 - a20_a00 * z0;
        std::get<2>(point) = std::get<2>(point) - a20_a00 * std::get<0>(point);

        // making a21 == 0
        auto a21_a11 = y2 / y1;
        // y2 = y2 - a21_a11 * y1; // y2 == 0
        z2 = z2 - a21_a11 * z1;
        std::get<2>(point) = std::get<2>(point) - a21_a11 * std::get<1>(point);

        float z = (float)std::get<2>(point) / z2;
        float y = (std::get<1>(point) - z1 * z) / y1;
        float x = (std::get<0>(point) - y0*y - z0*z) / x0;
        
        return std::tuple{x, y, z};
    }

    /// interactive method
    /// x0*x +y0*y + z0*z = p0
    /// x1*x +y1*y + z1*z = p1
    /// x2*x +y2*y + z2*z = p2
    /// becomes:
    /// f1 = 1/x0 * (p0 - y0*y -z0*z)
    /// f2 = 1/y1 * (p1 - x1*x - z1*z)
    /// f3 = 1/z2 *( p2 - x2*x - y2*y)
    /// with n interactions:
    /// x = f1(xn, yn, zn) // xn is not used so can be 0
    /// y = f2(xn, yn, zn) // yn is not used so can be 0
    /// z = f3(xn, yn, zn) // zn is not used so can be 0
    /// criteria of convergence:
    /// 1) strictly dominant diagonal:
    ///   | a00 a01 a02 |
    ///   | a10 a11 a12 |
    ///   | a20 a21 a22 |
    ///
    ///   |a00| > |a01| + |a02|
    ///   |a11| > |a10| + |a12|
    ///   |a22| > |a20| + |a21|
    ///   if true we have a solution. (the matrix converges)
    ///
    /// 2) dominant diagonal:
    ///   | a00 a01 a02 |
    ///   | a10 a11 a12 |
    ///   | a20 a21 a22 |
    ///
    ///   |a00| >= |a01| + |a02|
    ///   |a11| >= |a10| + |a12|
    ///   |a22| > |a20| + |a21|
    ///   if atleast 01 line is strictly dominant we have a solution. (the matrix converges)
    auto gauss_jacobi_method( auto&& f1, auto&& f2, auto&& f3, auto&& point, auto&& initialStimative, int interations = 6){
        // terms multipling x, y and z
        auto [x0, y0, z0] = process_func_arg_pairs(f1, std::tuple{1, 0, 0}, f1, std::tuple{0, 1, 0}, f1, std::tuple{0, 0, 1});

        auto [x1, y1, z1] = process_func_arg_pairs(f2, std::tuple{1, 0, 0}, f2, std::tuple{0, 1, 0}, f2, std::tuple{0, 0, 1});

        auto [x2, y2, z2] = process_func_arg_pairs(f3, std::tuple{1, 0, 0}, f3, std::tuple{0, 1, 0}, f3, std::tuple{0, 0, 1});
        
        auto [x, y, z] = jf::var::variables<3>();
            
        auto xx = (std::get<0>(point) -y0*y - z0*z) / x0;
        auto yy = (std::get<1>(point) -x1*x - z1*z) / y1;
        auto zz = (std::get<2>(point) -x2*x - y2*y) / z2;

        float inter1 = std::get<0>(initialStimative);
        float inter2 = std::get<1>(initialStimative);
        float inter3 = std::get<2>(initialStimative);
        float res1{};
        float res2{};
        float res3{};

        for (int i{}; i < interations; ++i) {
            res1 = xx(0, inter2, inter3);
            res2 = yy(inter1, 0, inter3);
            res3 = zz(inter1, inter2, 0);
            inter1 = res1;
            inter2 = res2;
            inter3 = res3;
        }
        
        // res1 = std::round(res1);
        // res2 = std::round(res2);
        // res3 = std::round(res3);
        return std::tuple{res1, res2, res3};
    }

    /// interactive method
    /// x0*x +y0*y + z0*z = p0
    /// x1*x +y1*y + z1*z = p1
    /// x2*x +y2*y + z2*z = p2
    /// becomes:
    /// f1 = 1/x0 * (p0 - y0*y -z0*z)
    /// f2 = 1/y1 * (p1 - x1*x - z1*z)
    /// f3 = 1/z2 *( p2 - x2*x - y2*y)
    /// with n interactions:
    /// xn+1 = f1(xn, yn, zn)     // xn is not used so can be 0
    /// yn+1 = f2(xn+1, yn, zn)   // yn is not used so can be 0
    /// zn+1 = f3(xn+1, yn+1, zn) // zn is not used so can be 0
    /// criteria of coonvergence:
    /// x0*x +y0*y + z0*z = p0
    /// x1*x +y1*y + z1*z = p1
    /// x2*x +y2*y + z2*z = p2
    /// B1 = (y0 + z0) / x0
    /// B2 = (x1*B1 + z1) / y1
    /// B3 = (x2*B1 + y2*B2) / z2
    /// if Bmax < 1, so the method generates a convergent sequence
    /// the smaller B, faster the convergence
    auto gauss_seidel_method( auto&& f1, auto&& f2, auto&& f3, auto&& point, auto&& initialStimative, int interations = 6){
        // terms multipling x, y and z
        auto [x0, y0, z0] = types::process_func_arg_pairs(f1, std::tuple{1, 0, 0}, f1, std::tuple{0, 1, 0}, f1, std::tuple{0, 0, 1});

        auto [x1, y1, z1] = types::process_func_arg_pairs(f2, std::tuple{1, 0, 0}, f2, std::tuple{0, 1, 0}, f2, std::tuple{0, 0, 1});

        auto [x2, y2, z2] = types::process_func_arg_pairs(f3, std::tuple{1, 0, 0}, f3, std::tuple{0, 1, 0}, f3, std::tuple{0, 0, 1});
        
        auto [x, y, z] = jf::var::variables<3>();
            
        auto xx = (std::get<0>(point) -y0*y - z0*z) / x0;
        auto yy = (std::get<1>(point) -x1*x - z1*z) / y1;
        auto zz = (std::get<2>(point) -x2*x - y2*y) / z2;

        float inter1 = std::get<0>(initialStimative);
        float inter2 = std::get<1>(initialStimative);
        float inter3 = std::get<2>(initialStimative);

        for (int i{}; i < interations; ++i) {
            inter1 = xx(0, inter2, inter3);
            inter2 = yy(inter1, 0, inter3);
            inter3 = zz(inter1, inter2, 0);
        }
        
        // inter1 = std::round(inter1);
        // inter2 = std::round(inter2);
        // inter3 = std::round(inter3);        
        return std::tuple{inter1, inter2, inter3};
    }


} // namespace jf::math


