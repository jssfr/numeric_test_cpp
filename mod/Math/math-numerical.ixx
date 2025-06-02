module;
#ifndef USING_IMPORT_STD_MOD
  #include"std_input.hpp"
#endif

export module math:Numerical;
#ifdef USING_IMPORT_STD_MOD
  import std;
#endif
import Config;

import :Variable;
import :Calculus;


namespace numConst{

    constexpr int stop = 15;
    constexpr auto signal = []<jf::types::number_c Argtype>( const Argtype& arg)
                { 
                    return (arg<0? -1 : 1);
                };


} //  namespace numConst


export namespace jf::math{

    namespace types = jf::types;

    template<types::number_c Number, class FuncType>
    constexpr auto bissec(Number interval1, Number interval2, FuncType&& func, double precision = 0.000005)
    {
        
        if(func(interval1) * func(interval2) > 0){ throw("can´t calculate because: func(interval1) * func(interval2) > 0. "); }
        
        double result = 0.0;
        for(int count{}; count < numConst::stop; ++count)
        {
            result = (interval1 + interval2) / 2;
            if(std::abs(func(result)) < precision){ return result; }

            if(numConst::signal(func(result)) == numConst::signal(func(interval1))){ interval1 = result; }
            else{interval2 = result; } 
        }

        return result;

    }// bissec
     
    template<types::number_c Number, class FuncType>
    constexpr auto false_position(Number interval1, Number interval2, FuncType&& func, double precision = 0.000005)
    {
       auto fun1 = func(interval1);
       auto fun2 = func(interval2);
        if(fun1 * fun2 > 0){ throw( "can´t calculate because: func(interval1) * func(interval2) > 0. "); }
        
        
        double result = 0.0;
        for(int count{}; count < numConst::stop; ++count)
        {
            result = (interval1*fun2 - interval2*fun1) / (fun2 - fun1);
            if(std::abs(func(result)) < precision){ return result; }

            if(numConst::signal(func(result)) == numConst::signal(func(interval1))){ interval1 = result; }
            else{interval2 = result; } 

            fun1 = func(interval1);
            fun2 = func(interval2);       
        }

        return result;

    }// false_position
    
    /*
     * take a function of type f(x) = 0 and  a x0 the first extimative 
     * @return the first repeate x decimal precision*/
    template<types::number_c Number, class FuncType>
    constexpr auto fixed_point(FuncType&& func, Number x0, int order = 4)
    {
     
        double result = 0.0;
        double x1 = func(x0);

        double trunc1 = std::trunc(x1 * std::pow(10, order)) / std::pow(10, order);
        double trunc2 = 0.0;

        for(int count{}; count < numConst::stop; ++count)
        {
            result = func(x1);
            trunc2 = std::trunc(result * std::pow(10, order)) / std::pow(10, order);
            if(trunc1 == trunc2){ return result; }

            x1 = result;
            trunc1 = trunc2;

        }

        return result;

    }// fixed_point

    template<types::number_c Number, class FuncType>
    constexpr auto newton_method(FuncType&& func, Number x0, double funcDelimiter = 0.000001, int order = 4)
    {

        double result = 0.0;
        double x1 = x0 - func(x0) / jf::math::five_point_stencil<1>(func, x0);

        double trunc1 = std::trunc(x1 * std::pow(10, order)) / std::pow(10, order);
        double trunc2 = 0.0;

        for(int count{}; count < numConst::stop; ++count)
        {
            result = x1 - func(x1) / jf::math::five_point_stencil<1>(func, x1);
            trunc2 = std::trunc(result * std::pow(10, order)) / std::pow(10, order);
            if(trunc1 == trunc2 || std::abs(func(result)) < funcDelimiter){ return result; }

            x1 = result;
            trunc1 = trunc2;

        }

        return result;
    } // newton_method

    template<types::number_c Number, class FuncType>
    constexpr auto secant_method(FuncType&& func, Number x0, Number x1, double funcDelimiter = 0.000001)
    {

        auto fun0 = func(x0);
        auto fun1 = func(x1);
        double result = 0.0;


        for(int count{}; count < numConst::stop; ++count)
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
    auto gauss_method( auto&& f1, auto&& f2, auto&& f3, types::tuple_or_array_c auto&& const_terms){
        // terms multipling x, y and z
        auto [x0, y0, z0] = types::process_func_arg_pairs(f1, std::tuple{1.f, 0.f, 0.f}, f1, std::tuple{0.f, 1.f, 0.f}, f1, std::tuple{0.f, 0.f, 1.f});

        auto [x1, y1, z1] = types::process_func_arg_pairs(f2, std::tuple{1.f, 0.f, 0.f}, f2, std::tuple{0.f, 1.f, 0.f}, f2, std::tuple{0.f, 0.f, 1.f});

        auto [x2, y2, z2] = types::process_func_arg_pairs(f3, std::tuple{1.f, 0.f, 0.f}, f3, std::tuple{0.f, 1.f, 0.f}, f3, std::tuple{0.f, 0.f, 1.f});
        
        // making a10 == 0
        float init0 = std::get<0>(const_terms);
        float init1 = std::get<1>(const_terms);
        float init2 = std::get<2>(const_terms);
        // making a10 == 0
        float a10_a00 = x1 / x0;
        // x1 == 0
        y1 = y1 - a10_a00 * y0;
        z1 = z1 - a10_a00 * z0;
        init1 = init1 - a10_a00 * init0;
        // making a20 == 0
        float a20_a00 = x2 / x0;
        // x2 == 0
        y2 = y2 - a20_a00 * y0;
        z2 = z2 - a20_a00 * z0;
        init2 = init2 - a20_a00 * init0;

        // making a21 == 0
        float a21_a11 = y2 / y1;
        // y2 == 0
        z2 = z2 - a21_a11 * z1;
        init2 = init2 - a21_a11 * init1;

        float z = init2 / z2;
        float y = (init1 - z1 * z) / y1;
        float x = (init0 - y0*y - z0*z) / x0;
        
        return std::tuple{x, y, z};
    }

    auto LU_method( auto&& f1, auto&& f2, auto&& f3, types::tuple_or_array_c auto&& const_terms){
        // terms multipling x, y and z

        auto [x0, y0, z0] = types::process_func_arg_pairs(f1, std::tuple{1.f, 0.f, 0.f}, f1, std::tuple{0.f, 1.f, 0.f}, f1, std::tuple{0.f, 0.f, 1.f});

        auto [x1, y1, z1] = types::process_func_arg_pairs(f2, std::tuple{1.f, 0.f, 0.f}, f2, std::tuple{0.f, 1.f, 0.f}, f2, std::tuple{0.f, 0.f, 1.f});

        auto [x2, y2, z2] = types::process_func_arg_pairs(f3, std::tuple{1.f, 0.f, 0.f}, f3, std::tuple{0.f, 1.f, 0.f}, f3, std::tuple{0.f, 0.f, 1.f});
        
        float init0 = std::get<0>(const_terms);
        float init1 = std::get<1>(const_terms);
        float init2 = std::get<2>(const_terms);
        
        // making a10 == 0
        float a10_a00 = x1 / x0;
        y1 = y1 - a10_a00 * y0;
        z1 = z1 - a10_a00 * z0;
        init1 = init1 - a10_a00 * init0;
        // making a20 == 0
        float a20_a00 = x2 / x0;
        y2 = y2 - a20_a00 * y0;
        z2 = z2 - a20_a00 * z0;
        init2 = init2 - a20_a00 * init0;

        // making a21 == 0
        float a21_a11 = y2 / y1;
        z2 = z2 - a21_a11 * z1;
        init2 = init2 - a21_a11 * init1;

        // first calculating using lower triangular
        float yy1 = std::get<0>(const_terms);
        float yy2 = std::get<1>(const_terms) - (yy1 * a10_a00);
        float yy3 = std::get<2>(const_terms) - (a20_a00 * yy1) - (a21_a11 * yy2);

        // using the result yy- for find the final result
        float z = yy3 / z2;
        float y = (yy2 - z1 * z) / y1;
        float x = (yy1 - y0*y - z0*z) / x0;
        
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
    auto gauss_jacobi_method( auto&& f1, auto&& f2, auto&& f3, 
            types::tuple_or_array_c auto&& const_terms, types::tuple_or_array_c auto&& initialStimative, int interations = 6){
        // terms multipling x, y and z
        auto [x0, y0, z0] = types::process_func_arg_pairs(f1, std::tuple{1.f, 0.f, 0.f}, f1, std::tuple{0.f, 1.f, 0.f}, f1, std::tuple{0.f, 0.f, 1.f});

        auto [x1, y1, z1] = types::process_func_arg_pairs(f2, std::tuple{1.f, 0.f, 0.f}, f2, std::tuple{0.f, 1.f, 0.f}, f2, std::tuple{0.f, 0.f, 1.f});

        auto [x2, y2, z2] = types::process_func_arg_pairs(f3, std::tuple{1.f, 0.f, 0.f}, f3, std::tuple{0.f, 1.f, 0.f}, f3, std::tuple{0.f, 0.f, 1.f});
        
       auto abs_ge = [](auto val1, auto val2, auto val3){
            if(std::abs(val1) > (std::abs(val2) + std::abs(val3))){
                return 1;
            }else if(std::abs(val1) == (std::abs(val2) + std::abs(val3))){
                return 0;
            }else{
                return -2;
            }
        };

       int row1 = abs_ge(x0, y0, z0);
       int row2 = abs_ge(y1, x1, z1);
       int row3 = abs_ge(z2, x2, y2);

       // if gjm == 3: criteria 1), if gjm == 1 or 2: criteria 2), if gjm < 1 no convergence
       if(int gjm = row1 + row2 + row3; gjm < 1){
            return LU_method(f1, f2, f3, const_terms);
       }
        auto [x, y, z] = jf::var::variables<3>();
            
        auto xx = (std::get<0>(const_terms) -y0*y - z0*z) / x0;
        auto yy = (std::get<1>(const_terms) -x1*x - z1*z) / y1;
        auto zz = (std::get<2>(const_terms) -x2*x - y2*y) / z2;

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
    auto gauss_seidel_method( auto&& f1, auto&& f2, auto&& f3, 
            types::tuple_or_array_c auto&& const_terms, types::tuple_or_array_c auto&& initialStimative, int interations = 6)
    {
        // terms multipling x, y and z
        auto [x0, y0, z0] = types::process_func_arg_pairs(f1, std::tuple{1.f, 0.f, 0.f}, f1, std::tuple{0.f, 1.f, 0.f}, f1, std::tuple{0.f, 0.f, 1.f});

        auto [x1, y1, z1] = types::process_func_arg_pairs(f2, std::tuple{1.f, 0.f, 0.f}, f2, std::tuple{0.f, 1.f, 0.f}, f2, std::tuple{0.f, 0.f, 1.f});

        auto [x2, y2, z2] = types::process_func_arg_pairs(f3, std::tuple{1.f, 0.f, 0.f}, f3, std::tuple{0.f, 1.f, 0.f}, f3, std::tuple{0.f, 0.f, 1.f});
        
        auto B1 = (y0 + z0) / x0;
        auto B2 = (x1*B1 + z1) / y1;
        auto B3 = (x2*B1 + y2*B2) / z2;
        auto Bmax = (B1 > B2 || B1 > B3)? B1 : (B2 > B3)? B2 : B3;
        if(Bmax >= 1){
             return LU_method( f1, f2, f3, const_terms);
        }

        auto [x, y, z] = jf::var::variables<3>();
            
        auto xx = (std::get<0>(const_terms) -y0*y - z0*z) / x0;
        auto yy = (std::get<1>(const_terms) -x1*x - z1*z) / y1;
        auto zz = (std::get<2>(const_terms) -x2*x - y2*y) / z2;

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

    class system3_solver
    {
        private:
            std::tuple<
                std::function<float(float, float, float)>,
                std::function<float(float, float, float)>,
                std::function<float(float, float, float)>
            > systems;
            std::tuple<float, float, float> rhds;
        public:

            system3_solver(auto&& x0, auto&& y0, auto&& z0, types::number_c auto rhs1,
                           auto&& x1, auto&& y1, auto&& z1, types::number_c auto rhs2,
                           auto&& x2, auto&& y2, auto&& z2, types::number_c auto rhs3
                          ) : systems{ x0 + y0 + z0, x1 + y1 + z1, x2 + y2 + z2 }, 
                rhds{ static_cast<float>(rhs1), static_cast<float>(rhs2), static_cast<float>(rhs3) }
            {  }

            auto solve() -> decltype(auto)
            {
                return LU_method(std::get<0>(systems), std::get<1>(systems), std::get<2>(systems), rhds);
            }
            auto solve_gauss() -> decltype(auto)
            {
                return gauss_method(std::get<0>(systems), std::get<1>(systems), std::get<2>(systems), rhds);
            }
            auto solve_gauss_seidel(types::tuple_or_array_c auto&& initialStimative, int interations = 6) -> decltype(auto)
            {
                return gauss_seidel_method(std::get<0>(systems), std::get<1>(systems),
                                        std::get<2>(systems), rhds, initialStimative, interations);
            }
            auto solve_gauss_jacobi(types::tuple_or_array_c auto&& initialStimative, int interations = 6) -> decltype(auto)
            {
                return gauss_jacobi_method(std::get<0>(systems), std::get<1>(systems),
                                        std::get<2>(systems), rhds, initialStimative, interations);
            }
    };


} // namespace jf::math


