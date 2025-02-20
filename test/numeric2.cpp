
#ifdef _MSC_VER
    #pragma warning(disable: 5050)
#endif


#include "std_input.hpp"
#include<print>

import Config;
import math;

/// Gauss method, gauss-jacob method, gauss-seidel method
/// solution method for solving linear systems equations


auto process_func_arg_pairs(auto ...args) -> decltype(auto)
    requires (sizeof...(args) % 2 == 0)
{
    constexpr auto N = sizeof...(args);

    auto process_func = [arguments = std::make_tuple(args...)]<auto ...i>(std::index_sequence<i...>){
        return std::tuple{ std::apply(std::get<i>(arguments), std::get<i+1>(arguments)) ... };
    };

    return jf::types::lambda_seq<0, N, 2>(process_func);
}

/// equations to make zero the terms under pricipal diagonal
/// under a00:
/// L1 = L1 - a10/a00 * L0 ==> L0 = line0 
/// L2 = L2 - a20/a00 * L0
/// under a11:
/// L2 = L2 - a21/a11 * L1
auto gauss_method( auto&& f1, auto&& f2, auto&& f3, auto&& point){
    // terms multipling x, y and z

    auto [x0, y0, z0] = process_func_arg_pairs(f1, std::tuple{1.f, 0, 0}, f1, std::tuple{0, 1.f, 0}, f1, std::tuple{0, 0, 1.f});

    auto [x1, y1, z1] = process_func_arg_pairs(f2, std::tuple{1.f, 0, 0}, f2, std::tuple{0, 1.f, 0}, f2, std::tuple{0, 0, 1.f});

    auto [x2, y2, z2] = process_func_arg_pairs(f3, std::tuple{1.f, 0, 0}, f3, std::tuple{0, 1.f, 0}, f3, std::tuple{0, 0, 1.f});
    
    float init0 = std::get<0>(point);
    float init1 = std::get<1>(point);
    float init2 = std::get<2>(point);
    // making a10 == 0
    float a10_a00 = x1 / x0;
    // x1 = x1 - a10_a00 * x0; // x1 == 0
    y1 = y1 - a10_a00 * y0;
    z1 = z1 - a10_a00 * z0;
    init1 = init1 - a10_a00 * init0;
    // making a20 == 0
    float a20_a00 = x2 / x0;
    //x2 = x2 - a20_a00 * x0; // x2 == 0
    y2 = y2 - a20_a00 * y0;
    z2 = z2 - a20_a00 * z0;
    init2 = init2 - a20_a00 * init0;

    // making a21 == 0
    float a21_a11 = y2 / y1;
    // y2 = y2 - a21_a11 * y1; // y2 == 0
    z2 = z2 - a21_a11 * z1;
    init2 = init2 - a21_a11 * init1;

    float z = init2 / z2;
    float y = (init1 - z1 * z) / y1;
    float x = (init0 - y0*y - z0*z) / x0;
    
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
///   with true we have a solution. (the matrix converges)
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

auto test_gauss_method() -> void
{
    auto [x, y, z] = jf::var::variables<3>();

    auto f1 = 2*x + 3*y -z;   //  = 5
    auto f2 = 4*x +4*y -3*z; //   = 3
    auto f3 = 2*x -3*y + z; //  = -1
    auto [xx, yy, zz] = gauss_method(f1, f2, f3, std::tuple{5.f, 3.f, -1.f});
    // auto f1 = 5*x - 2*y +1*z;   //  = 5
    // auto f2 = 1*x +4*y +2*z; //   = 3
    // auto f3 = -2*x +1*y + 6*z; //  = -1
    // 
    // auto [xx, yy, zz] = gauss_method(f1, f2, f3, std::tuple{12.f, 2.f, 9.f});

    std::print("x = {0}, y = {1}, z = {2}\n\n", xx, yy, zz);
}

auto test_gauss_jacobi_method() -> void{
    auto [x, y, z] = jf::var::variables<3>();

    auto f1 = 10*x +2*y + z;
    auto f2 = x + 5*y + z;
    auto f3 = 2*x +3*y + 10*z;

    auto [xx, yy, zz] = gauss_jacobi_method(f1, f2, f3, std::tuple{7, -8, 6}, std::tuple{0.5, -1, 0.5});

    // auto f1 = 5*x - 2*y +z;   //  = 5
    // auto f2 = 1*x +4*y +2*z; //   = 3
    // auto f3 = -2*x +1*y + 6*z; //  = -1

    // auto [xx, yy, zz] = gauss_jacobi_method(f1, f2, f3, std::tuple{12, 2, 9}, std::tuple{0.5, -1, 0.5}, 2);
    std::print("x = {}, y = {}, z = {}\n\n", xx, yy, zz);
}

auto test_gauss_seidel_method() -> void{
    auto [x, y, z] = jf::var::variables<3>();

    auto f1 = 10*x +2*y + z;
    auto f2 = x + 5*y + z;
    auto f3 = 2*x +3*y + 10*z;
    auto [xx, yy, zz] = gauss_seidel_method(f1, f2, f3, std::tuple{7, -8, 6}, std::tuple{0.5, -1, 0.5});
    // auto f1 = 5*x - 2*y +z;   //  = 5
    // auto f2 = 1*x +4*y +2*z; //   = 3
    // auto f3 = -2*x +1*y + 6*z; //  = -1
    //
    // auto [xx, yy, zz] = gauss_seidel_method(f1, f2, f3, std::tuple{12, 2, 9}, std::tuple{0.5, -1, 0.5}, 2);
    std::print("x = {}, y = {}, z = {}\n\n", xx, yy, zz);
}

int abs_ge(auto val1, auto val2, auto val3){
    if(std::abs(val1) > (std::abs(val2) + std::abs(val3))){
        return 1;
    }else if(std::abs(val1) == (std::abs(val2) + std::abs(val3))){
        return 0;
    }else{
        return -2;
    }
}


auto linear_solver(auto&& x0, auto&& y0, auto&& z0, jf::types::number_c auto&& rhs1,
                   auto&& x1, auto&& y1, auto&& z1, jf::types::number_c auto&& rhs2,
                   auto&& x2, auto&& y2, auto&& z2, jf::types::number_c auto&& rhs3
                    ) -> decltype(auto)
{
   return gauss_method(x0 + y0 + z0, x1 + y1 + z1, x2 + y2 + z2, std::tuple{rhs1, rhs2, rhs3} );
}

template<int GjorGs = 1>
auto linear_solver(auto&& x0, auto&& y0, auto&& z0, jf::types::number_c auto&& rhs1,
                   auto&& x1, auto&& y1, auto&& z1, jf::types::number_c auto&& rhs2,
                   auto&& x2, auto&& y2, auto&& z2, jf::types::number_c auto&& rhs3,
                   jf::types::tuple_or_array_c auto&& initialStimative, int inter = 6) -> decltype(auto)
{

   auto [xx0, yy0, zz0] = process_func_arg_pairs(x0, std::tuple{1, 0, 0}, y0, std::tuple{0, 1, 0}, z0, std::tuple{0, 0, 1});

   auto [xx1, yy1, zz1] = process_func_arg_pairs(x1, std::tuple{1, 0, 0}, y1, std::tuple{0, 1, 0}, z1, std::tuple{0, 0, 1});

   auto [xx2, yy2, zz2] = process_func_arg_pairs(x2, std::tuple{1, 0, 0}, y2, std::tuple{0, 1, 0}, z2, std::tuple{0, 0, 1});
    
   int row1 = abs_ge(xx0, yy0, zz0);
   int row2 = abs_ge(yy1, xx1, zz1);
   int row3 = abs_ge(zz2, xx2, yy2);

   if(int gjm = row1 + row2 + row3; gjm >= 1){
       std::println("gjm called.");
        return gauss_jacobi_method(x0 + y0 + z0, x1 + y1 + z1, x2 + y2 + z2, std::tuple{rhs1, rhs2, rhs3}, initialStimative, inter);
   }
   auto B1 = (yy0 + zz0) / xx0;
   auto B2 = (xx1*B1 + zz1) / yy1;
   auto B3 = (xx2*B1 + yy2*B2) / zz2;
   auto Bmax = (B1 > B2 || B1 > B3)? B1 : (B2 > B3)? B2 : B3;
   if(Bmax < 1){
       std::println("gsm called.");
        return gauss_seidel_method(x0 + y0 + z0, x1 + y1 + z1, x2 + y2 + z2, std::tuple{rhs1, rhs2, rhs3}, initialStimative, inter);
   }
    else{
        std::println("gm called.");
        return gauss_method( x0 + y0 + z0, x1 + y1 + z1, x2 + y2 + z2, std::tuple{rhs1, rhs2, rhs3});
   }
}
void test_solver(){

    auto [x, y, z] = jf::var::variables<3>();

    // calls gauss_method
    auto [xx, yy, zz] = linear_solver(2*x, 3*y, -1*z, 5,
                                      4*x, 4*y, -3*z, 3,
                                      2*x, -3*y, z, -1);
    std::print("xx = {}, yy = {}, zz = {}\n\n", xx, yy, zz);

    // calls gauss_jacobi_method
    auto [xj, yj, zj] = linear_solver(10*x, 2*y,   z, 7,
                                         x, 5*y,   z, -8,
                                       2*x, 3*y, 10*z, 6,
                                      std::tuple{0.5, -1, 0.5});
    std::print("xj = {}, yj = {}, zj = {}\n\n", xj, yj, zj);

    // 
    auto [xs, ys, zs] = linear_solver(21*x, -9*y, -12*z, -33,
                                      -3*x, 6*y, -2*z, 3,
                                      -8*x, -4*y, 22*z, 50);
    std::print("xs = {}, ys = {}, zs = {}\n\n", xs, ys, zs);

    auto [xr, yr, zr] = linear_solver(21*x, -9*y, -12*z, -33,
                                      -3*x, 6*y, -2*z, 3,
                                      -8*x, -4*y, 22*z, 50, 
                                      std::tuple{0.5, -1, 0.5}, 8);
    std::print("xr = {}, yr = {}, zr = {}\n\n", xr, yr, zr);

}


auto LU_method( auto&& f1, auto&& f2, auto&& f3, auto&& point){
    // terms multipling x, y and z

    auto [x0, y0, z0] = process_func_arg_pairs(f1, std::tuple{1.f, 0, 0}, f1, std::tuple{0, 1.f, 0}, f1, std::tuple{0, 0, 1.f});

    auto [x1, y1, z1] = process_func_arg_pairs(f2, std::tuple{1.f, 0, 0}, f2, std::tuple{0, 1.f, 0}, f2, std::tuple{0, 0, 1.f});

    auto [x2, y2, z2] = process_func_arg_pairs(f3, std::tuple{1.f, 0, 0}, f3, std::tuple{0, 1.f, 0}, f3, std::tuple{0, 0, 1.f});
    
    float init0 = std::get<0>(point);
    float init1 = std::get<1>(point);
    float init2 = std::get<2>(point);
    
    // upper triangular
    // | a00 a01 a02 |
    // |  0  a11 a12 |
    // |  0   0  a22 |
    
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
    
    // lower triangular
    // | 1   0    0 |
    // |a10  1    0 |
    // |a20  a21  1 |
    // a10 = a10_a00,
    // a20 = a20_a00,
    // a21 = a21_a11

    // first calculating using lower triangular
    float yy1 = std::get<0>(point);
    float yy2 = std::get<1>(point) - (yy1 * a10_a00);
    float yy3 = std::get<2>(point) - (a20_a00 * yy1) - (a21_a11 * yy2);

    // using the result yy- for find the final result
    float z = yy3 / z2;
    float y = (yy2 - z1 * z) / y1;
    float x = (yy1 - y0*y - z0*z) / x0;
    
    return std::tuple{x, y, z};
}

auto test_LU_method() -> void
{
    auto [x, y, z] = jf::var::variables<3>();

    auto f1 = 2*x + 3*y -z;   //  = 5
    auto f2 = 4*x +4*y -3*z; //   = 3
    auto f3 = 2*x -3*y + z; //  = -1
    auto [xx, yy, zz] = LU_method(f1, f2, f3, std::tuple{5.f, 3.f, -1.f});

    std::print("x = {0}, y = {1}, z = {2}\n\n", xx, yy, zz);

    // auto g1 = 10*x +2*y + z;
    // auto g2 = x + 5*y + z;
    // auto g3 = 2*x +3*y + 10*z;
    // auto [gx, gy, gz] = LU_method(g1, g2, g3, std::tuple{7.f, -8.f, 6.f});
    //
    // std::print("x = {0}, y = {1}, z = {2}\n\n", gx, gy, gz);
    auto h1 = 21*x -9*y -12*z;
    auto h2 = -3*x + 6*y -2*z;
    auto h3 = -8*x -4*y + 22*z;
    auto [hx, hy, hz] = LU_method(h1, h2, h3, std::tuple{-33.f, 3.f, 50.f});

    std::print("x = {0}, y = {1}, z = {2}\n\n", hx, hy, hz);
}

void test_system_solver(){

    auto [x, y, z] = jf::var::variables<3>();

    // calls gauss_method
    jf::math::system3_solver s1(2*x, 3*y, -1*z, 5,
                                      4*x, 4*y, -3*z, 3,
                                      2*x, -3*y, z, -1);
    auto [xx, yy, zz] = s1.solve();
    std::print("xx = {}, yy = {}, zz = {}\n\n", xx, yy, zz);

    // calls gauss_jacobi_method
    jf::math::system3_solver s2(10*x, 2*y,   z, 7,
                                         x, 5*y,   z, -8,
                                       2*x, 3*y, 10*z, 6);
    auto [xj, yj, zj] = s2.solve_gauss_jacobi(std::tuple{0.5, -1, 0.5});
    std::print("xj = {}, yj = {}, zj = {}\n\n", xj, yj, zj);

    // 
    jf::math::system3_solver s3(21*x, -9*y, -12*z, -33,
                                -3*x, 6*y, -2*z, 3,
                                -8*x, -4*y, 22*z, 50);
    auto [xs, ys, zs] = s3.solve_gauss();
    std::print("xs = {}, ys = {}, zs = {}\n\n", xs, ys, zs);

    jf::math::system3_solver s4(21*x, -9*y, -12*z, -33,
                                -3*x, 6*y, -2*z, 3,
                                -8*x, -4*y, 22*z, 50);
    auto [xr, yr, zr] = s4.solve_gauss_seidel(std::tuple{1.0, 2, 1.8}, 18);
    std::print("xr = {}, yr = {}, zr = {}\n\n", xr, yr, zr);

}
int main()
{
    // std::println("test gauss method");
    // test_gauss_method();
    // std::println("test gauss jacobi method");
    // test_gauss_jacobi_method();
    // std::println("test gauss seidel method");
    // test_gauss_seidel_method();
    // std::println("test solver");
    // test_solver();
    // std::println("test LU_method");
    // test_LU_method();
    test_system_solver();

    return 0;
}
