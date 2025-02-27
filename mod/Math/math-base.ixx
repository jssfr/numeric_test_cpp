module;

#include "std_input.hpp"
#include<tuple>

export module math:Base;
import Config;
import :matrix2d;


 export namespace jf::math{
    namespace types = jf::types;

   constexpr auto point_cart_cyl(types::number_c auto x, types::number_c auto y, types::number_c auto z)
    {
        auto rho = std::hypot(x, y);//std::sqrt(x*x + y*y);
        auto phi = std::atan2(y, x);

        return std::tuple{rho, phi, z};
    }
    constexpr auto point_cyl_cart(types::number_c auto rho, types::number_c auto phi, types::number_c auto z)
    {
        auto x = rho*std::cos(phi);
        auto y = rho*std::sin(phi);

        return std::tuple{x, y, z};
    }

    constexpr auto point_cart_spher(types::number_c auto x, types::number_c auto y, types::number_c auto z)
    {
        auto r     = std::hypot(x, y, z);//std::sqrt(x*x + y*y + z*z);
        auto theta = std::atan2(std::hypot(x, y), z);
        auto phi   = std::atan2(y, x);

        return std::tuple{r, theta, phi};
    }
    constexpr auto point_spher_cart(types::number_c auto r, types::number_c auto theta, types::number_c auto phi)
    {
        auto x = r*std::sin(theta)*std::cos(phi);
        auto y = r*std::sin(theta)*std::sin(phi);
        auto z = r*std::cos(theta);

        return std::tuple{x, y, z};
    }

    constexpr auto point_cyl_spher(types::number_c auto rho, types::number_c auto phi, types::number_c auto z)
    {
        auto r = std::hypot(rho, z);
        auto theta = std::atan2(rho, z);

        return std::tuple{r, theta, phi};
    }

    constexpr auto point_spher_cyl(types::number_c auto r, types::number_c auto theta, types::number_c auto phi)
    {
        auto rho = r*std::sin(theta);
        auto z   = r*std::cos(theta);

        return std::tuple{rho, phi, z};
    }

    // -----   vector in cartesian to a vector in cylindrical ----- //
    auto cart_cyl(types::number_c auto&& Ax, types::number_c auto&& Ay, types::number_c auto&& Az,
            types::number_c auto phi  )
    {
        jf::matrix::dmat var_cyl(3, 3);
        var_cyl.set_value({std::cos(phi), std::sin(phi),   0,
                        -std::sin(phi), std::cos(phi), 0,
                                    0,              0, 1 });

        jf::matrix::dmat coord_cart(3, 1);
        coord_cart.set_value({static_cast<double>(Ax), static_cast<double>(Ay), static_cast<double>(Az)});

        jf::matrix::dmat coord_cyl = var_cyl * coord_cart;

        return std::tuple{coord_cyl(0, 0), coord_cyl(1, 0), coord_cyl(2, 0)};

    }

    template<class Point = std::tuple<double, double, double>>
    auto cart_cyl(Point&& point, auto&&... funcs)
    {
        auto [ rho, phi, zz ] = point_cart_cyl(std::get<0>(point), std::get<1>(point), std::get<2>(point));

        jf::matrix::dmat var_cyl(3, 3);
        var_cyl.set_value({std::cos(phi), std::sin(phi),   0,
                        -std::sin(phi), std::cos(phi), 0,
                                    0,              0, 1 });

        auto [ Ax, Ay, Az ] = types::process_func_arg( point, std::forward_as_tuple(funcs...));
        jf::matrix::dmat coord_cart(3, 1);
        
        coord_cart.set_value({static_cast<double>(Ax), static_cast<double>(Ay), static_cast<double>(Az)});
        jf::matrix::dmat coord_cyl = var_cyl * coord_cart;

        return std::tuple{coord_cyl(0, 0), coord_cyl(1, 0), coord_cyl(2, 0)};

    }

    // ------- vector in cylindrical to a vector in cartesian  -----//
    auto cyl_cart(double rho, double phi, double z )
    {
        jf::matrix::dmat var_cyl(3, 3); // inverse
        var_cyl.set_value({std::cos(phi), -std::sin(phi),   0,
                         std::sin(phi), std::cos(phi), 0,
                                    0,              0, 1 });

        jf::matrix::dmat coord_cyl(3, 1);
        coord_cyl.set_value({rho, phi, z});

        jf::matrix::dmat coord_cart = var_cyl * coord_cyl;

        return std::tuple{coord_cart(0, 0), coord_cart(1, 0), coord_cart(2, 0)};

    }

    template<class Point = std::tuple<double, double, double>>
    auto cyl_cart(Point&& point, auto&&... funcs)
    {
        double phi = std::get<1>(point);

        jf::matrix::dmat var_cyl(3, 3);
        var_cyl.set_value({std::cos(phi), -std::sin(phi),   0,
                         std::sin(phi), std::cos(phi), 0,
                                    0,              0, 1 });

        auto [ Arho, Aphi, Az ] = types::process_func_arg( point, std::forward_as_tuple(funcs...));
        jf::matrix::dmat coord_cyl(3, 1);
        
        coord_cyl.set_value({static_cast<double>(Arho), static_cast<double>(Aphi), static_cast<double>(Az)});
        jf::matrix::dmat coord_cart = var_cyl * coord_cyl;

        return std::tuple{coord_cart(0, 0), coord_cart(1, 0), coord_cart(2, 0)};

    }

    // ------- vector in cartesian to a vector in spherical  -----//
    auto cart_spher(double x, double y, double z,
           double theta, double phi  )
    {
        jf::matrix::dmat var_spher(3, 3);
        var_spher.set_value({std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta),
                       std::cos(theta)*std::cos(phi), std::cos(theta)*std::sin(phi), -std::sin(theta),
                                      -std::sin(phi),                std::cos(phi),  0 });

        jf::matrix::dmat coord_cart(3, 1);
        coord_cart.set_value({x, y, z});

        jf::matrix::dmat coord_spher = var_spher * coord_cart;

        return std::tuple{coord_spher(0, 0), coord_spher(1, 0), coord_spher(2, 0)};

    }

   
    auto cart_spher(std::tuple<double, double, double>&& point, auto&&... funcs)
    {
        auto [ r, theta, phi ] = point_cart_cyl(std::get<0>(point), std::get<1>(point), std::get<2>(point));

        jf::matrix::dmat var_spher(3, 3);
        var_spher.set_value({std::sin(theta)*std::cos(phi), std::sin(theta)*std::sin(phi), std::cos(theta),
                       std::cos(theta)*std::cos(phi), std::cos(theta)*std::sin(phi), -std::sin(theta),
                                      -std::sin(phi),                std::cos(phi),  0 });

        auto [ Ax, Ay, Az ] = types::process_func_arg( point, std::forward_as_tuple(funcs...));
        jf::matrix::dmat coord_cart(3, 1);
                                      
        coord_cart.set_value({static_cast<double>(Ax), static_cast<double>(Ay), static_cast<double>(Az)});
        jf::matrix::dmat coord_spher = var_spher * coord_cart;

        return std::tuple{coord_spher(0, 0), coord_spher(1, 0), coord_spher(2, 0)};

    }

    // ------- vector in spherical to a vector in cartesian  -----//
    auto spher_cart(double Ar, double Atheta, double Aphi,
                            double theta, double phi )
    {
        jf::matrix::dmat var_spher(3, 3); // inverse
        var_spher.set_value({std::sin(theta)*std::cos(phi), std::cos(theta)*std::cos(phi), -std::sin(phi),
                             std::sin(theta)*std::sin(phi), std::cos(theta)*std::sin(phi), std::cos(phi),
                                           std::cos(theta),              -std::sin(theta),  0 });

        jf::matrix::dmat coord_spher(3, 1);
        coord_spher.set_value({Ar, Atheta, Aphi});

        jf::matrix::dmat coord_cart = var_spher * coord_spher;

        return std::tuple{coord_cart(0, 0), coord_cart(1, 0), coord_cart(2, 0)};

    }

   
    auto spher_cart(std::tuple<double, double, double>&& point, auto&&... funcs)
    {
        double theta = std::get<1>(point);
        double phi   = std::get<2>(point);

        jf::matrix::dmat var_spher(3, 3); // inverse
        var_spher.set_value({std::sin(theta)*std::cos(phi), std::cos(theta)*std::cos(phi), -std::sin(phi),
                             std::sin(theta)*std::sin(phi), std::cos(theta)*std::sin(phi), std::cos(phi),
                                           std::cos(theta),              -std::sin(theta),  0 });

        auto [ Ar, Atheta, Aphi ] = types::process_func_arg( point, std::forward_as_tuple(funcs...));
        jf::matrix::dmat coord_spher(3, 1);
        
        coord_spher.set_value({static_cast<double>(Ar), static_cast<double>(Atheta), static_cast<double>(Aphi)});
        jf::matrix::dmat coord_cart = var_spher * coord_spher;

        return std::tuple{coord_cart(0, 0), coord_cart(1, 0), coord_cart(2, 0)};
    }

    
    // ------- vector in spherical to a vector in cylindrical  -----//
    auto spher_cyl(double Ar, double Atheta, double Aphi, double theta )
    {
        jf::matrix::dmat var_spher(3, 3); // inverse
        var_spher.set_value({std::sin(theta), std::cos(theta), 0,
                                            0,              0, 1,
                            std::cos(theta), -std::sin(theta), 0 });

        jf::matrix::dmat coord_spher(3, 1);
        coord_spher.set_value({Ar, Atheta, Aphi});

        jf::matrix::dmat coord_cyl = var_spher * coord_spher;

        return std::tuple{coord_cyl(0, 0), coord_cyl(1, 0), coord_cyl(2, 0)};

    }

    auto spher_cyl(std::tuple<double, double, double>&& point, auto&&... funcs)
    {
        double theta = std::get<1>(point);
      
        jf::matrix::dmat var_spher(3, 3); // inverse
        var_spher.set_value({std::sin(theta), std::cos(theta), 0,
                                            0,              0, 1,
                            std::cos(theta), -std::sin(theta), 0 });

        auto [ Ar, Atheta, Aphi ] = types::process_func_arg( point, std::forward_as_tuple(funcs...));
        jf::matrix::dmat coord_spher(3, 1);
        
        coord_spher.set_value({static_cast<double>(Ar), static_cast<double>(Atheta), static_cast<double>(Aphi)});
        jf::matrix::dmat coord_cyl = var_spher * coord_spher;

        return std::tuple{coord_cyl(0, 0), coord_cyl(1, 0), coord_cyl(2, 0)};
    }

     // ------- vector in cylindrical to a vector in spherical  -----//
     auto cyl_spher(double Arho, double Aphi, double z, double theta )
     {
         jf::matrix::dmat var_cyl(3, 3); // inverse
         var_cyl.set_value({std::sin(theta),  0, std::cos(theta),
                            std::cos(theta),    0, -std::sin(theta),
                                        0,      1,   0 });
 
         jf::matrix::dmat coord_cyl(3, 1);
         coord_cyl.set_value({Arho, Aphi, z});
 
         jf::matrix::dmat coord_spher = var_cyl * coord_cyl;
 
         return std::tuple{coord_spher(0, 0), coord_spher(1, 0), coord_spher(2, 0)};
 
     }
 
     auto cyl_spher(std::tuple<double, double, double>&& point, auto&&... funcs)
     {
        auto [ r, theta, phi ] = point_cyl_spher(std::get<0>(point), std::get<1>(point), std::get<2>(point));
       
        jf::matrix::dmat var_cyl(3, 3); // inverse
        var_cyl.set_value({std::sin(theta),  0, std::cos(theta),
                           std::cos(theta),    0, -std::sin(theta),
                                       0,      1,   0 });

 
         auto [ Arho, Aphi, Az ] = types::process_func_arg( point, std::forward_as_tuple(funcs...));
         jf::matrix::dmat coord_cyl(3, 1);
         
         coord_cyl.set_value({static_cast<double>(Arho), static_cast<double>(Aphi), static_cast<double>(Az)});
         jf::matrix::dmat coord_spher = var_cyl * coord_cyl;
 
         return std::tuple{coord_spher(0, 0), coord_spher(1, 0), coord_spher(2, 0)};
     }

}



