module;

#include "std_input.hpp"

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
    constexpr auto cart_cyl(types::number_c auto&& Ax, types::number_c auto&& Ay, types::number_c auto&& Az,
            /* types::number_c auto rho, */ types::number_c auto phi /* types::number_c auto z */ )
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

    template<class Vec = std::tuple<
                std::function<double(double, double, double)>,
                std::function<double(double, double, double)>,
                std::function<double(double, double, double)>>,
            class Point = std::tuple<double, double, double>>
    constexpr auto cart_cyl(Vec&& vec, Point&& point)
    {
        double x = std::get<0>(point);
        double y = std::get<1>(point);
        double z = std::get<2>(point);
        auto [ rho, phi, zz ] = point_cart_cyl(x, y, z);

        jf::matrix::dmat var_cyl(3, 3);
        var_cyl.set_value({std::cos(phi), std::sin(phi),   0,
                        -std::sin(phi), std::cos(phi), 0,
                                    0,              0, 1 });

        jf::matrix::dmat coord_cart(3, 1);
        coord_cart.set_value({std::get<0>(vec)(x, y, z), std::get<1>(vec)(x, y, z), std::get<2>(vec)(x, y, z)});

        jf::matrix::dmat coord_cyl = var_cyl * coord_cart;

        return std::tuple{coord_cyl(0, 0), coord_cyl(1, 0), coord_cyl(2, 0)};

    }

    // ------- vector in cylindrical to a vector in cartesian  -----//
    constexpr auto cyl_cart(double rho, double phi, double z )
    {
        jf::matrix::dmat var_cyl(3, 3); // inverse
        var_cyl.set_value({std::cos(phi), std::sin(phi),   0,
                        -std::sin(phi), std::cos(phi), 0,
                                    0,              0, 1 });

        jf::matrix::dmat coord_cyl(3, 1);
        coord_cyl.set_value({rho, phi, z});

        jf::matrix::dmat coord_cart = var_cyl * coord_cyl;

        return std::tuple{coord_cart(0, 0), coord_cart(1, 0), coord_cart(2, 0)};

    }

    // ------- vector in cartesian to a vector in spherical  -----//
    constexpr auto cart_spher(double x, double y, double z,
           /* types::number_c auto r */ double theta, double phi  )
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

    // ------- vector in spherical to a vector in cartesian  -----//
    constexpr auto spher_cart(double Ar, double Atheta, double Aphi,
                           /* types::number_c auto r, */ double theta, double phi )
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
    
    // ------- vector in spherical to a vector in cylindrical  -----//
    constexpr auto spher_cyl(double Ar, double Atheta, double Aphi, double theta )
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

    // ------- vector in cylindrical to a vector in cartesian  -----//
    // constexpr auto cyl_cart(types::number_c auto rho, types::number_c auto phi, types::number_c auto z )
    // {
    //     jf::matrix::fmat var_cyl(3, 3); // inverse
    //     arg.set_value({std::cos(phi), std::sin(phi),   0,
    //                     -std::sin(phi), std::cos(phi), 0,
    //                                 0,              0, 1 });
    //
    //     jf::matrix::fmat coord_cyl(3, 1);
    //     coord_cart.set_value({rho, phi, z});
    //
    //     jf::matrix::fmat coord_cart = var_cyl * coord_cyl;
    //
    //     return std::tuple{coord_cart(0, 0), coord_cart(1, 0), coord_cart(2, 0)};
    //
    // }

}



