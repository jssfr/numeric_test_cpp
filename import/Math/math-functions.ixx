module;

#include<cmath>

export module math:Functions;

import Config;
import :Variable;

export namespace jf::math
{
    namespace functions{
        using namespace jf::var;
        template<class FuncType>
        inline constexpr auto ln(FuncType&& func){
            return [func](auto arg){ return std::log(func(arg)); };
        }
        template<class FuncType>
        inline constexpr auto cbrt(FuncType&& func)
        {
           return [func](auto arg){ return std::cbrt(func(arg)); }; 
        }


        template<class FuncType>
        inline constexpr auto sin(FuncType&& func){
                return [func](auto arg){ return std::sin(func(arg)); };
        }

        template<class FuncType>
        inline constexpr auto cos(FuncType&& func){
            return [func](auto arg){ return std::cos(func(arg)); };
        }
        
        // ---------------  operator + ---------------------

        

        template <class F, jf::types::integral_or_floating_point_c Gtype>
        auto operator+(F&& f, Gtype&& g) noexcept {
            return [=]<typename T, typename... ts>(T x0, ts... args) {
                return (validate(f, x0, args...) + g);
            };
        }

        template <jf::types::integral_or_floating_point_c Ftype, class G>
        auto operator+(Ftype&& f, G&& g) noexcept {
            return [=]<typename T, typename... ts>(T x0, ts... args) {
                return (f + validate(g, x0, args...));
            };
        }

        // --------------- operator -  ----------------------

        template <class F, jf::types::integral_or_floating_point_c Gtype>
        auto operator-(F&& f, Gtype&& g) noexcept {
            return [=]<typename T, typename... ts>(T x0, ts... args) {
                return (validate(f, x0, args...) - g);
            };
        }

        template <jf::types::integral_or_floating_point_c Ftype, class G>
        auto operator-(Ftype&& f, G&& g) noexcept {
            return [=]<typename T, typename... ts>(T x0, ts... args) {
                return (f - validate(g, x0, args...));
            };
        }
        // ---------------  operator * ---------------------


        template <class F, jf::types::integral_or_floating_point_c Gtype>
        auto operator*(F&& f, Gtype&& g) noexcept {
            return [=]<typename T, typename... ts>(T x0, ts... args) {
                return (validate(f, x0, args...) * g);
            };
        }

        template <jf::types::integral_or_floating_point_c Ftype, class G>
        auto operator*(Ftype&& f, G&& g) noexcept {
            return [=]<typename T, typename... ts>(T x0, ts... args) {
                return (f * validate(g, x0, args...));
            };
        }

    }


    template<class FuncType>
    inline constexpr auto ln(FuncType&& func){
        return functions::ln(func);
    }

    template<class FuncType>
    inline constexpr auto cbrt(FuncType&& func){
            return functions::cbrt(func);
    }


    template<class FuncType>
    inline constexpr auto sin(FuncType&& func){
            return functions::sin(func);
    }

    template<class FuncType>
    inline constexpr auto cos(FuncType&& func){
        return functions::cos(func);
    }

} // namespace jf::math
