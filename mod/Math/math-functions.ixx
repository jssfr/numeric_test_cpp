module;

#ifndef USING_IMPORT_STD_MOD
  #include<cmath>
#endif

export module math:Functions;
#ifdef USING_IMPORT_STD_MOD
  import std;
#endif
import Config;
import :Variable;


namespace jf::var{
    // using namespace jf::var;
    template<class FuncType>
    inline constexpr auto ln(FuncType&& func){
        //return [&func](auto arg){ return std::log(func(arg)); };
        return [=]<typename T, typename... ts>(T x0, ts... args) {
            return std::log(func(x0, args...));
        };
    }
    template<class FuncType>
    inline constexpr auto cbrt(FuncType&& func)
    {
       // return [&func](auto arg){ return std::cbrt(func(arg)); }; 
        return [=]<typename T, typename... ts>(T x0, ts... args) {
            return std::cbrt(func(x0, args...));
        };
    }

    template<class FuncType>
    inline constexpr auto sqrt(FuncType&& func)
    {
       // return [&func](auto arg){ return std::cbrt(func(arg)); }; 
        return [=]<typename T, typename... ts>(T x0, ts... args) {
            return std::sqrt(func(x0, args...));
        };
    }

    template<class FuncType>
    inline constexpr auto sin(FuncType&& func){
         //   return [&func](auto arg){ return std::sin(func(arg)); };
        return [=]<typename T, typename... ts>(T x0, ts... args) {
            return std::sin(func(x0, args...));
        };
    }

    template<class FuncType>
    inline constexpr auto cos(FuncType&& func){
        //return [&func](auto arg){ return std::cos(func(arg)); };
        return [=]<typename T, typename... ts>(T x0, ts... args) {
            return std::cos(func(x0, args...));
        };
    }
    
    // // ---------------  operator + ---------------------
    //
    // 
    //
    // template <class F, jf::types::number_c Gtype>
    // auto operator+(F&& f, Gtype&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (validate(f, x0, args...) + g);
    //     };
    // }
    //
    // template <jf::types::number_c Ftype, class G>
    // auto operator+(Ftype&& f, G&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (f + validate(g, x0, args...));
    //     };
    // }
    //
    // // --------------- operator -  ----------------------
    //
    // template <class F, jf::types::number_c Gtype>
    // auto operator-(F&& f, Gtype&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (validate(f, x0, args...) - g);
    //     };
    // }
    //
    // template <jf::types::number_c Ftype, class G>
    // auto operator-(Ftype&& f, G&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (f - validate(g, x0, args...));
    //     };
    // }
    // // ---------------  operator * ---------------------
    //
    //
    // template <class F, jf::types::number_c Gtype>
    // auto operator*(F&& f, Gtype&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (validate(f, x0, args...) * g);
    //     };
    // }
    //
    // template <jf::types::number_c Ftype, class G>
    // auto operator*(Ftype&& f, G&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (f * validate(g, x0, args...));
    //     };
    // }
    // // ---------------  operator / ---------------------
    //
    //
    // template <class F, jf::types::number_c Gtype>
    // auto operator/(F&& f, Gtype&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (validate(f, x0, args...) / g);
    //     };
    // }
    //
    // template <jf::types::number_c Ftype, class G>
    // auto operator/(Ftype&& f, G&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (f / validate(g, x0, args...));
    //     };
    // }
    // template <class Ftype, class G>
    // auto operator/(Ftype&& f, G&& g) noexcept {
    //     return [=]<typename T, typename... ts>(T x0, ts... args) {
    //         return (validate(f, x0, args...) / validate(g, x0, args...));
    //     };
    // }

}




export namespace jf::math
{

    template<class FuncType>
    inline constexpr auto ln(FuncType&& func){
        return jf::var::ln(func);
    }

    template<class FuncType>
    inline constexpr auto cbrt(FuncType&& func){
            return jf::var::cbrt(func);
    }

    template<class FuncType>
    inline constexpr auto sqrt(FuncType&& func){
            return jf::var::sqrt(func);
    }

    template<class FuncType>
    inline constexpr auto sin(FuncType&& func){
            return jf::var::sin(func);
    }

    template<class FuncType>
    inline constexpr auto cos(FuncType&& func){
        return jf::var::cos(func);
    }

} // namespace jf::math
