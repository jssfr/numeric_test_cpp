module;

#ifdef _MSC_VER
    #define __FUNCTION_NAME__ __FUNCSIG__
    #define USING_MSVC 1
#else
    #define __FUNCTION_NAME__ __PRETTY_FUNCTION__
    #define USING_MSVC 0
#endif

#include<string>
#include<cstring>

#include<typeinfo>

export module Config:Name;
import :Types;

template<typename Type>
std::string type_to_string()
{
    std::string fname(__FUNCTION_NAME__);//name of function with the type value
    
    const char* ftext0 = "jf::types::st_type_list";
    if(USING_MSVC){
        const char* ftext1 = "type_to_string<";
        
        if(auto pos = fname.find(ftext0); pos != std::string::npos){
            pos += std::strlen(ftext0);
            fname = fname.substr(pos+1);
            if(pos = fname.find("void"); pos != std::string::npos){
               fname = fname.substr(0, fname.size()-8);
            }
        }
        else if(auto pos = fname.find(ftext1); pos != std::string::npos){
            pos += std::strlen(ftext1);
            fname = fname.substr(pos);
            if(pos = fname.find("void"); pos != std::string::npos){
                 fname = fname.substr(0, fname.size()-8);
            }
        }
    }
    else{
        const char* ftext1 = "Type = ";
        if(auto pos = fname.find(ftext0); pos != std::string::npos){
            pos += std::strlen(ftext0);
            fname = fname.substr(pos);
            fname = fname.substr(0, fname.size()-2);
        }
        else if(auto pos = fname.find(ftext1); pos != std::string::npos){
            pos += std::strlen(ftext1);
            fname = fname.substr(pos+1);
            fname = fname.substr(0, fname.size()-1);
        }

    }
    return fname;// returns only the type
}

export{
    template<typename type_arg>
    std::string GetTypeName(){
        return type_to_string<type_arg>();
    }

    template<typename instance_arg>
    std::string GetTypeCategory(instance_arg&&){
        return type_to_string<instance_arg>();
    }

}

namespace hidden
    {
    
    template <size_t StartIndex, size_t EndIndex>
    struct st_static_loop 
    {
            //// \brief get all type in list
        template <typename Type, typename... Types>
        static void type_list(std::string& out, const jf::types::type_list_t<Type, Types...>&) {
            if constexpr (StartIndex < EndIndex) {
                out += GetTypeName<Type>();
                out += (StartIndex + 1 != EndIndex ? ", " : "");
                st_static_loop<StartIndex + 1, EndIndex>::type_list(out,
                                                                        jf::types::type_list_t<Types...>{});
            }
        }
    
        static void type_list(...) {}  // fallback
    };
    template<typename... types>
        struct st_type_id;
    template<>
        struct st_type_id<>
        {
            static std::string type(){
                return "";
            }
        };

    template<typename T, typename... Types>
        struct st_type_id<T, Types...> {
            static std::string type() {
                    return GetTypeName<T>() + (sizeof...(Types) == 0? "": ", " + 
                                                                st_type_id<Types...>::type());
            }
        };
    template<typename... types>
        struct st_type_id<jf::types::type_list_t<types...> >
        {
            static std::string type(){
                return st_type_id<types...>::type();
            }
        };

    template <typename Type, typename... Types>
    std::string type_func() {
        std::string out = "<";
    
        if constexpr (jf::types::is_template_v<Type> && sizeof...(Types) == 0) {
            using ty = jf::types::template_t<Type>;
            std::string helper = GetTypeName<ty>();
            out += helper;
        } else {
            out += hidden::st_type_id<jf::types::type_list_t<Type, Types...>>::type();
        }
        out += ">";
    
        return out;
    }; 
}  // namespace hidden


export namespace jf::types{
    //// @brief hold types and convert to a string
//// @return list of types
template <typename Type, typename... Types>
std::string type_list_v(){
    return ::hidden::type_func<Type, Types...>();  // -----
    }
} // namespace jf::types
 

