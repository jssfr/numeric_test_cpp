module;

#ifdef _MSC_VER
    #define __FUNCTION_NAME__ __FUNCSIG__
#else
    #define __FUNCTION_NAME__ __PRETTY_FUNCTION__
#endif

#include<string>
#include<cstring>

export module Config:Name;

template<typename Type>
std::string type_to_string()
{
    std::string fname(__FUNCTION_NAME__);//name of function with the type value
    const char* ftext = "[with Type = ";
    auto pos = fname.find(ftext) + std::strlen(ftext);
    fname = fname.substr(pos);
    pos = fname.find_first_of(';');
    return fname.substr(0, pos);// returns only the type
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
 

