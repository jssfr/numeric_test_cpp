#include <vector>
#include <string>
#include <iostream>
import Config;
//import <iostream>;
//import std;

int main()
{
    std::cout<<"Teste Teste\nisso aqui e uma porcaria!!\n\n";
   int a = 5;
   std::cout<<GetTypeCategory(a)<<"\n";
////    
   using type_d = double;
   std::cout<<GetTypeName<type_d>()<<"\n";
   
    using typelist = jf::type_list_t<int, double, short>;
    std::cout<<GetTypeName<typelist>()<<"\n";
  // std::cout<<\n", sizeof(type_d));
    std::cout<<jf::type_list_v<int, double, short, bool>() <<"\n";
    std::cout<<jf::type_list_v<typelist>() <<"\n";
    std::cout<<jf::type_list_v<std::vector<int>, int, short>() <<"\n";
    using typelist2 = jf::type_list_t<int, double, short, std::string, char, bool>;
    std::cout<< jf::type_list_v<typelist2>() << "\n";
    return 0 ;
}
// some modification is needed to port to module, but  it works fine.

