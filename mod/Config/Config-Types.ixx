module;

#include<string>
#include<typeinfo>

#ifndef _MSC_VER
#define USING_MSVC 0
#else
#define USING_MSVC 1
#endif

#include "std_input.hpp"

export module Config:Types;

import :Name;


export namespace jf::types
{
    namespace hidden {
    
    // primary template
    template <typename T>
    struct st_remove_reference {
        using type = T;
    };

    // especializations

    template <typename T>
    struct st_remove_reference<T&> {
        using type = T;
    };

    template <typename T>
    struct st_remove_reference<T&&> {
        using type = T;
    };

    }  // end namespace hidden

    template <typename T>
    using remove_reference_t =
        hidden::st_remove_reference<T>::type;  // same as std::remove_reference_t< T >

    template <typename T>
    using remove_const_reference_t = std::remove_const_t<std::remove_reference_t<T>>;

    template <typename Type>
    using remove_cv_ref_t = std::remove_cv_t<std::remove_reference_t<Type>>;

    struct NoType {};  ///// \brief if no type matches, this struct type is used
                       

    namespace hidden {
    // struct or class that doesn't have a body definition is called incomplete type
    template <typename... Types>
    struct st_type_count;  // primary template;

    //especialization
    template <>             // doesnt have any template parameter
    struct st_type_count<>  // <> specialization
    {
        static constexpr size_t value = 0;
    };

    template <typename Type,
              typename... Types>          // parameter count can be diferent from primary template one
    struct st_type_count<Type, Types...>  // type arguments mean one or more
    {
        static constexpr size_t value = 1 + st_type_count<Types...>::value;  // recursive
    };

    template <template <typename...> class TmpType, typename... Types>
    struct st_type_count<TmpType<Types...>> {
        static constexpr size_t value = st_type_count<Types...>::value;
    };
    }  // end namespace hidden

    template <typename... Types>
    constexpr auto type_count_v = hidden::st_type_count<Types...>::value;

    //============================

    template <typename... Types>
    struct st_type_list {};
    
    template <typename... Types>
    using type_list_t = st_type_list<Types...>;  //// \brief hold types

// end jf
    namespace hidden {
     template <typename T>
     struct st_is_template  // primary
     {
         using type = T;
         static constexpr bool value = false;
     };
     
     template <template <typename...> class TmpType, typename... Types>
     struct st_is_template<
         TmpType<Types...>>  // especialization , if argument is a template container this is called
     {
         using type = TmpType<Types...>;
         static constexpr bool value = true;
     };
    
    }  // end namespace hidden
    
    template <typename T>
    constexpr auto is_template_v = hidden::st_is_template<T>::value;  //// \return true or false
    
    template <typename T>
    using template_t = hidden::st_is_template<T>::type;  //// \return type or a template type

 
    namespace hidden {  //       element_type
    template <typename T>
    struct st_element_type {
        using type = T;
    };

    template <template <typename, typename...> class TmpType, typename Type, typename... Types>
    struct st_element_type<TmpType<Type, Types...>> {
        using type = Type;
    };
    }  // end namespace hidden

    template <typename T>
    using element_type_t = hidden::st_element_type<T>::type;

    namespace hidden {  // is std container
    // primary template
    template <typename T>
    struct st_is_container {
        static constexpr bool value = false;
    };

    // especializations

    template <typename Type, typename... Types>
    struct st_is_container<std::vector<Type, Types...>> {
        static constexpr bool value = true;
    };

    template <typename Type, typename... Types>
    struct st_is_container<std::deque<Type, Types...>> {
        static constexpr bool value = true;
    };

    template <typename Type, typename... Types>
    struct st_is_container<std::set<Type, Types...>> {
        static constexpr bool value = true;
    };

    template <class keyType, typename ValueType, typename... Types>
    struct st_is_container<std::map<keyType, ValueType, Types...>> {
        static constexpr bool value = true;
    };

    // TODO : add more containers
    }  // namespace hidden

    template <typename T>
    constexpr auto is_container_v = hidden::st_is_container<T>::value;

    //////////////  get<> implementation
    template <auto Index, typename CntrType>
    decltype(auto) get(CntrType&& cntr) {
        return std::get<Index>(std::forward<CntrType>(cntr));
    }

    template <auto Index, auto... Indices, typename CntrType>
    decltype(auto) get(CntrType&& cntr) {
        if constexpr (sizeof...(Indices) > 0)
            return get<Indices...>(get<Index>(std::forward<CntrType>(cntr)));
        else
            return get<Index>(std::forward<CntrType>(cntr));
    }

    template <typename IndexType, auto Index, auto... Indices, typename CntrType>
    decltype(auto) get(CntrType&& cntr, IndexType index) {
        if constexpr (sizeof...(Indices) > 0)
            return get<Indices...>(get<Index>(std::forward<CntrType>(cntr)))[(size_t)index];
        else
            return get<Index>(std::forward<CntrType>(cntr))[(size_t)index];
    }

    // -----

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

    //template<typename t>
    //    struct st_type_id<t>
    //    {
    //        static constexpr std::string type = std::string(typeid(t).name());
    //    };
    template<typename T, typename... Types>
        struct st_type_id<T, Types...> {
            static std::string type() {
                if(USING_MSVC){
                    return std::string(typeid(T).name()) + (sizeof...(Types) == 0? "": ", " + st_type_id<Types...>::type());
                }
                else{ 
                    return GetTypeName<T>() + (sizeof...(Types) == 0? "": ", " + 
                                                                st_type_id<Types...>::type());
                }       
            }
        };
    template<typename... types>
        struct st_type_id<jf::types::type_list_t<types...> >
        {
            static std::string type(){
                return st_type_id<types...>::type();
            }
        };
    }// end hidden

     namespace hidden {
    template <typename Type, typename... Types>
    std::string type_func() {
        std::string out = "<";
    
        if constexpr (is_template_v<Type> && sizeof...(Types) == 0) {
            using ty = template_t<Type>;
            std::string helper = GetTypeName<ty>().substr(16);
            // helper = helper.substr(1, helper.size()-2);
            out += helper.substr(1, helper.size() - 2);
            // out += helper;
        } else {
            //hidden::st_static_loop<0, sizeof...(Types) + 1>::template type_list(
            //    out, jf::types::type_list_t<Type, Types...>{});
            out += hidden::st_type_id<jf::types::type_list_t<Type, Types...>>::type();
        }
        out += ">";
    
        return out;
    };  // just  print types
    }  // namespace hidden



    //// @brief hold types and convert to a string
//// @return list of types
template <typename Type, typename... Types>
std::string type_list_v(){
   return hidden::type_func<Type, Types...>();  // -----
    
    }


    namespace hidden {
    template <typename S, typename T>
    auto common_type_func(S&& s, T&& t) -> remove_const_reference_t<decltype(true ? s : t)>;
    NoType common_type_func(...);

    template <typename S, typename T>
    using common_type_of_s_and_t = decltype(common_type_func(std::declval<S>(), std::declval<T>()));

    }  // namespace hidden

    template <typename S, typename T>
    constexpr auto exist_common_type_v = !std::is_same_v<NoType, hidden::common_type_of_s_and_t<S, T>>;

    template <typename S, typename T>
    constexpr auto are_types_operable_v = exist_common_type_v<S, T>;

    /// @brief A check to figure out if a type is integer
    /// @param T
    /// @return true or false
    template <typename T>
    constexpr auto is_integer_v =
        !std::is_same_v<bool, T> && std::is_integral_v<remove_const_reference_t<T>> &&
        std::is_arithmetic_v<remove_const_reference_t<T>>;

    namespace hidden {
    template <typename T, typename... Types>
    struct st_common_type;

    template <typename T>
    struct st_common_type<T> {
        using type = T;
    };

    template <typename S, typename T>
    struct st_common_type<S, T> {
        using type = common_type_of_s_and_t<S, T>;
    };

    template <typename S, typename T, typename... Types>
    struct st_common_type<S, T, Types...> {
        using type = common_type_of_s_and_t<typename st_common_type<S, T>::type,
                                            typename st_common_type<T, Types...>::type>;
    };
    }  // namespace hidden

    /// @brief Computes a common type between 1 or more types
    /// @tparam T
    /// @tparam ...Types
    /// @return A common type or Notype
    /// @
    template <typename... Types>
    using common_type_t = hidden::st_common_type<Types...>::type;

    /// @brief type check to see if exist a common type
    /// @return true or false
    /// @tparam T
    /// @tparam ...Types
    template <typename T, typename... Types>
    constexpr auto common_type_v = !std::is_same_v<NoType, common_type_t<T, Types...>>;

    template <typename T, typename... Types>
    concept CommonTypeExists_c = common_type_v<T, Types...>;

    namespace hidden {
    template <typename T, bool bInteger = false>
    struct st_make_signed {
        using type = T;
    };

    template <typename T>
    struct st_make_signed<T, true> {
        using type = std::make_signed_t<T>;
    };
    }  // namespace hidden

    template <typename T, bool bInteger>
    using signed_type_t = hidden::st_make_signed<T, bInteger>::type;
    
    template<typename T>
    using make_signed_t = signed_type_t<T, true>;
    
    template<typename... types>
    using common_signed_t = make_signed_t<common_type_t<std::remove_cvref_t<types>...>>;
    

    template<typename T>
    concept signed_c = std::is_signed_v<std::remove_cvref_t<T> >;

    template<typename T>
    concept unsigned_c = !signed_c<T>;

    template<typename T>
    concept arithmetic_c = std::is_arithmetic_v< std::remove_cvref_t<T> >;

    ///////////////////////// FIND TYPE

    namespace hidden {
    template <template <typename, typename> class BinaryOpr, typename ArgType, typename Type,
              typename... Types>
    struct st_find_type {
        static constexpr bool value = st_find_type<BinaryOpr, ArgType, Type>::value
                                          ? true
                                          : st_find_type<BinaryOpr, ArgType, Types...>::value;
        using type = std::conditional_t<st_find_type<BinaryOpr, ArgType, Type>::value, Type,
                                        typename st_find_type<BinaryOpr, ArgType, Types...>::type>;
    };

    template <template <typename, typename> class BinaryOpr, typename ArgType, typename Type>
    struct st_find_type<BinaryOpr, ArgType, Type> {
        static constexpr bool value = BinaryOpr<Type, ArgType>::value;
        using type = std::conditional_t<value, Type, NoType>;
    };

    template <template <typename, typename> class BinaryOpr, typename ArgType, typename Type,
              typename... Types>
    struct st_find_type<BinaryOpr, ArgType, jf::types::type_list_t<Type, Types...>> {
        static constexpr bool value = st_find_type<BinaryOpr, ArgType, Type, Types...>::value;
        using type = st_find_type<BinaryOpr, ArgType, Type, Types...>::type;
    };

    template <template <typename, typename> class BinaryOpr, typename ArgType, typename Type>
    struct st_find_type<BinaryOpr, ArgType, jf::types::type_list_t<Type>> {
        static constexpr bool value = st_find_type<BinaryOpr, ArgType, Type>::value;
        using type = st_find_type<BinaryOpr, ArgType, Type>::type;
    };

    }  // namespace hidden

    template <typename ArgType, typename Type, typename... Types>
    constexpr auto same_type_v = hidden::st_find_type<std::is_same, ArgType, Type, Types...>::value;
    
    template <typename ArgType, typename Type, typename... Types>
    concept same_type_c = same_type_v<ArgType, Type, Types...>;
    
    template <typename ArgType, typename Type, typename... Types>
    using same_type_t = hidden::st_find_type<std::is_same, ArgType, Type, Types...>::type;

    template <typename ArgType, typename Type, typename... Types>
    constexpr auto constructible_v =
        hidden::st_find_type<std::is_constructible, ArgType, Type, Types...>::value;
    
    template <typename ArgType, typename Type, typename... Types>
    constexpr auto constructible_c = constructible_v<ArgType, Type, Types...>;
    
    template <typename ArgType, typename Type, typename... Types>
    using constructible_t = hidden::st_find_type<std::is_constructible, ArgType, Type, Types...>::type;

    //------------------
    namespace hidden {
    template <typename ArgType, typename Type, typename... Types>
    struct st_find_best_type {
        static constexpr bool value =
            same_type_v<ArgType, Type, Types...> ? true : constructible_v<ArgType, Type, Types...>;
        using type = std::conditional_t<same_type_v<ArgType, Type, Types...>,
                                        same_type_t<ArgType, Type, Types...>,
                                        constructible_t<ArgType, Type, Types...>>;
    };
    }  // namespace hidden
    template <typename ArgType, typename Type, typename... Types>
    constexpr auto best_type_v = hidden::st_find_best_type<ArgType, Type, Types...>::value;

    template <typename ArgType, typename Type, typename... Types>
    using best_type_t = hidden::st_find_best_type<ArgType, Type, Types...>::type;

    ////////////////////////

    namespace hidden {
    template <typename... Types>
    struct st_select_first_type;

    template <>
    struct st_select_first_type<> {
        using type = NoType;
        using typelist = jf::types::type_list_t<>;
    };

    template <typename Type>
    struct st_select_first_type<Type> {
        using type = Type;
        using typelist = jf::types::type_list_t<>;
    };

    template <typename Type, typename... Types>
    struct st_select_first_type<Type, Types...> {
        using type = Type;
        using typelist = jf::types::type_list_t<Type>;
    };
    //////////////////

    template <>
    struct st_select_first_type<jf::types::type_list_t<>> {
        using type = NoType;
        using typelist = jf::types::type_list_t<>;
    };

    template <typename Type>
    struct st_select_first_type<jf::types::type_list_t<Type>> {
        using type = Type;
        using typelist = jf::types::type_list_t<Type>;
    };

    template <typename Type, typename... Types>
    struct st_select_first_type<jf::types::type_list_t<Type, Types...>> {
        using type = Type;
        using typelist = jf::types::type_list_t<Type>;
    };

    template<auto Ele1, auto ...Eles>
    struct st_select_first_element;

    template<auto Ele1>
    struct st_select_first_element<Ele1>{
        static constexpr auto value = Ele1;
    };

    template<auto Ele1, auto Ele2>
    struct st_select_first_element<Ele1, Ele2>{
        static constexpr auto value = Ele1;
    };

    template<auto Ele1, auto Ele2, auto ...Eles>
    struct st_select_first_element<Ele1, Ele2, Eles...>{
        static constexpr auto value = Ele1;
    };

    }  // namespace hidden

    template <typename... Types>
    using select_first_type_t = hidden::st_select_first_type<Types...>::type;
    
    template<auto Ele1, auto ... Eles>
    constexpr auto select_first_element_v = hidden::st_select_first_element<Ele1, Eles...>::value;
    
    template <typename... Types>
    using select_first_type_list_t = hidden::st_select_first_type<Types...>::typelist;

    namespace hidden {
    template <typename... Types>
    struct st_select_last_type;

    template <>
    struct st_select_last_type<> {
        using type = jf::types::NoType;
        using typelist = jf::types::type_list_t<>;
    };

    template <typename Type>
    struct st_select_last_type<Type> {
        using type = Type;
        using typelist = jf::types::type_list_t<type>;
    };

    template <typename Type, typename... Types>
    struct st_select_last_type<Type, Types...> {
        using type = std::conditional_t<sizeof...(Types) == 0, Type,
                                        typename st_select_last_type<Types...>::type>;
        using typelist = jf::types::type_list_t<type>;
    };
    //////////////////

    template <>
    struct st_select_last_type<jf::types::type_list_t<>> {
        using type = jf::types::NoType;
        using typelist = jf::types::type_list_t<>;
    };

    template <typename Type>
    struct st_select_last_type<jf::types::type_list_t<Type>> {
        using type = Type;
        using typelist = jf::types::type_list_t<type>;
    };

    template <typename Type, typename... Types>
    struct st_select_last_type<jf::types::type_list_t<Type, Types...>> {
        using type = st_select_last_type<Type, Types...>::type;
        using typelist = jf::types::type_list_t<type>;
    };

    template <auto... Eles>
    struct st_select_last_element;

    template <auto Ele>
    struct st_select_last_element<Ele> {
        static constexpr auto value = Ele;
    };

    template <auto Ele, auto... Eles>
    struct st_select_last_element<Ele, Eles...> {
        static constexpr auto value = (sizeof...(Eles) == 0)? Ele :
                                        st_select_last_element<Eles...>::value;
    };

    }  // namespace hidden

    template <typename... Types>
    using select_last_type_t = hidden::st_select_last_type<Types...>::type;

    template <auto... Eles>
    constexpr auto select_last_element_v = hidden::st_select_last_element<Eles...>::value;

    template <typename... Types>
    using select_last_type_list_t = hidden::st_select_last_type<Types...>::typelist;
// ---------------------------------------


    namespace hidden {
    template <auto SelectIndex, typename... Types>
    struct st_select_nth_type;

    template <auto SelectIndex>
    struct st_select_nth_type<SelectIndex> {
        using type = jf::types::NoType;
        using typelist = jf::types::type_list_t<>;
    };

    template <auto SelectIndex, typename Type, typename... Types>
    struct st_select_nth_type<SelectIndex, Type, Types...> {
        using type = std::conditional_t<SelectIndex == 0, Type,
                                        typename st_select_nth_type<SelectIndex - 1, Types...>::type>;
        using typelist = jf::types::type_list_t<type>;
    };

    /////////

    template <auto SelectIndex>
    struct st_select_nth_type<SelectIndex, jf::types::type_list_t<>> {
        using type = jf::types::NoType;
        using typelist = jf::types::type_list_t<>;
    };

    template <typename Type, typename... Types>
    struct st_select_nth_type<0, jf::types::type_list_t<Type, Types...>> {
        using type = Type;
        using typelist = jf::types::type_list_t<type>;
    };

    template <auto SelectIndex, typename Type, typename... Types>
    struct st_select_nth_type<SelectIndex, jf::types::type_list_t<Type, Types...>> {
        using type = st_select_nth_type<SelectIndex, Type, Types...>::type;
        using typelist = jf::types::type_list_t<type>;
    };

    template <auto SelectIndex, auto... Eles>
    struct st_select_nth_element;

    template <auto SelectIndex, auto Ele, auto... Eles>
    struct st_select_nth_element<SelectIndex, Ele, Eles...> {
        static constexpr auto value = (SelectIndex == 0)? Ele :
                                         st_select_nth_element<SelectIndex - 1, Eles...>::value;
    };
    // template <auto SelectIndex, auto Ele, auto... Eles>
    // struct st_select_nth_element<SelectIndex, Ele, std::vector{Eles...}> {
    //     constexpr auto value = (SelectIndex == 0)? Ele :
    //                                      st_select_nth_element<SelectIndex - 1, Eles...>::value;
    // };
    }  // namespace hidden

    template <std::size_t SelectIndex, typename... Types>
        requires (SelectIndex <= sizeof...(Types))
    using select_nth_type_t = hidden::st_select_nth_type<SelectIndex, Types...>::type;
    
    template <std::size_t SelectIndex, auto... Eles>
        requires (SelectIndex <= sizeof...(Eles))
    constexpr auto select_nth_element_v = hidden::st_select_nth_element<SelectIndex, Eles...>::value;

    template <auto SelectIndex, typename... Types>
    using select_nth_type_list_t = hidden::st_select_nth_type<SelectIndex, Types...>::typelist;

    namespace hidden {
    template <typename ArgType, typename ListType>
    struct st_push_front_type;

    template <typename ArgType, typename... Types>
    struct st_push_front_type<ArgType, jf::types::type_list_t<Types...>> {
        using type = jf::types::type_list_t<ArgType, Types...>;
    };

    template <typename ArgType, typename ListType>
    using push_front_type_t = st_push_front_type<ArgType, ListType>::type;

    template <typename ArgType, typename ListType>
    struct st_push_back_type;

    template <typename ArgType, typename... Types>
    struct st_push_back_type<ArgType, jf::types::type_list_t<Types...>> {
        using type = jf::types::type_list_t<Types..., ArgType>;
    };

    template <typename ArgType, typename ListType>
    using push_back_type_t = st_push_back_type<ArgType, ListType>::type;

    template <typename ListType>
    struct st_pop_front_type;

    template <>
    struct st_pop_front_type<jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<>;
    };

    template <typename ArgType, typename... Types>
    struct st_pop_front_type<jf::types::type_list_t<ArgType, Types...>> {
        using type = jf::types::type_list_t<Types...>;
    };

    template <typename ListType>
    using pop_front_type_t = st_pop_front_type<ListType>::type;

    template <typename LeftList, typename RightList>
    struct st_pop_back_type;

    template <typename... LeftTypes>
    struct st_pop_back_type<jf::types::type_list_t<LeftTypes...>, jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<LeftTypes...>;
    };

    template <typename... LeftTypes, typename Type, typename... RightTypes>
    struct st_pop_back_type<jf::types::type_list_t<LeftTypes...>, jf::types::type_list_t<Type, RightTypes...>> {
        using leftlist = push_back_type_t<Type, jf::types::type_list_t<LeftTypes...>>;
        using rightlist = jf::types::type_list_t<RightTypes...>;

        using type = std::conditional_t<sizeof...(RightTypes) == 0, jf::types::type_list_t<LeftTypes...>,
                                        typename st_pop_back_type<leftlist, rightlist>::type>;
    };

    template <typename ListType>
    using pop_back_type_t = st_pop_back_type<jf::types::type_list_t<>, ListType>::type;

    }  // namespace hidden

    template <typename ArgType, typename ListType>
    using push_front_type_t = hidden::push_front_type_t<ArgType, ListType>;

    template <typename ArgType, typename ListType>
    using push_back_type_t = hidden::push_back_type_t<ArgType, ListType>;

    template <typename ListType>
    using pop_front_type_t = hidden::pop_front_type_t<ListType>;

    template <typename ListType>
    using pop_back_type_t = hidden::pop_back_type_t<ListType>;

    namespace hidden {
    template <typename ArgType, typename... Types>
    struct st_is_type_in_list;

    template <typename ArgType>
    struct st_is_type_in_list<ArgType> {
        static constexpr bool value = false;
    };

    template <typename ArgType, typename Type, typename... Types>
    struct st_is_type_in_list<ArgType, Type, Types...> {
        static constexpr bool value =
            std::is_same_v<ArgType, Type> ? true : st_is_type_in_list<ArgType, Types...>::value;
    };

    template <typename ArgType>
    struct st_is_type_in_list<ArgType, jf::types::type_list_t<>> {
        static constexpr bool value = false;
    };

    template <typename ArgType, typename Type, typename... Types>
    struct st_is_type_in_list<ArgType, jf::types::type_list_t<Type, Types...>> {
        static constexpr bool value = st_is_type_in_list<ArgType, Type, Types...>::value;
    };

    template <typename ArgType, typename... Types>
    constexpr auto is_type_in_list_v = st_is_type_in_list<ArgType, Types...>::value;

    template <typename ListType, typename... Types>
    struct st_unique_type;

    template <typename... Types>
    struct st_unique_type<jf::types::type_list_t<Types...>> {
        using type = jf::types::type_list_t<Types...>;
    };

    template <typename... Types, typename Type>
    struct st_unique_type<jf::types::type_list_t<Types...>, Type> {
        using list = jf::types::type_list_t<Types...>;
        using type =
            std::conditional_t<is_type_in_list_v<Type, list>, list, jf::types::push_back_type_t<Type, list>>;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_unique_type<jf::types::type_list_t<Types...>, Type, RightTypes...> {
        using list = st_unique_type<jf::types::type_list_t<Types...>, Type>::type;
        using type = st_unique_type<list, RightTypes...>::type;
    };

    template <typename... Types, typename... RightTypes>
    struct st_unique_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<RightTypes...>> {
        using type = st_unique_type<jf::types::type_list_t<Types...>, RightTypes...>::type;
    };

    template <typename... Types>
    using unique_type_t = st_unique_type<jf::types::type_list_t<>, Types...>::type;

    template <typename LeftList, typename RightList>
    struct st_prepend_type;

    template <typename... LeftTypes, typename... RightTypes>
    struct st_prepend_type<jf::types::type_list_t<LeftTypes...>, jf::types::type_list_t<RightTypes...>> {
        using type = jf::types::type_list_t<LeftTypes..., RightTypes...>;
    };

    template <typename LeftList, typename RightList>
    struct st_append_type;

    template <typename... LeftTypes, typename... RightTypes>
    struct st_append_type<jf::types::type_list_t<LeftTypes...>, jf::types::type_list_t<RightTypes...>> {
        using type = jf::types::type_list_t<RightTypes..., LeftTypes...>;
    };
    }  // namespace hidden

    template <typename ArgType, typename... Types>
    constexpr auto is_type_in_list_v = hidden::is_type_in_list_v<ArgType, Types...>;

    template <typename... Types>
    using unique_type_t = hidden::unique_type_t<Types...>;

    template <typename LeftList, typename RightList>
    using prepend_type_t = hidden::st_prepend_type<LeftList, RightList>::type;

    template <typename LeftList, typename RightList>
    using append_type_t = hidden::st_append_type<LeftList, RightList>::type;

    namespace hidden {
    template <typename LeftList, typename... Types>
    struct st_union_type;

    template <typename... Types>
    struct st_union_type<jf::types::type_list_t<Types...>> {
        using type = jf::types::type_list_t<Types...>;
    };
    template <typename... Types, typename Type>
    struct st_union_type<jf::types::type_list_t<Types...>, Type> {
        using list = jf::types::type_list_t<Types...>;
        using type =
            std::conditional_t<is_type_in_list_v<Type, list>, list, push_back_type_t<Type, list>>;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_union_type<jf::types::type_list_t<Types...>, Type, RightTypes...> {
        using list = st_union_type<jf::types::type_list_t<Types...>, Type>::type;
        using type = std::conditional_t<sizeof...(RightTypes) == 0, list,
                                        typename st_union_type<list, RightTypes...>::type>;
    };
    //////////////////

    template <typename... Types>
    struct st_union_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<Types...>;
    };
    template <typename... Types, typename Type>
    struct st_union_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<Type>> {
        using type = st_union_type<jf::types::type_list_t<Types...>, Type>::type;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_union_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<Type, RightTypes...>> {
        using type = st_union_type<jf::types::type_list_t<Types...>, Type, RightTypes...>::type;
    };

    template <typename LeftList, typename RightList>
    using union_type_t = st_union_type<LeftList, RightList>::type;

    }  // namespace hidden

    template <typename LeftList, typename RightList>
    using union_type_t = hidden::union_type_t<unique_type_t<LeftList>, unique_type_t<RightList>>;

    namespace hidden {
    template <typename LeftList, typename... Types>
    struct st_intersection_type;

    template <typename... Types>
    struct st_intersection_type<jf::types::type_list_t<Types...>> {
        using type = jf::types::type_list_t<Types...>;
    };

    template <typename... Types, typename Type>
    struct st_intersection_type<jf::types::type_list_t<Types...>, Type> {
        using list = jf::types::type_list_t<Types...>;
        using type = std::conditional_t<jf::types::is_type_in_list_v<Type, list>, jf::types::type_list_t<Type>,
                                        jf::types::type_list_t<>>;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_intersection_type<jf::types::type_list_t<Types...>, Type, RightTypes...> {
        using leftlist = jf::types::type_list_t<Types...>;
        using list = st_intersection_type<leftlist, Type>::type;

        using type = std::conditional_t<
            sizeof...(RightTypes) == 0, list,
            jf::types::prepend_type_t<list, typename st_intersection_type<leftlist, RightTypes...>::type>>;
    };
    /////////////////

    template <typename... Types>
    struct st_intersection_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<Types...>;
    };

    template <typename... Types, typename Type>
    struct st_intersection_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<Type>> {
        using type = st_intersection_type<jf::types::type_list_t<Types...>, Type>::type;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_intersection_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<Type, RightTypes...>> {
        using type = st_intersection_type<jf::types::type_list_t<Types...>, Type, RightTypes...>::type;
    };

    template <typename LeftList, typename RightList>
    using intersection_type_t = st_intersection_type<LeftList, RightList>::type;

    }  // namespace hidden

    template <typename LeftList, typename RightList>
    using intersection_type_t =
        hidden::intersection_type_t<unique_type_t<LeftList>, unique_type_t<RightList>>;

    namespace hidden {
    template <typename ArgType, typename... Types>
    struct st_remove_type;

    template <typename ArgType>
    struct st_remove_type<ArgType> {
        using type = jf::types::type_list_t<>;
    };

    template <typename ArgType, typename Type>
    struct st_remove_type<ArgType, Type> {
        using type =
            std::conditional_t<std::is_same_v<ArgType, Type>, jf::types::type_list_t<>, jf::types::type_list_t<Type>>;
    };

    template <typename ArgType, typename Type, typename... RightTypes>
    struct st_remove_type<ArgType, Type, RightTypes...> {
        using list = st_remove_type<ArgType, Type>::type;
        using type = std::conditional_t<
            sizeof...(RightTypes) == 0, list,
            prepend_type_t<list, typename st_remove_type<ArgType, RightTypes...>::type>>;
    };
    /////////////////////////

    template <typename ArgType>
    struct st_remove_type<ArgType, jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<>;
    };

    template <typename ArgType, typename Type>
    struct st_remove_type<ArgType, jf::types::type_list_t<Type>> {
        using type = st_remove_type<ArgType, Type>::type;
    };

    template <typename ArgType, typename Type, typename... RightTypes>
    struct st_remove_type<ArgType, jf::types::type_list_t<Type, RightTypes...>> {
        using type = st_remove_type<ArgType, Type, RightTypes...>::type;
    };

    template <typename ArgType, typename TypeList>
    using remove_type_t = st_remove_type<ArgType, TypeList>::type;

    }  // namespace hidden

    template <typename ArgType, typename TypeList>
    using remove_type_t = hidden::remove_type_t<ArgType, TypeList>;

    namespace hidden {
    template <typename LeftList, typename... Types>
    struct st_difference_type;

    template <typename... Types>
    struct st_difference_type<jf::types::type_list_t<Types...>> {
        using type = jf::types::type_list_t<Types...>;
    };

    template <typename... Types, typename Type>
    struct st_difference_type<jf::types::type_list_t<Types...>, Type> {
        using type = remove_type_t<Type, jf::types::type_list_t<Types...>>;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_difference_type<jf::types::type_list_t<Types...>, Type, RightTypes...> {
        using list = st_difference_type<jf::types::type_list_t<Types...>, Type>::type;
        using type = std::conditional_t<sizeof...(RightTypes) == 0, list,
                                        typename st_difference_type<list, RightTypes...>::type>;
    };

    /////////////////////////////////

    template <typename... Types>
    struct st_difference_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<Types...>;
    };

    template <typename... Types, typename Type>
    struct st_difference_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<Type>> {
        using type = st_difference_type<jf::types::type_list_t<Types...>, Type>::type;
    };

    template <typename... Types, typename Type, typename... RightTypes>
    struct st_difference_type<jf::types::type_list_t<Types...>, jf::types::type_list_t<Type, RightTypes...>> {
        using type = st_difference_type<jf::types::type_list_t<Types...>, Type, RightTypes...>::type;
    };

    template <typename LeftList, typename RightList>
    using difference_type_t = st_difference_type<LeftList, RightList>::type;

    }  // namespace hidden

    /// @return types that are not in the RightList
    /// @tparam LeftList
    /// @tparam RightList
    template <typename LeftList, typename RightList>
    using difference_type_t = hidden::difference_type_t<LeftList, RightList>;

    namespace hidden {
    template <typename NewType, typename OldType, typename... Types>
    struct st_replace_type;

    template <typename NewType, typename OldType>
    struct st_replace_type<NewType, OldType> {
        using type = jf::types::type_list_t<>;
    };

    template <typename NewType, typename OldType, typename Type>
    struct st_replace_type<NewType, OldType, Type> {
        using type = std::conditional_t<std::is_same_v<OldType, Type>, jf::types::type_list_t<NewType>,
                                        jf::types::type_list_t<Type>>;
    };

    template <typename NewType, typename OldType, typename Type, typename... RightTypes>
    struct st_replace_type<NewType, OldType, Type, RightTypes...> {
        using list = st_replace_type<NewType, OldType, Type>::type;
        using type = std::conditional_t<
            sizeof...(RightTypes) == 0, list,
            prepend_type_t<list, typename st_replace_type<NewType, OldType, RightTypes...>::type>>;
    };

    /////////////////////

    template <typename NewType, typename OldType>
    struct st_replace_type<NewType, OldType, jf::types::type_list_t<>> {
        using type = jf::types::type_list_t<>;
    };

    template <typename NewType, typename OldType, typename Type>
    struct st_replace_type<NewType, OldType, jf::types::type_list_t<Type>> {
        using type = st_replace_type<NewType, OldType, Type>::type;
    };

    template <typename NewType, typename OldType, typename Type, typename... RightTypes>
    struct st_replace_type<NewType, OldType, jf::types::type_list_t<Type, RightTypes...>> {
        using type = st_replace_type<NewType, OldType, Type, RightTypes...>::type;
    };

    template <typename NewType, typename OldType, typename TypeList>
    using replace_type_t = st_replace_type<NewType, OldType, TypeList>::type;

    }  // namespace hidden

    /// @brief Replace all ocurrences of OldType in list with NewType
    /// @tparam NewType
    /// @tparam OldType
    /// @tparam TypeList
    template <typename NewType, typename OldType, typename TypeList>
    using replace_type_t = hidden::replace_type_t<NewType, OldType, TypeList>;

    template <typename T, T... Ints>
    struct integer_sequence {
        static constexpr std::size_t size() { return sizeof...(Ints); }
    };

    namespace hidden {
    template <typename T, typename... Ts>
    struct st_add_sequence;

    template <typename T, T... Ints1, T... Ints2>
    struct st_add_sequence<integer_sequence<T, Ints1...>, integer_sequence<T, Ints2...>> {
        using type = integer_sequence<T, Ints1..., Ints2...>;
    };

    template <typename T, T... Ints>
    struct st_add_sequence<integer_sequence<T, Ints...>> {
        using type = integer_sequence<T, Ints...>;
    };

    template <typename SeqType1, typename SeqType2>
    using add_sequence_t = st_add_sequence<SeqType1, SeqType2>::type;

    template <typename T, T N, T S>
    struct st_make_sequence {
        using type =
            add_sequence_t<integer_sequence<T, S>, typename st_make_sequence<T, N, S + 1>::type>;
    };

    template <typename T, T N>
    struct st_make_sequence<T, N, N> {
        using type = integer_sequence<T, N>;
    };

    template <typename T, T N>
    using make_integer_sequence_t = st_make_sequence<T, N - 1, T{0}>::type;

    }  // namespace hidden

    template <typename SeqType1, typename SeqType2>
    using add_sequence_t = hidden::add_sequence_t<SeqType1, SeqType2>;

    template <typename T, T N>
    using make_integer_sequence_t = hidden::make_integer_sequence_t<T, N>;

    namespace hidden {
    template <typename Type1, typename Type2, typename... Types>
    struct st_is_same {
        static constexpr bool value =
            st_is_same<Type1, Type2>::value && st_is_same<Type1, Types...>::value;
    };

    template <typename Type1, typename Type2>
    struct st_is_same<Type1, Type2> {
        static constexpr bool value = false;
    };

    template <typename T>
    struct st_is_same<T, T> {
        static constexpr bool value = true;
    };

    template <typename Type1, typename Type2, typename... Types>
    constexpr auto is_same_v = st_is_same<Type1, Type2, Types...>::value;
    }  // namespace hidden

    template <typename Type1, typename Type2, typename... Types>
    constexpr auto is_same_v = hidden::is_same_v<Type1, Type2, Types...>;

    template <typename Type>
    concept integral_or_floating_point_c = std::is_integral_v<remove_const_reference_t<Type>> ||
                                           std::is_floating_point_v<remove_const_reference_t<Type>>;

    template <typename Type>
    concept floating_c = std::is_floating_point_v<remove_const_reference_t<Type>> ;
    namespace hidden {

    // ============ element many

    // primary class template
    template <typename T>
    struct st_container_element {
        using type = T;
        using element_type = T;
        static constexpr bool is_character_array = false;
        static constexpr size_t count = 0;
    };

    template <size_t N>
    struct st_container_element<char[N]>  // non refference
    {
        using type = char;
        using element_type = std::string;
        static constexpr bool is_character_array = true;
        static constexpr size_t count = N;
    };

    template <size_t N>
    struct st_container_element<wchar_t[N]>  // non refference
    {
        using type = wchar_t;  // wchar works with wcout, wcin ... , if not use wont work
        using element_type = std::string;
        static constexpr bool is_character_array = true;
        static constexpr size_t count = N;
    };
    template <>
    struct st_container_element<
        const wchar_t*>  // non refference wchar_t apparently doesnt work in modern cpp
    {
        using type = wchar_t;
        using element_type = std::wstring;
        static constexpr bool is_character_array = true;
        static constexpr size_t array_count = 0;
    };

    }  // namespace hidden
    template <typename T>
    constexpr auto is_character_array_v =
        hidden::st_container_element<jf::types::remove_const_reference_t<T>>::is_character_array;

    template <typename T>
    constexpr auto array_count_v = hidden::st_container_element<jf::types::remove_const_reference_t<T>>::count;

    namespace hidden {

    template <typename T, size_t N>
    struct st_container_element<T[N]> {
        using type = T;
        using element_type = std::vector<T>;
        static constexpr bool is_character_array = false;
    };

    template <typename T, size_t N>
    struct st_container_element<T (&)[N]> {
        using type = T;
        using element_type = std::vector<T>;
        static constexpr bool is_character_array = false;
    };

    template <typename T>
    struct st_container_element<std::initializer_list<T>> {
        using type = T;
        using element_type = std::vector<T>;
        static constexpr bool is_character_array = false;
    };

    }  // namespace hidden

    template <typename T>
    using container_element_flat_t =
        hidden::st_container_element<jf::types::remove_const_reference_t<T>>::type;

    template <typename T>
    using container_element_t =
        hidden::st_container_element<jf::types::remove_const_reference_t<T>>::element_type;

    // decay array to pointer
    template <typename T, size_t N>
    auto decay_array(T (&array)[N]) {
        return array;
    }

    template <typename Type, typename... Types>
    auto make_vctr(Type&& first, Types&&... args) {
        if constexpr (std::is_array_v<jf::types::remove_const_reference_t<Type>>) {
            static_assert(is_character_array_v<Type>, "should be character array");

            using element_t = jf::types::container_element_t<Type>;
            using container_t = std::vector<element_t>;

            if (auto last_character = first[jf::types::array_count_v<Type> - 1];
                jf::types::container_element_flat_t<Type>{} == last_character)
                return container_t{jf::types::decay_array(first), jf::types::decay_array(args)...};
            else {
                return container_t{element_t{std::cbegin(first), std::cend(first)},
                                   element_t{std::cbegin(args), std::cend(args)}...};
            }

        } else {
            using element_t = jf::types::remove_const_reference_t<Type>;
            return std::vector<element_t>{std::forward<Type>(first), std::forward<Types>(args)...};
        }
    }

    template <template <typename, typename...> class Cntr, typename Type, typename... Types>
    auto create_container(Type&& first, Types&&... args) {
        if constexpr (std::is_array_v<jf::types::remove_const_reference_t<Type>>) {
            static_assert(is_character_array_v<Type>, "should be character array");

            using element_t = container_element_t<Type>;
            using container_t = Cntr<element_t>;

            if (auto last_character = first[array_count_v<Type> - 1];
                container_element_flat_t<Type>{} == last_character)
                return container_t{decay_array(first), decay_array(args)...};
            else {
                return container_t{element_t{std::cbegin(first), std::cend(first)},
                                   element_t{std::cbegin(args), std::cend(args)}...};
            }

        } else {
            using element_t = jf::types::remove_const_reference_t<Type>;
            return Cntr<element_t>{std::forward<Type>(first), std::forward<Types>(args)...};
        }
    }


    template<typename T>
        concept function_c = std::is_function_v<remove_cv_ref_t<T>>;
    ////////////////////////////////////
    // convert const char* to std::string
    // convert const w_char* to std::string
    // convert array int[] to std::vector<int>
    // convert std::intialize_list to std::vector<>
    template <typename Type>
    auto element_to_container(Type&& arg) {
        //array ..
        if constexpr (std::is_array_v<jf::types::remove_const_reference_t<Type>>) {
            if constexpr (jf::types::is_character_array_v<Type>) {
                using element_t = jf::types::container_element_t<Type>;

                if (const auto& last_character = arg[array_count_v<Type> - 1];
                    jf::types::container_element_flat_t<Type>{} == last_character)
                    return element_t{jf::types::decay_array(arg)};
                else
                    return element_t{std::cbegin(arg), std::cend(arg)};

            } else  // non character array
            {
                using element_t = jf::types::container_element_t<Type>;
                return element_t{std::cbegin(arg), std::cend(arg)};
            }
        } else  //non array
        {
            using element_t = jf::types::container_element_t<Type>;
            return element_t{std::forward<Type>(arg)};
        }
    }

    // belongs to non deduced context
    template <typename Type>
    auto element_to_container(std::initializer_list<Type>& lst) {
        // both are the same
        // using element_t = container_element_t<decltype(lst)>;
        using element_t = jf::types::container_element_t<std::initializer_list<Type>>;
        return element_t{lst};
    }
    template <typename Type>
    auto element_to_container(std::initializer_list<Type>&& lst) {
        // both are the same
        // using element_t = container_element_t<decltype(lst)>;
        using element_t = jf::types::container_element_t<std::initializer_list<Type>>;
        return element_t{std::forward<std::initializer_list<Type>>(lst)};
    }
    template <typename Type, typename... Types>
    auto make_tuple(Type&& first, Types&&... args) {
        using container_t =
            std::tuple<jf::types::container_element_t<Type>, jf::types::container_element_t<Types>...>;
        return container_t{jf::types::element_to_container(std::forward<Type>(first)),
                           jf::types::element_to_container(std::forward<Types>(args))...};
    }


    namespace hidden {
    template <typename T, typename... Types>
    struct st_are_all_same;

    template <typename T1, typename T2>
    struct st_are_all_same<T1, T2> {
        static constexpr bool value = false;
    };

    template <typename T>
    struct st_are_all_same<T, T> {
        static constexpr bool value = true;
    };

    template <typename T>
    struct st_are_all_same<T> {
        static constexpr bool value = true;
    };

    template <typename T1, typename T2, typename... Types>
    struct st_are_all_same<T1, T2, Types...> {
        static constexpr bool value =
            st_are_all_same<T1, T2>::value && st_are_all_same<T2, Types...>::value;
    };
    }  // namespace hidden

    template <typename T, typename... Types>
    constexpr auto are_all_same_v = hidden::st_are_all_same<T, Types...>::value;

    template <typename T, typename... Types>
    constexpr auto are_all_same_flat_v =
        are_all_same_v<jf::types::remove_const_reference_t<T>, jf::types::remove_const_reference_t<Types>...>;

    namespace hidden {
    // such as std::conditional<>
    template <bool bWhich, typename TrueType, typename FalseType>
    struct st_conditional;

    template <typename TrueType, typename FalseType>
    struct st_conditional<true, TrueType, FalseType> {
        using type = TrueType;
    };

    template <typename TrueType, typename FalseType>
    struct st_conditional<false, TrueType, FalseType> {
        using type = FalseType;
    };

    }  // namespace hidden
    template <bool bWhich, typename TrueType, typename FalseType>
    using conditional_t = hidden::st_conditional<bWhich, TrueType, FalseType>::type;
    
    
    namespace hidden{
        template<typename T>
            struct st_is_tuple_or_array{
                static constexpr auto value = false;
            };
        template<typename ...Types>
            struct st_is_tuple_or_array<std::tuple<Types...>>
            {
                static constexpr auto value = true;
            };
        template<typename T, size_t Nm>
            struct st_is_tuple_or_array<std::array<T, Nm>>
            {
                static constexpr auto value = true;
            };

        template<typename T>
        struct st_is_pair : std::false_type {  };
        template<typename T1, typename T2>
        struct st_is_pair<std::pair<T1, T2>> : std::true_type {  };
    } // namespace hidden
    template<typename T>
        constexpr auto is_tuple_or_array_v = hidden::st_is_tuple_or_array<std::remove_cvref_t<T> >::value;
    template<typename T>
        concept tuple_or_array_c = is_tuple_or_array_v<T>;

    template<typename T>
        constexpr auto is_pair_v = hidden::st_is_pair<std::remove_cvref_t<T> >::value;
    template<typename T>
        concept pair_c = is_pair_v<T>;

    namespace hidden{
        template<typename T, T startVal, T endVal, T stepVal, T... Indices>
        constexpr auto make_sequences(std::integer_sequence<T, Indices...> seq)
        {
            if constexpr ( stepVal > 0 && startVal < endVal )  
            {
                return make_sequences<T, startVal + stepVal, endVal, stepVal>
                    ( std::integer_sequence<T, Indices..., startVal>{ }); 
            }
            else if constexpr( stepVal < 0 && startVal > endVal ) 
            {
                return make_sequences<T, startVal + stepVal, endVal, stepVal>
                    ( std::integer_sequence<T, Indices..., startVal>{}); 
            }
            else {
                return seq;
            }
        }
    } // namespace hidden
      
    template<typename T, T startVal, T endVal, T stepVal>
    using make_sequence = decltype(hidden::make_sequences<T, startVal, endVal, stepVal>( std::integer_sequence<T>{} ) );

    template<std::size_t... Indices>
    using sequence = std::index_sequence<Indices...>;

    template< typename TypeContainer, auto  endVal = std::tuple_size_v<std::remove_cvref_t<TypeContainer> >>
        requires (tuple_or_array_c<TypeContainer> || pair_c<TypeContainer>)
    auto lambda_seq( auto&& func) -> decltype(auto)
        requires requires {func(make_sequence<decltype(endVal), 0, endVal, 1>{});}
    {
        return func(make_sequence<decltype(endVal), 0, endVal, 1>{});
    }

    template< auto  endVal>
    auto lambda_seq( auto&& func) -> decltype(auto)
        requires requires {func(make_sequence<decltype(endVal), 0, endVal, 1>{});}
    {
        return func(make_sequence<decltype(endVal), 0, endVal, 1>{});
    }

    template< auto startVal, auto  endVal, auto stepVal>
    auto lambda_seq( auto&& func) -> decltype(auto)
        requires requires {func(make_sequence<decltype(endVal), startVal, endVal, stepVal>{});}
    {
        return func(make_sequence<decltype(endVal), startVal, endVal, stepVal>{});
    }

    template< auto startVal, auto  endVal>
    auto lambda_seq( auto&& func) -> decltype(auto)
        requires requires {func(make_sequence<decltype(endVal), startVal, endVal, (startVal < endVal)? 1 : -1>{});}
    {
        return func(make_sequence<decltype(endVal), startVal, endVal, (startVal < endVal)? 1 : -1>{});
    }
    template<auto ... Indices>
    auto lambda_seq( auto&& func, sequence<Indices...> seq) -> decltype(auto)
        requires requires{ func(seq); }
    {
        return func(seq);
    }

    template<auto ... Indices>
    auto lambda_seq( auto&& func, sequence<Indices...>) -> decltype(auto)
    {
        return [&func, args = std::forward_as_tuple(Indices...)] {return std::apply(func, args);}();
    }

}// namespace jf::types
