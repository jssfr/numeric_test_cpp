module;

export module math:Compose;

export namespace jf::math::compose
{
// binary function composition for arbitrary types
template<class Func, class Gfunc> auto compose(Func fun, Gfunc gfun)
{
  return [fun, gfun](auto param) { return f(g(param)); };
}

// for convienience
template<class Func, class Gfunc> auto operator*(Func fun, Gfunc gfun)
{
  return compose(fun, gfun);
}

// composition for n arguments
template<class Func, typename... Fns> auto compose(Func fun, Fns&... fns)
{
  return fun * compose(fns...);
}

// composition for n copies of f
template<int Var, class Func>
// must wrap chain in a struct to allow partial template specialization
struct multi
{
  static Func chain(Func fun) { return fun * multi<Var - 1, Func>::chain(fun); }
};
template<class Func> struct multi<2, Func>
{
  static Func chain(Func fun) { return fun * fun; }
};

/*
        EXAMPLE
            function<double(int)> f = [](auto i) { return i + 3; };
            function<int(double)> g = [](auto i) { return i * 2; };
            function<int(int)> h = [](auto i) { return i + 1; };

            std::cout << '\n'
            << "9   == " << compose(f, g, f)(0) << '\n'
            << "5   == " << (f * g * h)(0) << '\n'
            << "100 == " << compose<100>(h)(0) << '\n';
*/
template<int Var, class Func> Func compose(Func fun)
{
  return multi<Var, Func>::chain(fun);
}
}  // namespace jf::math::compose
