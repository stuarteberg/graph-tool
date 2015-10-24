// graph-tool -- a general graph modification and manipulation thingy
//
// Copyright (C) 2006-2015 Tiago de Paula Peixoto <tiago@skewed.de>
//
// This program is free software; you can redistribute it and/or
// modify it under the terms of the GNU General Public License
// as published by the Free Software Foundation; either version 3
// of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#ifndef NESTED_FOR_LOOP_HH
#define NESTED_FOR_LOOP_HH

#include <boost/mpl/for_each.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/empty.hpp>
#include <boost/any.hpp>

namespace boost
{
namespace mpl
{
// The following is a implementation of a nested for_each loop, which runs a
// given Action functor for each combination of its arguments, given by the type
// ranges, as such:
//
//     struct foo
//     {
//         template<class T1, class T2, class T3>
//         void operator()(T1, T2, T3) const
//         {
//             ...
//         }
//     };
//
//     ...
//
//     typedef mpl::vector<int,float,long> r1;
//     typedef mpl::vector<string,double> r2;
//     typedef mpl::vector<size_t,char> r3;
//
//     any x = float(2);
//     any y = string("foo");
//     any z = size_t(42);
//
//     bool found = nested_for_each<r1,r2,r3>(foo(), x, y, z);
//
// The code above will run iterate through all combinations of foo::operator(T1,
// T2, T3) and call the one that corresponds to the actual types stored in x, y,
// and z. If the types are not found during iteration, we have found = true,
// otherwise found = false. This provides a more general compile-time to
// run-time bridge than the simpler mpl::for_each().

template <class Action, class T>
struct bind_arg
{
    bind_arg(Action a, any& arg, bool& found)
        : _a(a), _arg(arg), _found(found) {}

    template <class... Ts>
    __attribute__((always_inline))
    void operator()(Ts&&... args) const
    {
        T* v = const_cast<T*>(any_cast<T>(&_arg));
        if (v != 0)
            _a(*v, args...);
    }

    __attribute__((always_inline))
    void operator()() const
    {
        T* v = const_cast<T*>(any_cast<T>(&_arg));
        if (v != 0)
        {
            _a(*v);
            _found = true;
        }
    }

    Action _a;
    any& _arg;
    bool& _found;
};


template <class Action>
struct dispatch
{
    dispatch(Action a, any* args, bool& found)
        : _a(a), _args(args), _found(found) {}

    template <class T>
    __attribute__((always_inline))
    auto get_next() const
    {
        bind_arg<Action, T> a(_a, *_args, _found);
        return dispatch<bind_arg<Action, T>>(a, _args + 1, _found);
    }

    template <class T>
    __attribute__((always_inline))
    void operator()(T) const
    {
        bind_arg<Action, T> a(_a, *_args, _found);
        a();
    }

    __attribute__((always_inline))
    void operator()() const
    {
        _a();
    }

    Action _a;
    any* _args;
    bool& _found;
};

template <class F>
inline void for_each_pack(F)
{};

template <class F, class T, class... Ts>
inline void for_each_pack(F f)
{
    f(T());
    for_each_pack<F, Ts...>(f);
};

template <class F, class Seq, class Iter, class... Ts>
inline void for_each_alt(F f, Iter)
{
    for_each_alt<F, Seq, typename next<Iter>::type, Ts...,
                 typename deref<Iter>::type>
        (f, typename next<Iter>::type());
}

template <class F, class Seq, class Iter, class... Ts>
inline void for_each_alt(F f, typename end<Seq>::type)
{
    for_each_pack<F, Ts...>(f);
}

template <class TR1, class... TRS, class Action>
void nested_for_each_imp(Action a);

template <class Action, class... TRS>
struct inner_loop
{
    inner_loop(Action a): _a(a) {}

    template <class T>
    __attribute__((always_inline))
    void operator()(T) const
    {
        nested_for_each_imp<TRS...>(_a.template get_next<T>());
    }
    Action _a;
};

template <class TR1, class... TRS, class Action>
void nested_for_each_imp(Action a)
{
    for_each_alt<inner_loop<Action, TRS...>,
                 TR1, typename begin<TR1>::type>
        (inner_loop<Action, TRS...>(a), typename begin<TR1>::type());
}

template <class Action>
void nested_for_each_imp(Action a)
{
    a();
}

template <class... TRS, class Action, class... Args>
bool nested_for_each(Action a, Args... args)
{
    bool found = false;
    std::array<any, sizeof...(args)> as{args...};
    auto b = dispatch<Action>(a, &as[0], found);
    nested_for_each_imp<TRS...>(b);
    return found;
}

} // mpl namespace
} // boost namespace

#endif //NESTED_FOR_LOOP_HH
