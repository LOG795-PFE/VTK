/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkImplicitArrayTraits.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef vtkImplicitArrayTraits_h
#define vtkImplicitArrayTraits_h

#include "vtkSystemIncludes.h"
#include <numeric>
#include <type_traits>

/**
 * This file contains the traits for the implicit array mechanism in VTK. These traits are very much
 * an internal to vtkImplicitArrays and normal developers looking to develop a new vtkImplicitArray
 * should (ideally) not have to open this file.
 *
 * In order to ensure that template parameters passed to the vtkImplicitArray share a common
 * interface without having to subclass all of them from the same abstract class, we have decided to
 * use a trait mechanism to statically dispatch the functionalities of types passed as template
 * parameters to the array.
 *
 * There is 1 mandatory traits that a template type to vtkImplicitArray must implement:
 * - has_map_trait || is_closure_trait: ensures an implementation of int -> value
 *
 * Potential improvements to implicit arrays which would allow for write access would include the
 * following 2 optional traits:
 * - has_insert_trait || is_reference_closure_trait: provides an implementation to update the
 * internals of the template type to add or set new values to the array
 * - has_remove_trait: provides an implementation to update the internals of the template to remove
 * values from the array
 *
 * All the traits defining the behavior of the implicit "function" or "backend" to the
 * vtkImplicitArray should be composited into the implicit_array_traits
 */

namespace vtk
{
namespace detail
{

VTK_ABI_NAMESPACE_BEGIN

template <typename... T>
using void_t = void;

///@{
/**
 * \struct has_map_trait
 * \brief used to check whether the template type has a method named map
 */
template <typename, typename = void>
struct has_map_trait : std::false_type
{
};

template <typename T>
struct has_map_trait<T, void_t<decltype(&std::remove_reference<T>::type::map)>>
  : public has_map_trait<decltype(&std::remove_reference<T>::type::map)>
{
  using type = T;
};

template <typename T>
struct has_map_trait<T*> : public has_map_trait<T>
{
};

template <typename T>
struct has_map_trait<const T> : public has_map_trait<T>
{
};

template <typename R, typename T, typename Arg>
struct has_map_trait<R (T::*)(Arg) const> : public has_map_trait<R(Arg)>
{
};

template <typename R, typename Arg>
struct has_map_trait<R(Arg)>
{
  static_assert(std::is_integral<Arg>::value, "Argument to map must be integral type");
  static constexpr bool value = true;
  using rtype = R;
};
///@}

///@{
/**
 * \struct is_closure_trait
 * \brief A trait determining whether an object acts like a mono-variable integer closure
 */
template <typename, typename = void>
struct is_closure_trait : std::false_type
{
};

template <typename Closure>
struct is_closure_trait<Closure, void_t<decltype(&Closure::operator())>>
  : public is_closure_trait<decltype(&Closure::operator())>
{
  using type = Closure;
};

template <typename Closure>
struct is_closure_trait<Closure*> : public is_closure_trait<Closure>
{
};

template <typename Closure>
struct is_closure_trait<const Closure> : public is_closure_trait<Closure>
{
};

template <typename Closure>
struct is_closure_trait<Closure&> : public is_closure_trait<Closure>
{
};

template <typename Closure, typename R, typename Arg>
struct is_closure_trait<R (Closure::*)(Arg) const> : public is_closure_trait<R(Arg)>
{
};

template <typename R, typename Arg>
struct is_closure_trait<R (*)(Arg)> : public is_closure_trait<R(Arg)>
{
};

template <typename R, typename Arg>
struct is_closure_trait<R(Arg)>
{
  static_assert(std::is_integral<Arg>::value, "Argument to closure must be integral type");
  static constexpr bool value = true;
  using rtype = R;
};
///@}

namespace iarrays
{
/**
 * \enum ReadOperatorCodes
 * \brief An enum for formalizing the different trait types accepted for defining a "readable"
 * object
 */
enum ReadOperatorCodes
{
  NONE,
  MAP,
  CLOSURE
};
}

///@{
/**
 * \struct can_map_trait
 * \brief An intermediate trait for exposing a unified trait interface
 */
template <typename T, typename = void>
struct can_map_trait
{
  using type = T;
  static constexpr bool value = false;
  using rtype = void;
  static constexpr iarrays::ReadOperatorCodes code = iarrays::NONE;
};

template <typename T>
struct can_map_trait<T, void_t<typename has_map_trait<T>::rtype>>
{
  using type = T;
  static constexpr bool value = true;
  using rtype = typename has_map_trait<T>::rtype;
  static constexpr iarrays::ReadOperatorCodes code = iarrays::MAP;
};
///@}

///@{
/**
 * \struct can_close_trait
 * \brief An intermediate trait for exposing a unified trait interface
 */
template <typename T, typename = void>
struct can_close_trait
{
  using type = T;
  static constexpr bool value = false;
  using rtype = void;
  static constexpr iarrays::ReadOperatorCodes code = iarrays::NONE;
};

template <typename T>
struct can_close_trait<T, void_t<typename is_closure_trait<T>::rtype>>
{
  using type = T;
  static constexpr bool value = true;
  using rtype = typename is_closure_trait<T>::rtype;
  static constexpr iarrays::ReadOperatorCodes code = iarrays::CLOSURE;
};
///@}

/**
 * \struct implicit_array_traits
 * \brief A composite trait for handling all the different capabilities a "backend" to an
 * implicit array can have
 */
template <typename T>
struct implicit_array_traits
{
  using type = T;
  using trait =
    typename std::conditional<can_map_trait<T>::value, can_map_trait<T>, can_close_trait<T>>::type;
  static constexpr bool can_read = trait::value;
  using rtype = typename trait::rtype;
  static constexpr iarrays::ReadOperatorCodes code = trait::code;
  static constexpr bool default_constructible = std::is_default_constructible<T>::value;
};

VTK_ABI_NAMESPACE_END

} // detail
} // vtk

#endif // vtkImplicitArrayTraits_h
