/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkSMPToolsImpl.txx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/

#ifndef STDThreadvtkSMPToolsImpl_txx
#define STDThreadvtkSMPToolsImpl_txx

#include <algorithm>  // For std::sort
#include <functional> // For std::bind

#include "SMP/Common/vtkSMPToolsImpl.h"
#include "SMP/Common/vtkSMPToolsInternal.h" // For common vtk smp class
#include "SMP/STDThread/vtkSMPThreadPool.h" // For vtkSMPThreadPool
#include "vtkCommonCoreModule.h"            // For export macro

namespace vtk
{
namespace detail
{
namespace smp
{
VTK_ABI_NAMESPACE_BEGIN

int VTKCOMMONCORE_EXPORT GetNumberOfThreadsSTDThread();
bool VTKCOMMONCORE_EXPORT GetSingleThreadSTDThread();
void VTKCOMMONCORE_EXPORT PushThreadId(std::thread::id id);
void VTKCOMMONCORE_EXPORT PopThreadId();

//--------------------------------------------------------------------------------
template <typename FunctorInternal>
void ExecuteFunctorSTDThread(void* functor, vtkIdType from, vtkIdType grain, vtkIdType last)
{
  const vtkIdType to = std::min(from + grain, last);

  FunctorInternal& fi = *reinterpret_cast<FunctorInternal*>(functor);
  fi.Execute(from, to);
}

//--------------------------------------------------------------------------------
template <>
template <typename FunctorInternal>
void vtkSMPToolsImpl<BackendType::STDThread>::For(
  vtkIdType first, vtkIdType last, vtkIdType grain, FunctorInternal& fi)
{
  vtkIdType n = last - first;
  if (n <= 0)
  {
    return;
  }

  if (grain >= n || (!this->NestedActivated && this->IsParallel))
  {
    fi.Execute(first, last);
  }
  else
  {
    int threadNumber = GetNumberOfThreadsSTDThread();

    if (grain <= 0)
    {
      vtkIdType estimateGrain = (last - first) / (threadNumber * 4);
      grain = (estimateGrain > 0) ? estimateGrain : 1;
    }

    // /!\ This behaviour should be changed if we want more control on nested
    // (e.g only the 2 first nested For are in parallel)
    bool fromParallelCode = this->IsParallel.exchange(true);

    vtkSMPThreadPool pool(threadNumber);
    PushThreadId((pool.GetThreads())->at(0).get_id());

    for (vtkIdType from = first; from < last; from += grain)
    {
      auto job = std::bind(ExecuteFunctorSTDThread<FunctorInternal>, &fi, from, grain, last);
      pool.DoJob(job);
    }
    pool.Join();

    PopThreadId();
    // Atomic contortion to achieve this->IsParallel &= fromParallelCode.
    // This compare&exchange basically boils down to:
    // if (IsParallel == trueFlag)
    //   IsParallel = fromParallelCode;
    // else
    //   trueFlag = IsParallel;
    // Which either leaves IsParallel as false or sets it to fromParallelCode (i.e. &=).
    // Note that the return value of compare_exchange_weak() is not needed,
    // and that no looping is necessary.
    bool trueFlag = true;
    this->IsParallel.compare_exchange_weak(trueFlag, fromParallelCode);
  }
}

//--------------------------------------------------------------------------------
template <>
template <typename InputIt, typename OutputIt, typename Functor>
void vtkSMPToolsImpl<BackendType::STDThread>::Transform(
  InputIt inBegin, InputIt inEnd, OutputIt outBegin, Functor transform)
{
  auto size = std::distance(inBegin, inEnd);

  UnaryTransformCall<InputIt, OutputIt, Functor> exec(inBegin, outBegin, transform);
  this->For(0, size, 0, exec);
}

//--------------------------------------------------------------------------------
template <>
template <typename InputIt1, typename InputIt2, typename OutputIt, typename Functor>
void vtkSMPToolsImpl<BackendType::STDThread>::Transform(
  InputIt1 inBegin1, InputIt1 inEnd, InputIt2 inBegin2, OutputIt outBegin, Functor transform)
{
  auto size = std::distance(inBegin1, inEnd);

  BinaryTransformCall<InputIt1, InputIt2, OutputIt, Functor> exec(
    inBegin1, inBegin2, outBegin, transform);
  this->For(0, size, 0, exec);
}

//--------------------------------------------------------------------------------
template <>
template <typename Iterator, typename T>
void vtkSMPToolsImpl<BackendType::STDThread>::Fill(Iterator begin, Iterator end, const T& value)
{
  auto size = std::distance(begin, end);

  FillFunctor<T> fill(value);
  UnaryTransformCall<Iterator, Iterator, FillFunctor<T>> exec(begin, begin, fill);
  this->For(0, size, 0, exec);
}

//--------------------------------------------------------------------------------
template <>
template <typename RandomAccessIterator>
void vtkSMPToolsImpl<BackendType::STDThread>::Sort(
  RandomAccessIterator begin, RandomAccessIterator end)
{
  std::sort(begin, end);
}

//--------------------------------------------------------------------------------
template <>
template <typename RandomAccessIterator, typename Compare>
void vtkSMPToolsImpl<BackendType::STDThread>::Sort(
  RandomAccessIterator begin, RandomAccessIterator end, Compare comp)
{
  std::sort(begin, end, comp);
}

//--------------------------------------------------------------------------------
template <>
void vtkSMPToolsImpl<BackendType::STDThread>::Initialize(int);

//--------------------------------------------------------------------------------
template <>
int vtkSMPToolsImpl<BackendType::STDThread>::GetEstimatedNumberOfThreads();

//--------------------------------------------------------------------------------
template <>
bool vtkSMPToolsImpl<BackendType::STDThread>::GetSingleThread();

VTK_ABI_NAMESPACE_END
} // namespace smp
} // namespace detail
} // namespace vtk

#endif
