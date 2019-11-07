#pragma once

#include "BufferStructure.h"

namespace {

   template<class T>
   void finish_futures(T& first_futures) {
      for (auto& future : first_futures)
         future.get();
   }

   template <class T, class... Args>
   void finish_futures(T& first_futures, Args&... args) {
      for (auto& future : first_futures)
         future.get();
      finish_futures(args...);
   }

} // namespace {}


template<class T>
magneto::BufferStructure<T>::BufferStructure(
   const std::function<T()>& generator_fun, const int max_rng_threads
)
   : m_buffer(generator_fun()) // Fill the buffer
   , m_buffer_filler(generator_fun)
{
   // Start computing thread(s) for next buffer content immediately
   m_futures.reserve(max_rng_threads);
   for (int i = 0; i < max_rng_threads; ++i) {
      m_futures.emplace_back(std::async(std::launch::async, m_buffer_filler));
   }
}


template<class T>
magneto::BufferStructure<T>::~BufferStructure() {
   // Wait for any running futures to finish
   finish_futures(m_futures);
}


template<class T>
const T& magneto::BufferStructure<T>::get_buffer() const {
   return m_buffer;
}


template<class T>
void magneto::BufferStructure<T>::refill() {
   // Find a future that is already finished, or wait for one to finish.
   // future.is_ready() is ironically not yet in the standard, but it is in MSVC... and pretty neat
   while (true) {
      for (auto& future : m_futures) {
         if (!future._Is_ready())
            continue;

         // fill the buffer with the result of the future
         m_buffer = std::move(future.get());

         // restart the future to start computing the next one
         future = std::async(std::launch::async, m_buffer_filler);
         return;
      }
   }
}
