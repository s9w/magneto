#pragma once

#include <future>


namespace magneto {

   /// <summary>Structure that bundles a data buffer, the futures that compute the next buffer
   /// content concurrently, and the function to compute them all in one.</summary>
   template<class T>
   class BufferStructure {
   public:
      BufferStructure(const std::function<T()>& generator_fun, const int max_rng_threads);
      ~BufferStructure();

      /// <summary>Returns reference to the buffer content</summary>
      [[nodiscard]] const T& get_buffer() const;

      /// <summary>Refills the internal buffer with result from finished computing threat and restarts that it.</summary>
      void refill();

   private:
      void move_future_results_and_restart(std::future<T>& future);

      T m_buffer;
      std::vector<std::future<T>> m_futures;
      std::function<T()> m_buffer_filler;
   };
}

#include "BufferStructure.hpp"
