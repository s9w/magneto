#include "LatticeAlgorithms.h"
#include "system.h"

#include <sstream>

#include "spdlog/spdlog.h"

namespace {

   constexpr int get_exp_buffer_offset(int J) {
      // std::abs is not constexpr :(
      return J > 0 ? 8 * J : -8 * J;
   }

   std::string thread_id_to_string(const std::thread::id& id) {
      std::stringstream ss;
      ss << id;
      return ss.str();
   }


   magneto::UniformRandomBufferReturn get_random_buffer(const size_t buffer_size) {
      const std::string thread_id = thread_id_to_string(std::this_thread::get_id());
      unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
      std::mt19937_64 rng(seed);
      std::uniform_real_distribution <double > dist_one(0.0, 1.0);
      std::vector<double> normal_random_vector;
      normal_random_vector.reserve(buffer_size);
      for (int i = 0; i < buffer_size; ++i)
         normal_random_vector.emplace_back(dist_one(rng));
      spdlog::get("basic_logger")->debug("get_random_buffer() done from thread {}", thread_id);
      return { normal_random_vector, thread_id };
   }


   magneto::IndexPairVectorReturn get_lattice_indices(const size_t buffer_size, const int lattice_size) {
      const std::string thread_id = thread_id_to_string(std::this_thread::get_id());
      unsigned int seed = static_cast<unsigned int>(std::chrono::system_clock::now().time_since_epoch().count());
      std::mt19937_64 rng(seed);
      std::uniform_int_distribution<> dist_lattice(0, lattice_size - 1);
      magneto::IndexPairVector indices;
      indices.reserve(buffer_size);
      for (int i = 0; i < buffer_size; ++i)
         indices.emplace_back(dist_lattice(rng), dist_lattice(rng));
      spdlog::get("basic_logger")->debug("get_lattice_indices() done from thread {}", thread_id);
      return { indices, thread_id };
   }


   template <class T, class Tf, typename... Args>
   void assign_target_buffer_and_relaunch_future(
      T& target_buffer,
      std::future<magneto::ResultWithThreadId<T>>& future,
      const Tf& thread_function,
      Args... future_args
   ) {
      future.wait();
      const auto fresult = future.get();
      spdlog::get("basic_logger")->debug("fetched results and restarting. thread: {}", fresult.m_thread_id);
      target_buffer = std::move(fresult.m_result);
      future = std::async(std::launch::async, thread_function, future_args...);
   }


   /// <summary>Fills the target buffer with the result of any of the futures. If no future is finished,
   /// it's being waited until one is. The future is restarted at the end.</summary>
   template <class T, class Tf, typename... Args>
   void refill_target_buffer_from_futures(
      std::vector<std::future<magneto::ResultWithThreadId<T>>>& futures,
      T& target_buffer,
      const Tf& thread_function,
      Args... future_args
   ) {
      // usually there should be a thread already finished
      while (true) {
         for (auto& future : futures) {
            // future.is_ready() is ironically not yet in the standard, but it is in MSVC... and way too convenient to ignore
            if (future._Is_ready()) {
               assign_target_buffer_and_relaunch_future(target_buffer, future, thread_function, future_args...);
               return;
            }
         }
      }
   }

} // namespace {}


magneto::Metropolis::Metropolis(const int J, const double T, const int L, const int max_rng_threads)
   : m_cached_exp_values(get_cached_exp_values(J, T))
   , m_J(J)
{
   const size_t random_buffer_size = L * L;

   auto get_random_buffer_baked = [&] {return get_random_buffer(random_buffer_size); };
   auto get_lattice_indices_baked = [&] {return get_lattice_indices(random_buffer_size, L); };
   m_random_buffer = get_random_buffer_baked();
   m_lattice_index_buffer = get_lattice_indices_baked();

   m_future_random_buffers.reserve(max_rng_threads);
   m_lattice_index_futures.reserve(max_rng_threads);
   for (int i = 0; i < max_rng_threads; ++i) {
      m_future_random_buffers.emplace_back(std::async(std::launch::async, get_random_buffer, random_buffer_size));
      m_lattice_index_futures.emplace_back(std::async(std::launch::async, get_lattice_indices, random_buffer_size, L));
   }
}

magneto::Metropolis::~Metropolis(){
   // letting the future threads finish
   spdlog::get("basic_logger")->info("Destructor: {} lattice_index_futures and {} future_random_buffers", m_lattice_index_futures.size(), m_future_random_buffers.size());
   for (auto& future : m_lattice_index_futures)
   	future.get();
   for (auto& future : m_future_random_buffers)
      future.get();
}


void magneto::Metropolis::run(LatticeType& lattice){
   const int L = static_cast<int>(lattice.size());
   std::uniform_int_distribution <int> dist_grid(0, L - 1);
   std::uniform_real_distribution <double > dist_one(0.0, 1.0);
   const int buffer_offset = get_exp_buffer_offset(m_J);
   int flip_i, flip_j;
   int dE;
   for (int i = 0; i < m_random_buffer.size(); ++i) {
      flip_i = m_lattice_index_buffer[i].first;
      flip_j = m_lattice_index_buffer[i].second;
      dE = m_J * get_dE(lattice, flip_i, flip_j);
      if (dE <= 0 || (m_random_buffer[i] < m_cached_exp_values[dE + buffer_offset]))
         lattice[flip_i][flip_j] *= -1;
   }

   const size_t random_buffer_size = L * L;
   refill_target_buffer_from_futures(m_future_random_buffers, m_random_buffer, get_random_buffer, random_buffer_size);
   refill_target_buffer_from_futures(m_lattice_index_futures, m_lattice_index_buffer, get_lattice_indices, random_buffer_size, L);
}


std::vector<double> magneto::get_cached_exp_values(const int J, const double T) {
   std::vector<double> exp_values;
   const int min_value = -8 * J;
   const int max_value = 8 * J;
   const int value_count = max_value - min_value + 1;
   const int buffer_offset = get_exp_buffer_offset(J);
   for (int dE = 0; dE < value_count; ++dE)
      exp_values.push_back(exp(-(dE - buffer_offset) / T));
   return exp_values;
}
