#ifndef HP_TIMER_HH
#define HP_TIMER_HH

#include <chrono>

#define HP_TIMER_OPCOUNTER oc::OpCounter<double>

#define HP_TIMER_DURATION(name) __hp_timer_##name##_duration
#define HP_TIMER_STARTTIME(name) __hp_timer_##name##_start
#define HP_TIMER_OPCOUNTERS_START(name) __hp_timer_##name##_counters_start
#define HP_TIMER_OPCOUNTERS(name) __hp_timer_##name##_counters
#define HP_TIMER_ELAPSED(name) std::chrono::duration_cast<std::chrono::duration<double> >( HP_TIMER_DURATION(name) ).count()

#ifdef ENABLE_HP_TIMERS

#ifdef ENABLE_COUNTER

#define HP_DECLARE_TIMER(name)                               \
  std::chrono::high_resolution_clock::duration HP_TIMER_DURATION(name);	\
  std::chrono::high_resolution_clock::time_point HP_TIMER_STARTTIME(name); \
  HP_TIMER_OPCOUNTER::Counters HP_TIMER_OPCOUNTERS_START(name); \
  HP_TIMER_OPCOUNTER::Counters HP_TIMER_OPCOUNTERS(name);

#define HP_TIMER_START(name) \
  do { \
  HP_TIMER_OPCOUNTERS_START(name) = HP_TIMER_OPCOUNTER::counters; \
  HP_TIMER_STARTTIME(name) = std::chrono::high_resolution_clock::now(); \
  } while(false)

#define HP_TIMER_STOP(name) \
  do { \
  std::chrono::high_resolution_clock::time_point __hp_end_time = std::chrono::high_resolution_clock::now(); \
  HP_TIMER_OPCOUNTERS(name) += HP_TIMER_OPCOUNTER::counters - HP_TIMER_OPCOUNTERS_START(name); \
  HP_TIMER_DURATION(name) += __hp_end_time - HP_TIMER_STARTTIME(name); \
  } while(false)

#define HP_TIMER_RESET(name) \
  do { \
    HP_TIMER_DURATION(name) = std::chrono::high_resolution_clock::duration::zero(); \
  HP_TIMER_OPCOUNTERS(name).reset(); \
  } while (false)

#else

#define HP_DECLARE_TIMER(name)                               \
  std::chrono::high_resolution_clock::duration HP_TIMER_DURATION(name);	\
  std::chrono::high_resolution_clock::time_point HP_TIMER_STARTTIME(name);

#define HP_TIMER_START(name) HP_TIMER_STARTTIME(name) = std::chrono::high_resolution_clock::now();
#define HP_TIMER_STOP(name) \
  do { \
  std::chrono::high_resolution_clock::time_point __hp_end_time = std::chrono::high_resolution_clock::now(); \
  HP_TIMER_DURATION(name) += __hp_end_time - HP_TIMER_STARTTIME(name); \
  } while(false)

#define HP_TIMER_RESET(name) HP_TIMER_DURATION(name) = std::chrono::high_resolution_clock::duration::zero();

#endif // ENABLE_COUNTER

#else // ENABLE_HP_TIMERS

#define HP_DECLARE_TIMER(name)
#define HP_TIMER_START(name)
#define HP_TIMER_STOP(name)
#define HP_TIMER_RESET(name)

#endif // ENABLE_HP_TIMERS

#endif // HP_TIMER_HH
