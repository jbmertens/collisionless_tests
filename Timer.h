#ifndef TIMER
#define TIMER

#include <ctime>
#include <iostream>
#include <map>
#include <iomanip>

/**
 * @brief Individual timer classes used by the TimerManager class.
 */
class Timer
{
private:
  double m_secs;
  struct timespec m_starttime;
  struct timespec m_stoptime;

public:
  Timer() : m_secs(0) {}

  inline double time() const { return m_secs; }

  inline void start() { clock_gettime(CLOCK_MONOTONIC, &m_starttime); }

  Timer operator+(const Timer &t2)
  {
    Timer t;
    t.m_secs = time() + t2.time();
    return t;
  }

  Timer operator-(const Timer &t2)
  {
    Timer t;
    t.m_secs = time() - t2.time();
    return t;
  }

  /**
   * @brief Stop timer
   */
  void stop()
  {
    clock_gettime(CLOCK_MONOTONIC, &m_stoptime);
    m_secs += (double)(m_stoptime.tv_sec - m_starttime.tv_sec);
    m_secs += (m_stoptime.tv_nsec - m_starttime.tv_nsec)*1e-9;
  }

  /**
   * @brief Reset timer
   */
  void reset()
  {
    m_secs = 0.;
    m_starttime.tv_sec  = 0;
    m_starttime.tv_nsec = 0;
    m_stoptime.tv_sec   = 0;
    m_stoptime.tv_nsec  = 0;
  }
};

/**
 * @brief TimerManager class containing multiple timers;
 * access individual timers via, eg, TM["my_timer"].start()
 */
class TimerManager
{
public:
  TimerManager() {};

  std::string getStateString()
  {
    std::map<std::string,Timer>::iterator it;

    std::string str = "==== TimerManager ====\n";
    
    for(it = m_timers.begin(); it != m_timers.end(); ++it)
    {
      std::string name = it->first;
      std::string seconds = std::to_string(it->second.time());
      str += "  " + name + ": " + seconds + "s\n";
    }

    return str;
  }

  inline Timer& operator[](std::string key) { return m_timers[key]; }

  friend std::ostream& operator<<(std::ostream &ostr, TimerManager T);

private:
  std::map<std::string, Timer> m_timers;
};

std::ostream& operator<<(std::ostream &ostr, const Timer &t)
{
  ostr << std::fixed << std::setprecision(3) << t.time() << "s";
  return ostr;
}

std::ostream& operator<<(std::ostream &ostr, TimerManager T)
{
  std::map<std::string,Timer>::iterator it;

  ostr << "==== TimerManager ====" << std::endl;
  for(it = T.m_timers.begin(); it != T.m_timers.end(); ++it) {
    ostr << "  " << it->first << ": " << it->second << std::endl;
  }

  return ostr;
}

#endif
