#ifndef TIMER_H_
#define TIMER_H_

#include <chrono> // std::chrono
#include <string> // std::string
#include <sstream> // std::stringstream
#include <iomanip> // std::put_time

struct LocalTimer
{
    using clock = std::chrono::high_resolution_clock;
    using duration = std::chrono::duration<double>;
    using time_point = std::chrono::time_point<clock, duration>;

    time_point start, end;

    void start_timer() { start = clock::now(); }
    void stop_timer()  { end = clock::now(); }
    double get_elapsed() { return duration(end - start).count(); }
};

std::string return_current_date_and_time()
{
    /*
     * Shamelessly copied from here: https://stackoverflow.com/questions/17223096/outputting-date-and-time-in-c-using-stdchrono
     */

    auto now = std::chrono::system_clock::now();
    auto in_time_t = std::chrono::system_clock::to_time_t(now);

    std::stringstream ss;
    ss << std::put_time(std::localtime(&in_time_t), "%Y/%m/%d %X");
    return ss.str();
}

#endif
