#ifndef MESSAGE_H
#define MESSAGE_H

#include <iostream>
#include <cstdlib>
#include <unistd.h>
#include <chrono>
#include <cmath>

namespace Stream
{
class Message
{
public:
    Message(const std::string& name);

    std::ostream& print( const uint32_t level=0 );

    void printerr( const std::string& message );

    void loadbar(const double fraction);

    void timerstart();

    void timerstop();

private:
    std::string m_name; ///< The name of the class calling
    uint32_t m_cols;

    std::chrono::time_point<std::chrono::high_resolution_clock> m_time;
    bool m_clockstarted;

    char* m_termtype;
    char m_term_buffer[2048];
};
}
#endif // MESSAGE_H

