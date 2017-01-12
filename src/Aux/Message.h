#ifndef MESSAGE_H
#define MESSAGE_H

/* ----------------------------------------------------------------------
 * Copyright (C) 2016 Daniel Beale. All rights reserved.
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * ---------------------------------------------------------------------- */

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

