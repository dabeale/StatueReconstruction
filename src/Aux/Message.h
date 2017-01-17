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

/**
 * \brief A collection of functions and classes to write to stdout
 */
namespace Stream
{

/**
 * @brief The Message class
 * @details This class is for printing messages to the command line. It provides a few basic methods for printing titles,
 * error messages and a loadbar. The object retains a name which is printed out before any text in order to identify
 * the class which is producing the messages.
 */
class Message
{
public:
    /**
     * @brief Message
     * @details A class for printing informative text to the output.
     * @param name The name of the class which is producing the output
     */
    Message(const std::string& name);

    /**
     * @brief print Print text to std out.
     * @details The method returns a reference to the output stream, with text already streamed to it.
     * The messages are of the format:
     *  [ClassName] --> Text output
     * The size of the arrow is determined by the level.
     * @param level The level of debug message.
     * @return The output stream.
     */
    std::ostream& print( const uint32_t level=0 );

    /**
     * @brief printerr Print an error message to stderr
     * @details The message is highlighted with stars
     * @param message The message to be printed
     */
    void printerr( const std::string& message );

    /**
     * @brief loadbar Display a loadbar
     * @details The load bar is of the format
     * [---------------------]
     * It has a fixed width, and the size is determined by the input parameter 'fraction'.
     * @param fraction A number between 0 and 1. The method must be called with a 1 in order to print the newline character.
     */
    void loadbar(const double fraction);

    /**
     * @brief timerstart Start the timer.
     */
    void timerstart();

    /**
     * @brief timerstop Stop the timer
     * @details This method will print the amount of time elapsed since timerstart was called.
     */
    void timerstop();
private:
    std::string m_name; ///< The name of the class calling
    uint32_t m_cols; ///< The number of columns to use for the load bar.

    std::chrono::time_point<std::chrono::high_resolution_clock> m_time;
    bool m_clockstarted; ///< True if the clock was started

    char* m_termtype; ///< This is not used
    char m_term_buffer[2048]; ///< This is not used
};
}
#endif // MESSAGE_H

