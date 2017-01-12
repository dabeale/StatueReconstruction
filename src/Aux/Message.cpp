
#include "Message.h"

namespace Stream
{
Message::Message(const std::string& name) : m_name(name),
    m_clockstarted(false)
{
    m_termtype = std::getenv("TERM");
    //int success = tgetent (m_term_buffer, m_termtype);
    // tgetnum ("li") << "," <<  tgetnum ("co")
/*
    if(m_termtype == NULL || success==0)
    {
        m_cols = 100;
    }
    else
    {
        m_cols = tgetnum ((char *)"co");
    }
    */
    m_cols = 30;
}

std::ostream& Message::print( const uint32_t level )
{
    std::cout << m_name << ": ";
    for(uint32_t k=0; k<level; ++k)
    {
        std::cout << "-";
    }
    if(level>0)
    {
        std::cout << ">";
    }
    std::cout <<" ";
    std::fflush(stdout);
    return std::cout;
}

void Message::printerr( const std::string& message )
{
    std::cerr << m_name << ": **" << message << "**" << std::endl;
}

void Message::loadbar(const double fraction)
{
    std::cout << static_cast<int>(fraction*100) << "% ";
    std::cout << "[";
    for(uint32_t k=0; k<m_cols; ++k )
    {
        if( k < fraction*m_cols)
        {
            std::cout << "-";
        }
        else
        {
            std::cout << " ";
        }
    }
    std::cout << "]\r";
    if(std::abs(fraction - 1.0) < 1e-6)
    {
        std::cout << "\n";
    }
    std::fflush(stdout);
}

void Message::timerstart()
{
    if(!m_clockstarted)
    {
        m_clockstarted=true;
        m_time=std::chrono::high_resolution_clock::now();
    }
    else
    {
        printerr("The clock is already started");
    }
}

void Message::timerstop()
{
    auto tend = std::chrono::high_resolution_clock::now();
    print(1) << std::chrono::duration<double, std::milli>(tend-m_time).count()/1000 << " seconds" << std::endl;
    m_clockstarted=false;
}

}
