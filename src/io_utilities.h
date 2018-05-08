#pragma once

#include "config.h"

#include <boost/lexical_cast.hpp>
#include <omp.h>
#include <string>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdio>

/************************************************************************/
/* print program header: description, version and contacts              */
/************************************************************************/
void print_header();

/************************************************************************/
/* timing functions                                                     */
/************************************************************************/
clock_t elapsed_seconds();

/************************************************************************/
/* messaging functions                                                  */
/************************************************************************/
template<class MsgT>
void print_tm(const MsgT& msg, std::ostream& os = std::cout)
{
    size_t time_steps_width = boost::lexical_cast<std::string>(PIC::Config::time_steps()).size();
    os << std::setw(time_steps_width)
       << elapsed_seconds() << " sec: "
       << msg
       << std::endl;
}

template<class MsgPrefixT, class MsgSuffixT>
void print_tm(const MsgPrefixT& msg_prefix, const MsgSuffixT& msg_suffix, std::ostream& os = std::cout)
{
    size_t time_steps_width = boost::lexical_cast<std::string>(PIC::Config::time_steps()).size();

    os << std::setw(time_steps_width)
       << elapsed_seconds() << " sec: "
       << msg_prefix << msg_suffix
       << std::endl;
}

template<class MsgT>
void print(const MsgT& msg, std::ostream& os = std::cout)
{
    os << "\n" << msg << std::endl;
}

/************************************************************************/
/* Create data file name for saving                                     */
/************************************************************************/
std::string create_out_file_name(const std::string& prefix,
                                 const std::string& type_name,
                                 index_t step);

/************************************************************************/
/* Logging function                                                     */
/************************************************************************/
template<class MsgT>
void log(const MsgT& msg, bool duplicate_on_stdout = true)
{
    const int thread_num = omp_get_thread_num();

    char log_file_name[50];
    sprintf(log_file_name, "opic_thread_%d.log", thread_num);
    std::ofstream ofs_log(log_file_name, std::ios_base::app);

    if (!ofs_log)
    {
        char error_msg[200];
        sprintf(error_msg, "Can't open \"%s\" file.", log_file_name);
        std::cout << error_msg;
        return;
    }

    ofs_log << setiosflags(std::ios_base::scientific);

    print(msg, ofs_log);

    if(duplicate_on_stdout)
    {
        print(msg);
    }
}
