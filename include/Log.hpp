#pragma once

#include <sstream>
#include <boost/format.hpp>
#include <iostream>

using namespace std;

enum log_level_t {
    LOG_NOTHING,
    LOG_ERROR,
    LOG_BASIC,
    LOG_VERBOSE,
    LOG_DEBUG
};
//Global variable.
log_level_t GLOBAL_LEVEL;

namespace log_impl {
    class formatted_log_t {
        public:
            formatted_log_t(log_level_t level, const wchar_t* msg) 
                :
                    fmt(msg),
                    level(level)
            {}
            ~formatted_log_t()
            {
                // GLOBAL_LEVEL is a global variable and could be changed at runtime
                // customization could be here.
                if (level <= GLOBAL_LEVEL) wcout << fmt << endl;
            }
            template <typename T>
            formatted_log_t& operator %(T value) {
                fmt % value;
                return *this;
            }
        protected:
            log_level_t level;
            boost::wformat fmt;
    };
}
//namespace log_impl
// Helper function. Class formatted_log_t will not be used directly
// This function needs to be called eg log<LOG_BASIC>(L"Basic result is %1 divided by %2.") % 1 % L"Two";
template <log_level_t level>
log_impl::formatted_log_t log(const wchar_t* msg) {
    return log_impl::formatted_log_t(level, msg);
}
