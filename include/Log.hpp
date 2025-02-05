#ifndef _LOG_H
#define _LOG_H

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

// Global variable.
// This is like a forward declaration for functions
// The actual declaration is in Main.cpp
extern log_level_t GLOBAL_VERBOSITY_LEVEL;

namespace log_impl {
    class formatted_log_t {
        public:
            formatted_log_t(log_level_t level, const char* msg) 
                :
                    fmt(msg),
                    level(level)
            {}
            ~formatted_log_t()
            {
                // GLOBAL_LEVEL is a global variable and could be changed at runtime
                // customization could be here.
                if (level <= GLOBAL_VERBOSITY_LEVEL) cout << fmt << endl;
            }
            template <typename T>
            formatted_log_t& operator %(T value) {
                fmt % value;
                return *this;
            }
        protected:
            log_level_t level;
            boost::format fmt;
    };
}
//namespace log_impl
// Helper function. Class formatted_log_t will not be used directly
// This function needs to be called eg log<LOG_BASIC>(L"Basic result is %1 divided by %2.") % 1 % L"Two";
template <log_level_t level>
log_impl::formatted_log_t log(const char* msg) {
    return log_impl::formatted_log_t(level, msg);
}

#endif
