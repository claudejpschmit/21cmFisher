#include "Log.hpp"

log_level_t GLOBAL_VERBOSITY_LEVEL = LOG_DEBUG;


void set_verbosity(log_level_t level)
{
    GLOBAL_VERBOSITY_LEVEL = level;
}


log_impl::formatted_log_t::formatted_log_t(log_level_t level, const wchar_t* msg) 
     :
        fmt(msg),
        level(level)
{}
log_impl::formatted_log_t::~formatted_log_t()
           {
                // GLOBAL_LEVEL is a global variable and could be changed at runtime
                // customization could be here.
                if (level <= GLOBAL_VERBOSITY_LEVEL) wcout << fmt << endl;
            }
template <log_level_t level>
log_impl::formatted_log_t log(const wchar_t* msg) {
    return log_impl::formatted_log_t(level, msg);
}

