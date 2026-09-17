#pragma once

#include <boost/version.hpp>

#if BOOST_VERSION <= 108000
#  define BOOST_USE_WINAPI_VERSION BOOST_WINAPI_VERSION_WIN7
#endif




using namespace std::literals;
namespace logging = boost::log;
namespace sinks = boost::log::sinks;
namespace keywords = boost::log::keywords;
namespace expr = boost::log ::expressions;
namespace attrs = boost::log::attributes;

#define __FILENAME__ \
  (strrchr(__FILE__, '\\') ? strrchr(__FILE__, '\\') + 1 : __FILE__)

BOOST_LOG_ATTRIBUTE_KEYWORD(line_id, "LineID", unsigned int)
BOOST_LOG_ATTRIBUTE_KEYWORD(line, "Line", int)
BOOST_LOG_ATTRIBUTE_KEYWORD(file, "File", std::string)
BOOST_LOG_ATTRIBUTE_KEYWORD(timestamp, "TimeStamp", boost::posix_time::ptime)



void MyFormatter(logging::record_view const& rec,
                 logging::formatting_ostream& strm);

void InitBoostLogFilter(const std::string& file_name,
                        boost::log::trivial::severity_level level);

#define BLOG(LEVEL)                                                       \
  BOOST_LOG_TRIVIAL(LEVEL) << boost::log::add_value("Line", __LINE__)     \
                           << boost::log::add_value("File", __FILENAME__) \
                           << boost::log::add_value("Function", __FUNCTION__)


