#include "Logger.h"

#include <filesystem>
#include <fstream>
#include <iostream>

void MyFormatter(logging::record_view const& rec,
                 logging::formatting_ostream& strm) {
  auto ts = *rec[timestamp];
  strm << to_iso_extended_string(ts) << " ";
  strm << rec[logging::trivial::severity] << " ";
  strm << "[" << rec[file] << " " << rec[line] << "] ";
  strm << rec[logging::expressions::smessage];
}

void InitBoostLogFilter(const std::string& file_name,
                        boost::log::trivial::severity_level level) {
  try {
    std::filesystem::path log_path(file_name);
    std::filesystem::path log_dir = log_path.parent_path();
    if (log_dir.empty()) {
      log_dir = std::filesystem::current_path();
    }

    if (!std::filesystem::exists(log_dir)) {
      std::filesystem::create_directories(log_dir);
    }

    std::filesystem::path test_file = log_dir / ".log_write_test";
    std::ofstream test_stream(test_file, std::ios::out | std::ios::app);

    if (!test_stream.is_open()) {
      throw std::runtime_error(
          "Директория недоступна для записи (отказано в доступе)");
    }

    test_stream << "test";
    test_stream.close();
    std::filesystem::remove(test_file);
  } catch (const std::exception& e) {
    std::cerr << "[CRITICAL LOG ERROR] Не удалось инициализировать логи на "
                 "диске: "
              << e.what() << std::endl;

    logging::add_common_attributes();
    boost::log::core::get()->set_filter(boost::log::trivial::severity >=
                                        level);
    logging::add_console_log(std::cout, keywords::format = &MyFormatter);
    return;
  }

  auto log_file_attribute = file_name + "%5N.log";
  logging::add_common_attributes();

  boost::shared_ptr<logging::core> core = logging::core::get();
  core->set_filter(boost::log::trivial::severity >= level);
  logging::add_console_log(std::cout, keywords::format = &MyFormatter);
  boost::shared_ptr<sinks::text_file_backend> backend =
      boost::make_shared<sinks::text_file_backend>(
          keywords::file_name = log_file_attribute.data(),
          keywords::rotation_size = 50ull * 1024ull * 1024ull,
          keywords::time_based_rotation =
              sinks::file::rotation_at_time_point(12, 0, 0),
          keywords::open_mode = std::ios_base::out | std::ios_base::app,
          keywords::auto_flush = true);

  typedef sinks::synchronous_sink<sinks::text_file_backend> SinkT;
  boost::shared_ptr<SinkT> sink(new SinkT(backend));
  sink->set_formatter(&MyFormatter);
  core->add_sink(sink);
}
