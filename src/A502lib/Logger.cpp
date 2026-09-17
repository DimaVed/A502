
#include "Logger.h"

void MyFormatter(logging::record_view const& rec,
                 logging::formatting_ostream& strm) {
  auto ts = *rec[timestamp];
  // Output time
  strm << to_iso_extended_string(ts) << " ";

  // Output severity
  strm << rec[logging::trivial::severity] << " ";
  // output file and line
  strm << "[" << rec[file] << " " << rec[line] << "] ";

  // Output message
  strm << rec[logging::expressions::smessage];
}

void InitBoostLogFilter(const std::string& file_name, 
                        boost::log::trivial::severity_level level) {

      // === НАЧАЛО ПРОВЕРКИ ДИСКА ===
  try {
    std::filesystem::path log_path(file_name);
    // Получаем директорию (если в file_name передан только файл, это будет текущая папка)
    std::filesystem::path log_dir = log_path.parent_path();
    if (log_dir.empty()) {
      log_dir = std::filesystem::current_path();
    }

    // 1. Создаем директорию, если её нет
    if (!std::filesystem::exists(log_dir)) {
      std::filesystem::create_directories(log_dir);
    }

    // 2. Тест на запись: создаем временный проверочный файл
    std::filesystem::path test_file = log_dir / ".log_write_test";
    std::ofstream         test_stream(test_file, std::ios::out | std::ios::app);

    if (!test_stream.is_open()) {
      throw std::runtime_error("Директория недоступна для записи (отказано в доступе)");
    }

    test_stream << "test";
    test_stream.close();
    std::filesystem::remove(test_file); // Удаляем тестовый файл

  } catch (const std::exception& e) {
    // Обработка ошибки: пишем в консоль и прекращаем инициализацию файлового лога
    std::cerr << "[CRITICAL LOG ERROR] Не удалось инициализировать логи на диске: " << e.what()
              << std::endl;

    // Опционально: инициализируем ТОЛЬКО консоль, чтобы приложение не падало вслепую
    logging::add_common_attributes();
    boost::log::core::get()->set_filter(boost::log::trivial::severity >= level);
    logging::add_console_log(std::cout, keywords::format = &MyFormatter);
    return;
  }
  // === КОНЕЦ ПРОВЕРКИ ДИСКА ===


  auto log_file_attribute = file_name + "%5N.log";
  logging::add_common_attributes();
  
  boost::shared_ptr<logging::core> core = logging::core::get();
  core->set_filter(boost::log::trivial::severity >= level);
  logging::add_console_log(std::cout, keywords::format = &MyFormatter);
  boost::shared_ptr<sinks::text_file_backend> backend =
      boost::make_shared<sinks::text_file_backend>(
          // file name pattern
          keywords::file_name = log_file_attribute.data(),
          // rotate the file upon reaching 50 MiB size...
          keywords::rotation_size = 50ull * 1024ull * 1024ull,
          // ...or at noon, whichever comes first
          keywords::time_based_rotation =
              sinks::file::rotation_at_time_point(12, 0, 0),
          // append old log file
          keywords::open_mode = std::ios_base::out | std::ios_base::app,
          // flush each string
          keywords::auto_flush = true);

  typedef sinks::synchronous_sink<sinks::text_file_backend> SinkT;
  boost::shared_ptr<SinkT> sink(new SinkT(backend));
  sink->set_formatter(&MyFormatter);
  core->add_sink(sink);
}

