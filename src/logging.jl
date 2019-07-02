# Logging helping functions based on ConsoleLogger
# It's goal is to hide the lines where the log is coming from.

function nosource_metafmt(level, _module, group, id, file, line)
    color = default_logcolor(level)
    prefix = (level == Logging.Warn ? "Warning" : string(level))*':'
    suffix = ""
    return color, prefix, suffix
end

nosource_debuglogger() = ConsoleLogger(stderr, Logging.Debug, meta_formatter=nosource_metafmt)
nosource_infologger() = ConsoleLogger(stderr, Logging.Info, meta_formatter=nosource_metafmt)
sourced_debuglogger() = ConsoleLogger(stderr, Logging.Debug)
sourced_infologger() = ConsoleLogger(stderr, Logging.Info)
