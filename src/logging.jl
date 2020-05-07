# Logging helping functions based on ConsoleLogger
# It's goal is to hide the lines where the log is coming from.

function ns_mfmt(level, _module, group, id, file, line)
    color = default_logcolor(level)
    prefix = (level == Logging.Warn ? "Warning" : string(level))*':'
    suffix = ""
    return color, prefix, suffix
end

nosource_debuglogger() = ConsoleLogger(stderr, Logging.Debug, meta_formatter=ns_mfmt)
nosource_infologger() = ConsoleLogger(stderr, Logging.Info, meta_formatter=ns_mfmt)
sourced_debuglogger() = ConsoleLogger(stderr, Logging.Debug)
sourced_infologger() = ConsoleLogger(stderr, Logging.Info)

# custom LogLevels

const IterInf = Logging.Info-5
const IterLow1 = Logging.Info-10
const IterLow2 = Logging.Info-15
const Compute1 = Logging.Info-20
const Compute2 = Logging.Info-25

# custom loggers for custom Loglevels

nosource_IterInflog() = ConsoleLogger(stderr, IterInf, meta_formatter=ns_mfmt)
nosource_IterLow1log() = ConsoleLogger(stderr, IterLow1, meta_formatter=ns_mfmt)
nosource_IterLow2log() = ConsoleLogger(stderr, IterLow2, meta_formatter=ns_mfmt)
nosource_Compute1log() = ConsoleLogger(stderr, Compute1, meta_formatter=ns_mfmt)
nosource_Compute2log() = ConsoleLogger(stderr, Compute2, meta_formatter=ns_mfmt)
nosource_Errorlog() = ConsoleLogger(stderr, Logging.Error, meta_formatter=ns_mfmt)
