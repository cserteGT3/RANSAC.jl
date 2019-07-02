#-------------------------------------------------------------------------------
# RANSACLogger based on SimpleLogger
# It's goal is to hide the lines where the log is coming.

"""
    RANSACLogger(stream::IO=stderr, level=Logging.Info)

A simple logger that doesn't prints the module, line number, etc.
Otherwise it behaves like [SimpleLogger](@ref).
"""
struct RANSACLogger <: AbstractLogger
    stream::IO
    min_level::LogLevel
    message_limits::Dict{Any,Int}
end

RANSACLogger(stream::IO=stderr, level=Logging.Info) = RANSACLogger(stream, level, Dict{Any,Int}())

shouldlog(logger::RANSACLogger, level, _module, group, id) = get(logger.message_limits, id, 1) > 0

min_enabled_level(logger::RANSACLogger) = logger.min_level

catch_exceptions(logger::RANSACLogger) = false

function handle_message(logger::RANSACLogger, level, message, _module, group, id,
                        filepath, line; maxlog=nothing, kwargs...)
    if maxlog != nothing && maxlog isa Integer
        remaining = get!(logger.message_limits, id, maxlog)
        logger.message_limits[id] = remaining - 1
        remaining > 0 || return
    end
    buf = IOBuffer()
    iob = IOContext(buf, logger.stream)
    levelstr = level == Warn ? "Warning" : string(level)
    msglines = split(chomp(string(message)), '\n')
    println(iob, "┌ ", levelstr, ": ", msglines[1])
    for i in 2:length(msglines)
        println(iob, "│ ", msglines[i])
    end
    for (key, val) in kwargs
        println(iob, "│   ", key, " = ", val)
    end
    write(logger.stream, take!(buf))
    nothing
end

ransacdebuglogger() = RANSACLogger(stderr, Logging.Debug)
ransacinfologger() = RANSACLogger(stderr, Logging.Info)
