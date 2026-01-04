#pragma once

#include <string>
#include "DATA_SOURCE.h"

namespace ridgelets {

struct ParseResult {
    bool ok;
    DATA_SOURCE::input_parse params;
    std::string error_message;
};

/**
 * Parse command-line arguments and fill DATA_SOURCE::input_parse.
 *
 * argc/argv are the standard program arguments.
 * Returns a ParseResult with ok==true on success, or ok==false and
 * error_message set on parse/validation error. If user requested help,
 * ok==false and error_message will be empty (help printed by parser).
 */
ParseResult parse_cli(int argc, char** argv);

} // namespace ridgelets
