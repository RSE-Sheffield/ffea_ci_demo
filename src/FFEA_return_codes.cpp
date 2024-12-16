// 
//  This file is part of the FFEA simulation package
//  
//  Copyright (c) by the Theory and Development FFEA teams,
//  as they appear in the README.md file. 
// 
//  FFEA is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
// 
//  FFEA is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU General Public License for more details.
// 
//  You should have received a copy of the GNU General Public License
//  along with FFEA.  If not, see <http://www.gnu.org/licenses/>.
// 
//  To help us fund FFEA development, we humbly ask that you cite 
//  the research papers on the package.
//

#include "FFEA_return_codes.h"

void FFEA_error_text() {
    printf("\x1b[31mERROR: \x1b[m");
}

void FFEA_caution_text() {
    printf("\x1b[33mCAUTION: \x1b[m");
}

FFEAException::FFEAException(const char* format, ...) {
    va_list argp;
    va_start(argp, format);
    err_message += parseArgs(format, argp);
    va_end(argp);
    // Print the exception message to stdout preceded by ERROR in ANSI Red
    printf("\x1b[33mERROR\x1b[m: %s\n", err_message.c_str());
}
FFEAException::FFEAException(const std::string& message)
    : err_message(message) { }

char const* FFEAException::what() const noexcept {
    return err_message.c_str();
}
std::string FFEAException::parseArgs(const char* format, va_list argp) {
    std::string rtn = format;
    // Create a copy of the va_list, as vsnprintf can invalidate elements of argp and find the required buffer length
    va_list argpCopy;
    va_copy(argpCopy, argp);
    const int buffLen = vsnprintf(nullptr, 0, format, argpCopy) + 1;
    va_end(argpCopy);
    char* buffer = reinterpret_cast<char*>(malloc(buffLen * sizeof(char)));
    // Populate the buffer with the original va_list
    int ct = vsnprintf(buffer, buffLen, format, argp);
    if (ct >= 0) {
        // Success!
        buffer[buffLen - 1] = '\0';
        rtn = std::string(buffer);
    }
    free(buffer);
    return rtn;
}