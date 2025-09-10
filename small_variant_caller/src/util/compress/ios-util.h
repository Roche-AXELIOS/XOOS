#pragma once

#include <istream>

#include <htslib/hts.h>

#include <xoos/types/fs.h>
#include <xoos/types/vec.h>

/**
 * @file ios-util.h
 * @brief Input/output stream utilities for file operations.
 *
 * This module provides convenient wrapper functions for opening files,
 * checking stream states, and performing buffered I/O operations.
 */

namespace xoos::svc {

/**
 * @brief Opens an output file stream with error checking.
 *
 * Creates and opens an output file stream with the specified mode.
 * Throws an exception if the file cannot be opened.
 *
 * @param output Path to the output file to open
 * @param mode Opening mode flags (e.g., std::ios::binary, std::ios::app)
 * @return An opened output file stream
 *
 * @throws std::runtime_error if the file cannot be opened for writing
 *
 * @note The parent directory must exist before calling this function.
 *       The function will create or overwrite the file depending on the mode.
 */
std::ofstream OpenOutput(const fs::path& output, std::ios::openmode mode);

/**
 * @brief Opens an output file stream with default mode.
 *
 * Creates and opens an output file stream with default flags (out | trunc).
 * Throws an exception if the file cannot be opened.
 *
 * @param output Path to the output file to open
 * @return An opened output file stream
 *
 * @throws std::runtime_error if the file cannot be opened for writing
 *
 * @note Uses default mode std::ios::out | std::ios::trunc which creates
 *       a new file or truncates an existing file.
 */
std::ofstream OpenOutput(const fs::path& output);

/**
 * @brief Opens an input file stream with error checking.
 *
 * Creates and opens an input file stream with the specified mode.
 * Throws an exception if the file cannot be opened.
 *
 * @param input Path to the input file to open
 * @param mode Opening mode flags (e.g., std::ios::binary)
 * @return An opened input file stream
 *
 * @throws std::runtime_error if the file cannot be opened for reading
 *
 * @note The file must exist and be readable before calling this function.
 */
std::ifstream OpenInput(const fs::path& input, std::ios::openmode mode);

/**
 * @brief Opens an input file stream with default mode.
 *
 * Creates and opens an input file stream with default flags (std::ios::in).
 * Throws an exception if the file cannot be opened.
 *
 * @param input Path to the input file to open
 * @return An opened input file stream
 *
 * @throws std::runtime_error if the file cannot be opened for reading
 *
 * @note The file must exist and be readable before calling this function.
 */
std::ifstream OpenInput(const fs::path& input);

/**
 * @brief Checks if an I/O stream is in a bad state and throws if so.
 *
 * Verifies that the stream is not in a bad state (indicating serious I/O errors).
 * Throws an exception with the file path if the stream has encountered errors.
 *
 * @param path The file path associated with the stream (for error reporting)
 * @param ios The I/O stream to check
 *
 * @throws std::runtime_error if the stream is in a bad state
 *
 * @note This function only checks the bad() flag, not fail() or eof().
 *       Use after I/O operations to detect hardware/system-level errors.
 */
void RequireIoNotBad(const fs::path& path, std::ios& ios);

/**
 * @brief Checks if an I/O stream is in a bad state and throws if so.
 *
 * Verifies that the stream is not in a bad state (indicating serious I/O errors).
 * Throws an exception if the stream has encountered errors.
 *
 * @param ios The I/O stream to check
 *
 * @throws std::runtime_error if the stream is in a bad state
 *
 * @note This function only checks the bad() flag, not fail() or eof().
 *       Use after I/O operations to detect hardware/system-level errors.
 */
void RequireIoNotBad(std::ios& ios);

/**
 * @brief Reads data from a stream into a buffer.
 *
 * Attempts to read up to buffer.size() bytes from the input stream into
 * the provided buffer. Returns the actual number of bytes read.
 *
 * @param in The input stream to read from
 * @param buffer The buffer to read data into (must be pre-allocated)
 * @return The number of bytes actually read (may be less than buffer size)
 *
 * @note The buffer must be pre-allocated to the desired read size.
 *       The function may read fewer bytes than requested if EOF is reached
 *       or if the stream encounters an error.
 *       Check the stream state after calling to detect errors or EOF.
 */
size_t Read(std::istream& in, vec<char>& buffer);

/**
 * @brief Reads entire content from a stream as a string.
 *
 * Reads all available data from the input stream and returns it as a string.
 * Uses iterator-based approach for efficient reading.
 *
 * @param in The input stream to read from
 * @return The entire stream content as a string
 *
 * @throws std::runtime_error if the stream encounters an error during reading
 *
 * @note Loads entire content into memory. Use with caution for large streams.
 *       The stream position will be at EOF after successful completion.
 */
std::string Read(std::istream& in);

/**
 * @brief Reads entire content from a file with specified mode.
 *
 * Opens the file with the specified mode and reads all content as a string.
 * Automatically handles file opening and closing with error checking.
 *
 * @param path Path to the file to read
 * @param mode Opening mode flags (e.g., std::ios::binary)
 * @return The entire file content as a string
 *
 * @throws std::runtime_error if the file cannot be opened or read
 *
 * @note Loads entire file content into memory. Use with caution for large files.
 */
std::string Read(const fs::path& path, std::ios::openmode mode);

/**
 * @brief Reads entire content from a file with default mode.
 *
 * Opens the file in binary mode and reads all content as a string.
 * Automatically handles file opening and closing with error checking.
 *
 * @param path Path to the file to read
 * @return The entire file content as a string
 *
 * @throws std::runtime_error if the file cannot be opened or read
 *
 * @note Uses binary mode for reading. Loads entire file content into memory.
 */
std::string Read(const fs::path& path);

/**
 * @brief Reads lines from a stream into a vector.
 *
 * Reads the stream line by line and returns each line as a separate string
 * in a vector. Uses std::getline for line parsing.
 *
 * @param in The input stream to read from
 * @return Vector of strings, each containing one line from the stream
 *
 * @throws std::runtime_error if the stream encounters an error during reading
 *
 * @note Empty lines are preserved as empty strings in the vector.
 *       Line endings are removed by std::getline.
 */
vec<std::string> ReadLines(std::istream& in);

/**
 * @brief Reads lines from a file into a vector.
 *
 * Opens the file and reads it line by line, returning each line as a separate
 * string in a vector. Automatically handles file opening and closing.
 *
 * @param path Path to the file to read
 * @return Vector of strings, each containing one line from the file
 *
 * @throws std::runtime_error if the file cannot be opened or read
 *
 * @note Empty lines are preserved as empty strings in the vector.
 *       Uses text mode for reading to handle line endings properly.
 */
vec<std::string> ReadLines(const fs::path& path);

/**
 * @brief Writes string content to an output stream.
 *
 * Writes the provided string content to the output stream and checks
 * for errors after the operation.
 *
 * @param out The output stream to write to
 * @param content The string content to write
 *
 * @throws std::runtime_error if the stream encounters an error during writing
 *
 * @note Does not add any line endings. Content is written exactly as provided.
 */
void Write(std::ostream& out, const std::string& content);

/**
 * @brief Writes string content to a file with specified mode.
 *
 * Opens the file with the specified mode and writes the content.
 * Automatically handles file opening and closing with error checking.
 *
 * @param path Path to the file to write
 * @param mode Opening mode flags (e.g., std::ios::binary, std::ios::app)
 * @param content The string content to write
 *
 * @throws std::runtime_error if the file cannot be opened or written
 *
 * @note The parent directory must exist before calling this function.
 */
void Write(const fs::path& path, std::ios::openmode mode, const std::string& content);

/**
 * @brief Writes string content to a file with default mode.
 *
 * Opens the file with default mode (out | trunc) and writes the content.
 * Automatically handles file opening and closing with error checking.
 *
 * @param path Path to the file to write
 * @param content The string content to write
 *
 * @throws std::runtime_error if the file cannot be opened or written
 *
 * @note Creates a new file or truncates an existing file.
 *       The parent directory must exist before calling this function.
 */
void Write(const fs::path& path, const std::string& content);

/**
 * @brief Writes lines to an output stream.
 *
 * Writes each string in the vector as a separate line to the output stream.
 * Automatically adds newline characters after each line.
 *
 * @param out The output stream to write to
 * @param lines Vector of strings to write as lines
 *
 * @throws std::runtime_error if the stream encounters an error during writing
 *
 * @note Each string in the vector becomes one line with '\n' appended.
 *       Empty strings in the vector will result in blank lines.
 */
void WriteLines(std::ostream& out, const vec<std::string>& lines);

/**
 * @brief Writes lines to a file.
 *
 * Opens the file and writes each string in the vector as a separate line.
 * Automatically handles file opening and closing with error checking.
 *
 * @param path Path to the file to write
 * @param lines Vector of strings to write as lines
 *
 * @throws std::runtime_error if the file cannot be opened or written
 *
 * @note Creates a new file or truncates an existing file.
 *       Each string becomes one line with '\n' appended.
 */
void WriteLines(const fs::path& path, const vec<std::string>& lines);

}  // namespace xoos::svc
