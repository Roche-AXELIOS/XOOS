#pragma once

#include <uuid.h>

#include <string>
#include <vector>

namespace xoos::uuid {

/**
 * Create a UUID from a random generator (uuid v4)
 */
uuids::uuid CreateRandomUuid();
/**
 * Return a series of uuids from a random generator (uuid v4)
 */
std::vector<uuids::uuid> CreateRandomUuids(uint32_t count);

// Create a UUID from a string using the name generator (uuid v5)

/**
 * Create a UUID from a string using the name generator (uuid v5) based on a random UUID.
 */
uuids::uuid CreateUuidFromString(const std::string& from_str);

/**
 * Create a UUID from a string using the name generator and base UUID (uuid v5)
 */
uuids::uuid CreateUuidFromString(const std::string& from_str, uuids::uuid uuid_base);

/**
 * Create a UUID from a string using the name generator and base UUID as a string (uuid v5)
 */
uuids::uuid CreateUuidFromString(const std::string& from_str, const std::string& uuid_base);

/**
 * UUID struct that contains the UUID and a timestamp
 */
struct UUID {
  uuids::uuid uuid;
  std::chrono::time_point<std::chrono::system_clock> timestamp;

  UUID();

  UUID(const std::string& from_str, const std::string& base_uuid_str);

  UUID(const std::string& from_str, uuids::uuid base_uuid);

  /**
   * Returns the UUID string
   */
  std::string String() const;

  /**
   * Returns the UUID string with a prefixed iso date and time.
   */
  std::string IsoDateTimeString() const;

  /**
   * Returns the UUID string with a prefixed iso date
   */
  std::string IsoDateString() const;
};

UUID ExecutionId();

}  // namespace xoos::uuid
