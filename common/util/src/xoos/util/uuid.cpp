#include "xoos/util/uuid.h"

#include <uuid.h>

#include <fmt/chrono.h>
#include <fmt/format.h>

namespace xoos::uuid {

uuids::uuid CreateRandomUuid() {
  std::random_device rd;
  auto seed_data = std::array<int, std::mt19937::state_size>{};
  std::generate(std::begin(seed_data), std::end(seed_data), std::ref(rd));
  std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
  std::mt19937 generator(seq);
  uuids::basic_uuid_random_generator gen{generator};
  return gen();
}

std::vector<uuids::uuid> CreateRandomUuids(uint32_t count) {
  std::random_device rd;
  auto seed_data = std::array<int, std::mt19937::state_size>{};
  std::generate(std::begin(seed_data), std::end(seed_data), std::ref(rd));
  std::seed_seq seq(std::begin(seed_data), std::end(seed_data));
  std::mt19937 generator(seq);
  const uuids::basic_uuid_random_generator gen{generator};
  auto vecs = std::vector<uuids::uuid>(count);
  std::generate(vecs.begin(), vecs.end(), gen);
  return vecs;
}

uuids::uuid CreateUuidFromString(const std::string& from_str) {
  uuids::uuid_name_generator gen(CreateRandomUuid());
  return gen(from_str);
}

uuids::uuid CreateUuidFromString(const std::string& from_str, uuids::uuid uuid_base) {
  uuids::uuid_name_generator gen(uuid_base);
  return gen(from_str);
}

uuids::uuid CreateUuidFromString(const std::string& from_str, const std::string& uuid_base) {
  return CreateUuidFromString(from_str, uuids::uuid::from_string(uuid_base).value());
}

static const UUID EXECUTION_ID = UUID();  // NOLINT

UUID ExecutionId() {
  return EXECUTION_ID;
}

std::string UUID::IsoDateString() const {
  return fmt::format("{:%Y%m%d}-{}", timestamp, uuids::to_string(uuid));
}

std::string UUID::String() const {
  return uuids::to_string(uuid);
}

std::string UUID::IsoDateTimeString() const {
  return fmt::format("{:%Y%m%d%H%M%S}-{}", timestamp, uuids::to_string(uuid));
}

UUID::UUID(const std::string& from_str, uuids::uuid base_uuid)
    : uuid(CreateUuidFromString(from_str, base_uuid)), timestamp(std::chrono::system_clock::now()) {
}

UUID::UUID(const std::string& from_str, const std::string& base_uuid_str)
    : uuid(CreateUuidFromString(from_str, base_uuid_str)), timestamp(std::chrono::system_clock::now()) {
}

UUID::UUID() : uuid(CreateRandomUuid()), timestamp(std::chrono::system_clock::now()) {
}
}  // namespace xoos::uuid
