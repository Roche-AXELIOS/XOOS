#include <xoos/cli/cli.h>

#include "alignment-metrics-cli.h"

int main(int argc, char** argv) {
  return xoos::alignment_metrics::RunAlignmentMetricsApp(argc, argv, PROGRAM_NAME, VERSION);
}
