#include <core/RandomNumberGenerator.h>
#include <utils/Options.h>

using namespace Minisat;

static const char* _cat = "CORE";
static DoubleOption opt_random_seed (_cat, "rnd-seed", "Used by the random variable selection", 91648253, DoubleRange(0, false, HUGE_VAL, false));

RandomNumberGenerator::RandomNumberGenerator()
    : random_seed(opt_random_seed)
{}