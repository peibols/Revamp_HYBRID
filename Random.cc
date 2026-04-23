#include "Random.h"

numrand::numrand()
    : _ir(0), rng_(std::random_device{}()), dist_(0.0, 1.0)
{
}

numrand::numrand(int ir)
    : _ir(ir), rng_(static_cast<std::uint64_t>(ir)), dist_(0.0, 1.0)
{
}

numrand::~numrand() = default;

double numrand::rando()
{
  return dist_(rng_);
}

void numrand::SetIr(int ir)
{
  _ir = ir;
  rng_.seed(static_cast<std::uint64_t>(ir));
}

int numrand::GetIr() const
{
  return _ir;
}
