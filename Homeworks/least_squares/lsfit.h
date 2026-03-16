#pragma once
#include "qr.h"
#include <functional>
#include <tuple>
#include <vector>

namespace pp {

using Func = std::function<double(double)>;

std::tuple<pp::vector, pp::matrix> lsfit(
    const std::vector<Func>& fs,
    const pp::vector& x,
    const pp::vector& y,
    const pp::vector& dy
);

} // namespace pp