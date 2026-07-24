#pragma once

#include <vector>

namespace bh {
struct TrajectoryPoint {
    double affine_parameter{};
    double radius{};
    double phi{};
    double coordinate_time{};
    double radial_velocity{};
};

enum class TrajectoryTermination {
    completed,
    crossed_horizon,
    reached_escape_radius,
    reached_target_radius,
    turning_point,
    invalid_state
};

struct Trajectory {
    std::vector<TrajectoryPoint> points;
    TrajectoryTermination termination{TrajectoryTermination::completed};
};
}
