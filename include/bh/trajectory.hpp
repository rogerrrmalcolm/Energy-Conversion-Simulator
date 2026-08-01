#pragma once

#include <cstddef>
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

struct IntegrationDiagnostics {
    std::size_t accepted_steps{};
    std::size_t rejected_steps{};
    // Dimensionless residual of the radial first integral at reported points.
    double maximum_normalized_radial_residual{};
    double maximum_normalized_error{};
    double final_step{};
};

struct Trajectory {
    std::vector<TrajectoryPoint> points;
    TrajectoryTermination termination{TrajectoryTermination::completed};
    IntegrationDiagnostics diagnostics{};
};
}
