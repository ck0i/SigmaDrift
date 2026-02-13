#pragma once

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <numbers>
#include <numeric>

namespace motor_synergy {

struct trajectory_point {
    double x, y;
    double t;
};

struct config {
    double fitts_a = 50.0;
    double fitts_b = 150.0;
    double target_width = 20.0;

    double undershoot_min = 0.92;
    double undershoot_max = 0.97;
    double peak_time_ratio = 0.35;
    double primary_sigma_min = 0.18;
    double primary_sigma_max = 0.28;

    double overshoot_prob = 0.15;
    double overshoot_min = 1.02;
    double overshoot_max = 1.08;
    double correction_sigma_min = 0.12;
    double correction_sigma_max = 0.20;
    double second_correction_prob = 0.25;

    double curvature_scale = 0.025;

    double ou_theta = 3.5;
    double ou_sigma = 1.2;

    double tremor_freq_min = 8.0;
    double tremor_freq_max = 12.0;
    double tremor_amp_min = 0.15;
    double tremor_amp_max = 0.55;

    double sdn_k = 0.04;

    double sample_dt_mean = 7.8;
    double gamma_shape = 3.5;
};

struct metrics {
    double movement_time;
    double path_length;
    double straight_distance;
    double path_efficiency;
    double peak_speed;
    double time_to_peak;
    int num_submovements;
    double endpoint_error;
    double fitts_predicted_mt;
};

namespace detail {

inline double normal_cdf(double x) {
    return 0.5 * (1.0 + std::erf(x / std::numbers::sqrt2));
}

inline double lognormal_cdf(double t, double t0, double mu, double sigma) {
    if (t <= t0) return 0.0;
    return normal_cdf((std::log(t - t0) - mu) / sigma);
}

inline double lognormal_pdf(double t, double t0, double mu, double sigma) {
    if (t <= t0) return 0.0;
    double dt = t - t0;
    double z = (std::log(dt) - mu) / sigma;
    return 1.0 / (sigma * std::sqrt(2.0 * std::numbers::pi) * dt)
         * std::exp(-0.5 * z * z);
}

// s^2*(1-s)^3 normalized to peak=1.0 at s=0.4 - can be tweaked.
// curvature is maximal during the acceleration phase
inline double curvature_profile(double s) {
    if (s <= 0.0 || s >= 1.0) return 0.0;
    double v = s * s * (1.0 - s) * (1.0 - s) * (1.0 - s);
    constexpr double norm = 0.4 * 0.4 * 0.6 * 0.6 * 0.6;
    return v / norm;
}

// vertical movements produce more curvature due to wrist/forearm geometry
inline double direction_factor(double angle) {
    double sa = std::abs(std::sin(angle));
    double ca = std::abs(std::cos(angle));
    return 0.5 + 0.8 * sa - 0.15 * ca;
}

} // namespace detail

inline std::vector<trajectory_point> generate(
    double x0, double y0, double x1, double y1,
    const config& cfg = {}, uint64_t seed = 0)
{
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    auto uniform = [&](double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(rng);
    };
    auto normal = [&](double m, double s) {
        return std::normal_distribution<double>(m, s)(rng);
    };
    auto gamma_dist = [&](double shape, double scale) {
        return std::gamma_distribution<double>(shape, scale)(rng);
    };

    double dx = x1 - x0, dy = y1 - y0;
    double distance = std::hypot(dx, dy);
    double direction = std::atan2(dy, dx);

    if (distance < 1.0)
        return {{x0, y0, 0.0}, {x1, y1, 50.0}};

    double tx = dx / distance, ty = dy / distance;
    double nx = -ty, ny = tx;

    double id = std::log2(distance / cfg.target_width + 1.0);
    double mt = (cfg.fitts_a + cfg.fitts_b * id) * std::exp(normal(0.0, 0.08));
    mt = std::max(mt, 80.0);

    bool overshoot = uniform(0.0, 1.0) < cfg.overshoot_prob;
    double reach = overshoot
        ? uniform(cfg.overshoot_min, cfg.overshoot_max)
        : uniform(cfg.undershoot_min, cfg.undershoot_max);

    double primary_D = distance * reach;
    double primary_sigma = uniform(cfg.primary_sigma_min, cfg.primary_sigma_max);

    // mu derived from mode = exp(mu - sigma^2) so that peak velocity lands at peak_t
    double peak_t = mt * uniform(cfg.peak_time_ratio - 0.03, cfg.peak_time_ratio + 0.03);
    double primary_mu = std::log(peak_t) + primary_sigma * primary_sigma;

    struct correction { double D, t0, mu, sigma, dir_x, dir_y; };
    std::vector<correction> corrections;

    double remaining = distance - primary_D;
    if (std::abs(remaining) > 0.5) {
        double dir = remaining > 0.0 ? 1.0 : -1.0;
        double cD = std::abs(remaining) * uniform(0.88, 1.02);
        double cS = uniform(cfg.correction_sigma_min, cfg.correction_sigma_max);
        double cPeak = mt * uniform(0.12, 0.18);
        corrections.push_back({
            cD, mt * uniform(0.55, 0.68),
            std::log(cPeak) + cS * cS, cS,
            tx * dir, ty * dir
        });

        double left = remaining - cD * dir;
        if (std::abs(left) > 0.3 && uniform(0.0, 1.0) < cfg.second_correction_prob) {
            double d2 = left > 0.0 ? 1.0 : -1.0;
            double cD2 = std::abs(left) * uniform(0.85, 1.05);
            double cS2 = uniform(0.10, 0.16);
            double cP2 = mt * uniform(0.08, 0.12);
            corrections.push_back({
                cD2, mt * uniform(0.78, 0.88),
                std::log(cP2) + cS2 * cS2, cS2,
                tx * d2, ty * d2
            });
        }
    }

    double curv_amp = distance * cfg.curvature_scale
        * detail::direction_factor(direction) * normal(0.0, 1.0);

    double tremor_freq = uniform(cfg.tremor_freq_min, cfg.tremor_freq_max);
    double tremor_amp = uniform(cfg.tremor_amp_min, cfg.tremor_amp_max);
    double tph_x = uniform(0.0, 2.0 * std::numbers::pi);
    double tph_y = uniform(0.0, 2.0 * std::numbers::pi);
    double ou_x = 0.0, ou_y = 0.0;

    double total_t = mt * 1.15;
    double g_scale = cfg.sample_dt_mean / cfg.gamma_shape;

    std::vector<double> times = {0.0};
    for (double t = 0.0; t < total_t;) {
        double dt = std::clamp(gamma_dist(cfg.gamma_shape, g_scale), 2.0, 25.0);
        t += dt;
        if (t <= total_t + 15.0) times.push_back(t);
    }

    std::vector<trajectory_point> result;
    result.reserve(times.size());

    for (size_t i = 0; i < times.size(); ++i) {
        double t = times[i];
        double dt_ms = (i > 0) ? (t - times[i - 1]) : cfg.sample_dt_mean;
        double dt_s = dt_ms / 1000.0;

        double s = detail::lognormal_cdf(t, 0.0, primary_mu, primary_sigma);

        double bx = x0 + tx * primary_D * s;
        double by = y0 + ty * primary_D * s;

        bx += nx * curv_amp * detail::curvature_profile(s);
        by += ny * curv_amp * detail::curvature_profile(s);

        for (const auto& c : corrections) {
            double cs = detail::lognormal_cdf(t, c.t0, c.mu, c.sigma);
            bx += c.dir_x * c.D * cs;
            by += c.dir_y * c.D * cs;
        }

        double speed = primary_D * detail::lognormal_pdf(t, 0.0, primary_mu, primary_sigma);
        for (const auto& c : corrections)
            speed += c.D * detail::lognormal_pdf(t, c.t0, c.mu, c.sigma);

        // euler-maruyama step for OU process (dt_s in seconds)
        ou_x += -cfg.ou_theta * ou_x * dt_s + cfg.ou_sigma * std::sqrt(dt_s) * normal(0.0, 1.0);
        ou_y += -cfg.ou_theta * ou_y * dt_s + cfg.ou_sigma * std::sqrt(dt_s) * normal(0.0, 1.0);

        // tremor gain drops with speed (proprioceptive suppression)
        double t_s = t / 1000.0;
        double trem_mod = 1.0 / (1.0 + speed * 0.3);
        double tr_x = tremor_amp * trem_mod
            * std::sin(2.0 * std::numbers::pi * tremor_freq * t_s + tph_x);
        double tr_y = tremor_amp * trem_mod
            * std::sin(2.0 * std::numbers::pi * tremor_freq * t_s + tph_y);

        // noise magnitude proportional to motor command (harris-wolpert SDN)
        double sdn_x = cfg.sdn_k * speed * normal(0.0, 1.0);
        double sdn_y = cfg.sdn_k * speed * normal(0.0, 1.0);

        result.push_back({bx + ou_x + tr_x + sdn_x, by + ou_y + tr_y + sdn_y, t});
    }

    return result;
}

inline metrics compute_metrics(
    const std::vector<trajectory_point>& path,
    double target_x, double target_y,
    double target_width, double straight_dist)
{
    metrics m{};
    if (path.size() < 2) return m;

    m.movement_time = path.back().t - path.front().t;
    m.straight_distance = straight_dist;

    double max_speed = 0.0;
    m.path_length = 0.0;
    std::vector<double> speeds(path.size(), 0.0);

    for (size_t i = 1; i < path.size(); ++i) {
        double dx = path[i].x - path[i - 1].x;
        double dy = path[i].y - path[i - 1].y;
        double dt = path[i].t - path[i - 1].t;
        double seg = std::hypot(dx, dy);
        m.path_length += seg;
        double spd = (dt > 0.0) ? seg / dt : 0.0;
        speeds[i] = spd;
        if (spd > max_speed) { max_speed = spd; m.time_to_peak = path[i].t; }
    }

    m.peak_speed = max_speed;
    m.path_efficiency = (m.path_length > 0.0) ? m.straight_distance / m.path_length : 1.0;

    // peaks above 15% of max in the speed signal → sub-movement count
    double threshold = max_speed * 0.15;
    int peaks = 0;
    for (size_t i = 2; i + 1 < speeds.size(); ++i) {
        if (speeds[i] > threshold && speeds[i] > speeds[i - 1] && speeds[i] > speeds[i + 1])
            ++peaks;
    }
    m.num_submovements = std::max(peaks, 1);

    m.endpoint_error = std::hypot(path.back().x - target_x, path.back().y - target_y);

    double id = std::log2(straight_dist / target_width + 1.0);
    m.fitts_predicted_mt = 50.0 + 150.0 * id;

    return m;
}

} // namespace motor_synergy

namespace windmouse { 
    // this is windmouse, for comparison. the graph will show both trajectories and metrics side by side,
	// but the main point is to show how windmouse's trajectory is more erratic and less efficient than sigma-drift's, 
    // due to the lack of an explicit internal model and reliance on noisy feedback and random exploration.
	// - this leads to detections in the metrics like lower path efficiency, more sub-movements, and higher endpoint error, 
    // especially on smaller targets.

inline std::vector<motor_synergy::trajectory_point> generate(
    double x0, double y0, double x1, double y1,
    double gravity = 9.0, double wind_str = 3.0,
    double max_step = 15.0, double target_area = 8.0,
    uint64_t seed = 0)
{
    std::mt19937_64 rng(seed ? seed : std::random_device{}());
    auto randf = [&](double lo, double hi) {
        return std::uniform_real_distribution<double>(lo, hi)(rng);
    };

    std::vector<motor_synergy::trajectory_point> result;
    double xs = x0, ys = y0;
    double vx = 0, vy = 0, wx = 0, wy = 0;
    double t = 0.0, step = max_step;

    result.push_back({xs, ys, 0.0});

    for (int iter = 0; iter < 5000; ++iter) {
        double dist = std::hypot(x1 - xs, y1 - ys);
        if (dist < 1.0) break;

        double w = std::min(wind_str, dist);
        if (dist >= target_area) {
            wx = wx / std::sqrt(3.0) + randf(-w, w) / std::sqrt(5.0);
            wy = wy / std::sqrt(3.0) + randf(-w, w) / std::sqrt(5.0);
        } else {
            wx /= std::sqrt(3.0);
            wy /= std::sqrt(3.0);
            if (step < 3.0) step = randf(3.0, 6.0);
            else step /= std::sqrt(5.0);
        }

        vx += wx + gravity * (x1 - xs) / dist;
        vy += wy + gravity * (y1 - ys) / dist;
        double vmag = std::hypot(vx, vy);
        if (vmag > step) {
            double r = step / 2.0 + randf(0.0, step / 2.0);
            vx = vx / vmag * r;
            vy = vy / vmag * r;
        }

        xs += vx;
        ys += vy;
        t += randf(5.0, 15.0);
        result.push_back({xs, ys, t});
    }

    return result;
}

} // namespace windmouse
