#include "core/constants.hpp"
#include "core/droplet.hpp"
#include "io/binning.hpp"
#include "physics/collision.hpp"
#include "physics/kernels.hpp"
#include "utils/math.hpp"
#include "utils/time.hpp"

#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>

namespace {

constexpr bool DEBUG = true;
constexpr bool PLOT = true;

/// Print interval for status updates (in seconds)
constexpr unsigned int STATUS_DT = 60;

/// Write CSV header
void write_csv_header(std::ofstream &csv_file,
                      const sd_cpp::io::BinGrid &bin_grid) {
  csv_file << "-9999";
  for (size_t i = 0; i < sd_cpp::io::NR; ++i) {
    csv_file << "," << std::scientific << bin_grid.mids[i];
  }
  csv_file << "\n";
}

/// Write CSV row
void write_csv_row(std::ofstream &csv_file, unsigned int time,
                   const std::vector<double> &values) {
  csv_file << time;
  for (double val : values) {
    csv_file << "," << std::scientific << val;
  }
  csv_file << "\n";
}

/// Compute rolling median (smoothing)
std::vector<double>
smooth_values(const std::array<double, sd_cpp::io::NR> &values,
              size_t window_size) {
  std::vector<double> smoothed;
  smoothed.reserve(sd_cpp::io::NR);

  const size_t half_window = window_size / 2;

  for (size_t i = 0; i < sd_cpp::io::NR; ++i) {
    const size_t start = (i >= half_window) ? (i - half_window) : 0;
    const size_t end = std::min(i + half_window + 1, sd_cpp::io::NR);

    std::vector<double> window_vals(values.begin() + start,
                                    values.begin() + end);
    smoothed.push_back(sd_cpp::utils::median(window_vals));
  }

  return smoothed;
}

} // anonymous namespace

int main() {
  using namespace sd_cpp;

  // SIMULATION PARAMETERS
  const double n0 = 27.0 * utils::pow2(23); // Total number concentration
  const double r0 = 30.531e-6 / 3.0;        // Total droplet radius
  const double x0 =
      (4.0 / 3.0) * constants::PI * std::pow(r0, 3); // Total droplet volume
  const double delta_v = 1e6;                        // Total parcel volume
  const unsigned int t_c = 1;                        // Model timestep (seconds)
  const physics::Kernel kernel = physics::Kernel::Golovin; // Collision kernel

  const size_t n_part =
      static_cast<size_t>(utils::pow2(21)); // Total number of superdroplets
  const unsigned int t_end = 3600;          // Total simulation time (seconds)
  const unsigned int plot_dt = 600;         // Output interval time
  const size_t smooth_window = 9;

  // CSV output configuration
  const std::string output_fn = "collision_output.csv";
  std::ofstream csv_file;

  // Random number generator
  std::mt19937_64 rng(0); // Seed with 0 for reproducibility
  std::exponential_distribution<double> exp_dist(1.0 / x0);

  std::cout << "\n=== SUPERDROPLET COLLISION-COALESCENCE MODEL ===\n\n";

  // SIMULATION SETUP
  physics::ModelConfig model_config{
      t_c, delta_v, static_cast<unsigned int>(n_part), kernel, DEBUG};

  std::cout << "MODEL SETUP\n";
  std::cout << "  step_seconds: " << model_config.step_seconds << "\n";
  std::cout << "  delta_v: " << model_config.delta_v << "\n";
  std::cout << "  num_droplets: " << model_config.num_droplets << "\n";
  std::cout << "  kernel: " << physics::kernel_name(model_config.kernel)
            << "\n";
  std::cout << "  debug: " << (model_config.debug ? "true" : "false") << "\n";

  // Initialize the droplet array
  const double total_droplets = delta_v * n0;
  const size_t xi_i = static_cast<size_t>(
      std::floor(total_droplets / static_cast<double>(n_part)));

  // Compute droplet masses according to an exponential distribution
  std::vector<double> x_grid;
  x_grid.reserve(n_part);
  for (size_t i = 0; i < n_part; ++i) {
    x_grid.push_back(exp_dist(rng));
  }

  // Sort the masses
  std::sort(x_grid.begin(), x_grid.end());

  // Compute the radius grid
  std::vector<double> r_grid;
  r_grid.reserve(n_part);
  for (double x : x_grid) {
    r_grid.push_back(std::cbrt(x * constants::THREE_FOURTH / constants::PI));
  }

  // Create the initial droplet vector
  std::vector<Droplet> droplets;
  droplets.reserve(n_part);
  for (double radius : r_grid) {
    droplets.emplace_back(xi_i, radius);
  }

  std::cout << "\nGRID SETUP\n";
  std::cout << "   radii: " << std::scientific << std::setprecision(3)
            << droplets[0].radius() << " - " << droplets[n_part - 1].radius()
            << " m\n";
  std::cout << "  volume: " << droplets[0].volume() << " - "
            << droplets[n_part - 1].volume() << " m^3\n";
  std::cout << "    mass: " << droplets[0].mass() << " - "
            << droplets[n_part - 1].mass() << " kg\n";

  std::cout << "\nSD SETUP\n";
  std::cout << "        N_s: " << n_part << "\n";
  std::cout << "       xi_i: " << xi_i << "\n";
  std::cout << "N per SD_xi: "
            << static_cast<size_t>(total_droplets) / xi_i / n_part << "\n";

  const double wm0 = total_water(droplets);
  std::cout << "\nInitial total water mass = " << std::scientific
            << std::setprecision(3) << wm0 << " kg\n\n";

  // Initialize CSV output
  if (PLOT) {
    csv_file.open(output_fn);
    if (!csv_file.is_open()) {
      std::cerr << "Error: Could not open output file " << output_fn << "\n";
      return 1;
    }

    auto bin_grid = io::bin_droplets(droplets);
    write_csv_header(csv_file, bin_grid);
  }

  // Main loop
  std::vector<double> masses = {wm0};
  utils::Stopwatch stopwatch(0);
  unsigned int step = 0;

  std::cout << "BEGINNING MAIN SIMULATION LOOP\n\n";

  auto start_time = std::chrono::steady_clock::now();

  while (stopwatch.total_seconds() <= t_end) {
    if (PLOT && step % plot_dt == 0) {
      std::cout << "Plotting. Stopwatch = " << stopwatch.minutes << " min "
                << stopwatch.seconds << " s\n";

      auto bin_grid = io::bin_droplets(droplets);

      // Compute rolling median of the values of the bin grid
      auto smoothed_values = smooth_values(bin_grid.values, smooth_window);

      write_csv_row(csv_file, stopwatch.total_seconds(), smoothed_values);
    }

    auto result = physics::collision_step(droplets, model_config, rng);

    // Only compute total water when needed for plotting
    if (PLOT && step % plot_dt == 0) {
      masses.push_back(result.total_water);
    }

    stopwatch.increment(t_c);

    // Print status update at regular intervals
    if (DEBUG && (step % STATUS_DT == 0)) {
      std::cout << "STEP: " << step << " (" << stopwatch.minutes << " min "
                << stopwatch.seconds << " s)"
                << " | Collisions: " << result.counter
                << " | Prob: " << std::fixed << std::setprecision(2)
                << result.min_prob << " - " << result.max_prob << " ["
                << result.big_probs << "]"
                << " | Water: " << std::scientific << std::setprecision(3)
                << result.total_water << " kg\n";
    }

    step++;
  }

  std::cout << "\n";

  auto end_time = std::chrono::steady_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
      end_time - start_time);

  if (PLOT) {
    csv_file.close();
  }

  std::cout << "Simulation completed successfully.\n";
  double elapsed_seconds = duration.count() / 1000.0;
  std::cout << "Runtime: " << std::fixed << std::setprecision(2)
            << elapsed_seconds << " seconds\n";

  const double final_mass = masses.back();
  std::cout << "Remaining water mass: " << std::scientific
            << std::setprecision(3) << final_mass << " kg (" << std::fixed
            << std::setprecision(1) << (final_mass / wm0 * 100.0)
            << "% of initial)\n";

  // Write timing to file
  std::ofstream time_file("time.out");
  time_file << std::fixed << std::setprecision(6) << elapsed_seconds << "\n";
  time_file.close();

  return 0;
}
