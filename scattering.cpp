#include "/home/yihanw/repositories/SpaceHub/src/spaceHub.hpp"
#include "/home/yihanw/repositories/SpaceHub/src/taskflow/taskflow.hpp"
using namespace hub;
using namespace unit;
using namespace callback;
using namespace force;
using namespace orbit;
using f = Interactions<NewtonianGrav>;
using Solver = methods::DefaultMethod<f>;
using Particle = Solver::Particle;

// Q_max: max closest approach
double calc_max_impact_parameter(double Q_max, double v_inf, double M_tot) {
    return Q_max * sqrt(1 + 2 * M_tot / v_inf / v_inf / Q_max);
}

void job(size_t thread_id, size_t scattering_num, double b_semi_major_axis) {
    double v_inf = 3.66_kms;
    double r_start = 50 * b_semi_major_axis;
    double a_j1 = 5_AU;
    double a_j2 = 15_AU;

    std::fstream post_flyby_file("post-flyby-" + std::to_string(thread_id) + ".txt", std::ios::out);

    print(post_flyby_file,
          "m0 [solar mass], mj1 [solar mass], a_j1 [au], e_j1, inc_j1 [arc], Omega_j1 [arc], omega_j1 [arc], "
          "nv_j1 [arc], m0+mj1 [solar mass], mj2 [solar mass], a_j2 [au], e_j2, inc_j2 [arc], Omega_j2 [arc], "
          "omega_j2 [arc], nv_j2 [arc], m1 [solar mass], m2 [solar mass], a_b [au], e_b, inc_b [arc], Omega_b "
          "[arc], omega_b [arc], nv_b [arc]\n post-flyby particles data [time, id, mass, x, y, z, vx, vy, "
          "vz](per line)\n");

    print(post_flyby_file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle star{1_Ms};
        Particle jupiter1{1_Mj};
        Particle jupiter2{1_Mj};
        Particle star1{1_Ms};
        Particle star2{1_Ms};

        // create planetary system with two giant planets
        double planet_inc = random::Uniform(0, consts::pi);

        auto jupiter1_orb = Elliptic(star.mass, jupiter1.mass, a_j1, 0.0, planet_inc, isotherm, isotherm, isotherm);

        auto jupiter2_orb =
            Elliptic(star.mass + jupiter1.mass, jupiter2.mass, a_j2, 0.0, planet_inc, isotherm, isotherm, isotherm);

        move_particles(jupiter1_orb, jupiter1);

        move_to_COM_frame(star, jupiter1);

        move_particles(jupiter2_orb, jupiter2);

        // create binary star
        auto bin_orb = Elliptic(star1.mass, star2.mass, b_semi_major_axis, 0.0, isotherm, isotherm, isotherm, isotherm);

        move_particles(bin_orb, star2);

        move_to_COM_frame(star1, star2);

        // create scattering hyperbolic orbit

        double b_max = calc_max_impact_parameter(a_j2 * 2, v_inf, M_tot(star, jupiter1, jupiter2, star1, star2)) +
                       b_semi_major_axis / 2;

        auto incident_orb =
            scattering::incident_orbit(M_tot(star, jupiter1, jupiter2), M_tot(star1, star2), v_inf, b_max, r_start);

        move_particles(incident_orb, star1, star2);

        move_to_COM_frame(star, jupiter1, jupiter2, star1, star2);

        double scattering_t_end = 2 * time_to_periapsis(group(star, jupiter1, jupiter2), group(star1, star2));

        Solver sim{0, star, jupiter1, jupiter2, star1, star2};

        Solver::RunArgs args;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            print(post_flyby_file, jupiter1_orb, ',', jupiter2_orb, ',', bin_orb, '\n', ptc, '\n');
        });

        sim.run(args);
    }
}

int main() {
    size_t n = 5000;  // total scattering number
    size_t job_num = 12;

    tf::Executor executor;

    double a_binary = 30_AU;

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, a_binary);
    }
    executor.wait_for_all();  // wait all jobs to be finished
    return 0;
}
