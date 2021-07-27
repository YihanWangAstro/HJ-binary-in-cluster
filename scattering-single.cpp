#include "/home/yihanw/repositories/SpaceHub/src/spaceHub.hpp"           //change to your path to spacehub
#include "/home/yihanw/repositories/SpaceHub/src/taskflow/taskflow.hpp"  //change to your path to spacehub
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

void job(size_t thread_id, size_t scattering_num, double sigma) {
    double v_inf = sigma * 3.66_kms;
    double a_j1 = 5_AU;
    double a_j2 = 15_AU;
    double r_start = 50 * a_j2;

    std::fstream post_flyby_file("post-flyby-single-" + std::to_string(thread_id) + ".txt", std::ios::out);

    print(post_flyby_file,
          "m0+mj1+mj2 [solar mass], m1 [solar mass], p_hyper [au], e_hyper, inc_hyper [arc], Omega_hyper [arc], "
          "omega_hyper [arc], nv_hyper [arc], "
          "m0 [solar mass], mj1 [solar mass], a_j1 [au], e_j1, inc_j1 [arc], Omega_j1 [arc], omega_j1 [arc], "
          "nv_j1 [arc], m0+mj1 [solar mass], mj2 [solar mass], a_j2 [au], e_j2, inc_j2 [arc], Omega_j2 [arc], "
          "omega_j2 [arc], nv_j2 [arc]\n post-flyby particles data [time, id, mass, x, y, z, vx, vy, "
          "vz](per line)\n");

    print(post_flyby_file, std::setprecision(16));

    print(std::cout, "starting job on thread ", thread_id, " \n");

    for (size_t i = 0; i < scattering_num; ++i) {
        Particle star{1_Ms};
        Particle jupiter1{1_Mj};
        Particle jupiter2{1_Mj};
        Particle star1{2_Ms};

        // create planetary system with two giant planets
        double planet_inc = random::Uniform(0, consts::pi);

        auto jupiter1_orb = Elliptic(star.mass, jupiter1.mass, a_j1, 0.0, planet_inc, isotherm, isotherm, isotherm);

        auto jupiter2_orb =
            Elliptic(star.mass + jupiter1.mass, jupiter2.mass, a_j2, 0.0, planet_inc, isotherm, isotherm, isotherm);

        move_particles(jupiter1_orb, jupiter1);

        move_to_COM_frame(star, jupiter1);

        move_particles(jupiter2_orb, jupiter2);

        // create scattering hyperbolic orbit

        double b_max = calc_max_impact_parameter(a_j2 * 3, v_inf, M_tot(star, jupiter1, jupiter2, star1));

        auto incident_orb =
            scattering::incident_orbit(M_tot(star, jupiter1, jupiter2), star1.mass, v_inf, b_max, r_start);

        move_particles(incident_orb, star1);

        move_to_COM_frame(star, jupiter1, jupiter2, star1);

        double scattering_t_end = 2 * time_to_periapsis(group(star, jupiter1, jupiter2), star1);

        Solver sim{0, star, jupiter1, jupiter2, star1};

        Solver::RunArgs args;

        args.add_stop_condition(scattering_t_end);

        args.add_stop_point_operation([&](auto& ptc, auto h) {
            print(post_flyby_file, incident_orb, ',', jupiter1_orb, ',', jupiter2_orb, '\n', ptc, '\n');
        });

        sim.run(args);
    }
}

int main(int argc, char** argv) {
    size_t n = 5000;  // total scattering number
    size_t job_num = 12;

    tf::Executor executor;

    double sigma = std::atof(argv[1]);

    for (size_t i = 0; i < job_num; ++i) {
        executor.silent_async(job, i, n / job_num, sigma);
    }
    executor.wait_for_all();  // wait all jobs to be finished
    return 0;
}
