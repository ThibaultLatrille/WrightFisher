#include "wright_fisher.hpp"
#include <fstream>
#include "io.hpp"
#include "tclap/CmdLine.h"

using namespace TCLAP;
using namespace std;


class ArgParse {
  protected:
    TCLAP::CmdLine& cmd;

  public:
    explicit ArgParse(TCLAP::CmdLine& cmd) : cmd{cmd} {}

    TCLAP::ValueArg<std::string> output_path{"o", "output", "output path", true, "", "string", cmd};
    TCLAP::ValueArg<u_long> seed{
        "", "seed", "Random number generation seed", false, 0, "u_long", cmd};
    TCLAP::ValueArg<double> mutation_rate_per_generation{"", "mutation_rate_per_generation",
        "Mutation rate (at the root)", false, 1e-6, "double", cmd};
    TCLAP::ValueArg<u_long> population_size{
        "", "population_size", "Effective population size", false, 100, "u_long", cmd};
    TCLAP::ValueArg<u_long> number_of_generations{"", "number_of_generations",
        "Number of generations for the simulation", false, 1000, "u_long", cmd};

    TCLAP::SwitchArg gamma_reflected{"", "gamma_reflected", "True if reflected gamma.", cmd, false};
    TCLAP::ValueArg<double> gamma_mean{
        "", "gamma_mean", "Mean value of selection coefficient effect", false, 1.0, "double", cmd};
    TCLAP::ValueArg<double> gamma_shape{
        "", "gamma_shape", "Shape of the gamma distribution", false, 1.0, "double", cmd};
    TCLAP::ValueArg<u_long> sequence_size{
        "", "sequence_size", "Sequence size", false, 30, "u_long", cmd};
};


int main(int argc, char* argv[]) {
    CmdLine cmd{"PolyDfe", ' ', "0.1"};
    ArgParse args(cmd);
    cmd.parse(argc, argv);

    u_long arg_seed = args.seed.getValue();
    cout << "Random generator seed: " << arg_seed << endl;
    generator.seed(arg_seed);

    auto dfe_param = DfeParameters(args.sequence_size.getValue(), args.gamma_reflected.getValue(),
        args.gamma_mean.getValue(), args.gamma_shape.getValue());
    auto population = Population(
        dfe_param, args.population_size.getValue(), args.mutation_rate_per_generation.getValue());
    for (u_long i = 0; i < args.population_size.getValue(); ++i) { population.forward(); }
    for (u_long i = 0; i < args.number_of_generations.getValue(); ++i) { population.forward(); }

    return 0;
}