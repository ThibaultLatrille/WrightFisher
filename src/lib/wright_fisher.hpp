#include <set>
#include "codon.hpp"
#include "io.hpp"
#include "random.hpp"
#include "statistic.hpp"

bool is_synonymous(char codon_from, char codon_to) {
    return codonLexico.codon_to_aa.at(codon_from) == codonLexico.codon_to_aa.at(codon_to);
}

class Substitution {
  public:
    Substitution(char codon_from, char codon_to, double time_event, double time_between)
        : codon_from{codon_from},
          codon_to{codon_to},
          time_event{time_event},
          time_between{time_between} {}

    bool is_synonymous() const { return ::is_synonymous(codon_from, codon_to); }

    bool is_non_synonymous() const { return !::is_synonymous(codon_from, codon_to); }

    void add_to_trace(Trace &trace) const {
        trace.add("codon_from", codonLexico.codon_string(codon_from));
        trace.add("codon_to", codonLexico.codon_string(codon_to));
        trace.add("aa_from", codonLexico.codon_aa_string(codon_from));
        trace.add("aa_to", codonLexico.codon_aa_string(codon_to));
        trace.add("synonymous", is_synonymous());
        trace.add("time_event", time_event);
        trace.add("time_between", time_between);
    }

    char codon_from;
    char codon_to;
    double time_event;
    double time_between;
};

class DfeParameters {
  public:
    u_long exon_size;
    bool reflected;
    double mean;
    double shape;
    double scale;

    ~DfeParameters() = default;
    explicit DfeParameters() = default;

    explicit DfeParameters(u_long exon_size, bool reflected, double mean, double shape)
        : exon_size{exon_size}, reflected{reflected}, mean{mean}, shape{shape} {
        if (shape != 0.0) { scale = mean / shape; }
    }

    double selection_coefficient(u_long site, char aa_to, double const &pop_size) const {
        auto bk_seed = generator();
        generator.seed(site * 20 + aa_to);
        double s{-mean};
        if (shape != 0.0) { s = -std::gamma_distribution<double>(shape, scale)(generator); }
        if (reflected) {
            // Estimating the distribution of fitness effects from DNA sequence data: Implications
            // for the molecular clock. Gwenaël Piganeau and Adam Eyre-Walker, PNAS, 2003.
            // https://doi.org/10.1073/pnas.1833064100
            double S = -4 * pop_size * s;
            double x = exp(S) / (exp(S) + 1);
            assert(0.5 <= x && x <= 1.0);
            double draw = std::uniform_real_distribution<double>(0.0, 1.0)(generator);
            if (draw > x) {
                s = -s;
                assert(s >= 0);
            }
        }
        generator.seed(bk_seed);
        return s;
    }

    u_long nbr_sites() const { return exon_size; }

    std::array<double, 20> aa_selection_coefficients(
        std::vector<char> const &codon_seq, u_long site, double const &pop_size) const {
        std::array<double, 20> profiles{};
        for (char aa = 0; aa < 20; ++aa) {
            profiles[aa] = selection_coefficient(site, aa, pop_size);
        }
        return profiles;
    }

    // Method computing the equilibrium frequencies for one site.
    static std::array<double, 64> codon_frequencies(
        std::array<double, 20> const &aa_selection_coefficient, double const &nuc_matrix,
        double const &pop_size) {
        std::array<double, 64> codon_frequencies{0};
        // For each site of the vector of the site frequencies.
        for (char codon{0}; codon < 64; codon++) {
            double codon_freq = 1.0;

            if (codonLexico.codon_to_aa[codon] != 20) {
                codon_frequencies[codon] =
                    codon_freq *
                    exp(4 * aa_selection_coefficient[codonLexico.codon_to_aa[codon]] * pop_size);
            } else {
                codon_frequencies[codon] = 0.;
            }
        }

        double sum_freq = std::accumulate(codon_frequencies.begin(), codon_frequencies.end(), 0.0);
        for (char codon{0}; codon < 64; codon++) { codon_frequencies[codon] /= sum_freq; }

        return codon_frequencies;
    }
};

class Haplotype {
  public:
    // The nbr_copies of sites in the sequence (each position is a codon, thus the DNA sequence is 3
    // times greater).
    u_long nbr_copies{0};
    double log_fitness{0.0};
    std::unordered_map<u_long, char> diff_sites;

    ~Haplotype() = default;
    Haplotype() = default;

    // Constructor of Reference_seq.
    // size: the size of the DNA sequence.
    Haplotype(u_long const nbr_copies, double const DfeParameters)
        : nbr_copies{nbr_copies}, log_fitness{DfeParameters} {};

    bool check_consistency() const {
        for (auto &diff : diff_sites) {
            if (codonLexico.codon_to_aa.at(diff.second) == 20) {
                std::cerr << "The haplotype contains a stop codon" << std::endl;
                return false;
            }
        }
        return true;
    }

    struct GreaterThan {
        bool operator()(Haplotype const &left, Haplotype const &right) {
            return left.nbr_copies > right.nbr_copies;
        }

        bool operator()(Haplotype const &left, u_long right) { return left.nbr_copies > right; }

        bool operator()(u_long left, Haplotype const &right) { return left > right.nbr_copies; }
    };
};

// Class representing a genetically linked sequences (exons are unlinked between them)
class Population {
  public:
    // Reference sequence
    u_long nbr_sites{};
    u_long nbr_nucleotides{};
    std::vector<char> codon_seq;

    double mutation_rate{};
    u_long population_size{};
    DfeParameters fitness_map;

    double time_current{0};
    // Haplotypes
    std::vector<Haplotype> haplotype_vector;

    // Keep track of mutations and fixations
    std::vector<Substitution> substitutions{};

    mutable std::binomial_distribution<u_long> binomial_distrib;

    ~Population() = default;
    explicit Population() = default;

    // Constructor
    Population(
        DfeParameters const &f_state, u_long const &population_size, double const &_mutation_rate)
        : nbr_sites{u_long(f_state.nbr_sites())},
          nbr_nucleotides{u_long(3 * f_state.nbr_sites())},
          codon_seq(f_state.nbr_sites(), 0),
          mutation_rate{_mutation_rate},
          population_size{population_size},
          fitness_map{f_state},
          haplotype_vector{} {
        for (u_long site{0}; site < nbr_sites; site++) {
            // Draw codon from codon frequencies
            std::array<double, 64> codon_freqs = fitness_map.codon_frequencies(
                fitness_map.aa_selection_coefficients(codon_seq, site, population_size),
                mutation_rate, population_size);
            std::discrete_distribution<char> freq_codon_distr(
                codon_freqs.begin(), codon_freqs.end());
            codon_seq[site] = freq_codon_distr(generator);
        }
        binomial_distrib = std::binomial_distribution<u_long>(
            2 * population_size * nbr_nucleotides, mutation_rate);

        haplotype_vector.emplace_back(2 * population_size, 0.0);
        assert(check_consistency());
    }

    bool check_consistency() const {
        if (haplotype_vector.size() > 2 * population_size) {
            std::cerr << "Too many haplotypes" << std::endl;
            return false;
        }
        if (haplotype_vector.empty()) {
            std::cerr << "No haplotype" << std::endl;
            return false;
        }

        u_long nbr_copies{0};
        for (auto &haplotype : haplotype_vector) {
            if (!haplotype.check_consistency()) {
                std::cerr << "The haplotype is not consistent" << std::endl;
                return false;
            }
            nbr_copies += haplotype.nbr_copies;
        }
        if (nbr_copies != 2 * population_size) {
            std::cerr << "The number of copies is not equal to the population size." << std::endl;
            return false;
        }
        return true;
    }

    void forward() {
        mutation();
        selection_and_drift();
        extinction();
        fixation();
    }

    u_long random_nuc_site() const {
        return std::uniform_int_distribution<u_long>(0, nbr_nucleotides - 1)(generator);
    }

    void mutation() {
        // Randomly choose the number of mutations at this generation (for this exon)
        u_long binomial_draw = binomial_distrib(generator);

        // Early break if 0 mutation are drawn
        if (binomial_draw == 0) { return; }

        // Compute the distribution of haplotype frequency in the population
        std::discrete_distribution<u_long> haplotype_distri;
        if (haplotype_vector.size() != 1) {
            std::vector<u_long> nbr_copies(haplotype_vector.size(), 0);
            std::transform(haplotype_vector.begin(), haplotype_vector.end(), nbr_copies.begin(),
                [](Haplotype const &h) { return h.nbr_copies; });
            haplotype_distri =
                std::discrete_distribution<u_long>(nbr_copies.begin(), nbr_copies.end());
        }

        for (u_long nbr_draws{0}; nbr_draws < binomial_draw; nbr_draws++) {
            // Randomly choose the haplotype
            u_long hap_id{0};
            if (haplotype_vector.size() != 1) { hap_id = haplotype_distri(generator); }

            // Early continue (next iteration) if the selected haplotype is not held by any
            // individual
            if (haplotype_vector.at(hap_id).nbr_copies == 0) { continue; }

            // Randomly choose the site
            u_long nuc_site = random_nuc_site();
            u_long codon_site = nuc_site / 3;
            auto nuc_pos = static_cast<char>(nuc_site % 3);

            // The corresponding codon given the selected haplotype and site
            char codon_from{0};
            if (haplotype_vector.at(hap_id).diff_sites.count(codon_site) > 0) {
                codon_from = haplotype_vector.at(hap_id).diff_sites.at(codon_site);
            } else {
                codon_from = codon_seq.at(codon_site);
            }

            // The selected nucleotide
            std::array<char, 3> triplet_nuc = codonLexico.codon_to_triplet.at(codon_from);

            // Randomly choose the target nucleotide
            triplet_nuc[nuc_pos] = std::uniform_int_distribution<char>(0, 3)(generator);
            char codon_to = codonLexico.triplet_to_codon(
                triplet_nuc.at(0), triplet_nuc.at(1), triplet_nuc.at(2));
            if (codonLexico.codon_to_aa.at(codon_to) == 20) { continue; }

            Haplotype haplotype = haplotype_vector.at(hap_id);

            if (!is_synonymous(codon_from, codon_to)) {
                // Update the fitness of the new haplotype
                auto bk = codon_seq.at(codon_site);
                codon_seq[codon_site] = codon_from;
                haplotype.log_fitness += fitness_map.selection_coefficient(
                    codon_site, codonLexico.codon_to_aa[codon_to], population_size);
                codon_seq[codon_site] = bk;
            }

            // Depending on whether the mutation is back to the reference sequence
            if (codon_to == codon_seq.at(codon_site)) {
                assert(haplotype.diff_sites.count(codon_site) > 0);
                haplotype.diff_sites.erase(codon_site);
            } else {
                haplotype.diff_sites[codon_site] = codon_to;
            }

            // Update the vector of haplotypes
            haplotype_vector.at(hap_id).nbr_copies--;
            haplotype.nbr_copies = 1;
            haplotype_vector.push_back(haplotype);

        }
    }

    void selection_and_drift() {
        u_long children_tot = 2 * population_size;

        // Early break if only 1 haplotype
        if (haplotype_vector.size() == 1) {
            haplotype_vector.front().nbr_copies = children_tot;
            return;
        }

        // The fitness associated to each haplotype (weigthed by the number of copies)
        std::vector<double> fitnesses(haplotype_vector.size(), 0);
        std::transform(haplotype_vector.begin(), haplotype_vector.end(), fitnesses.begin(),
            [](Haplotype const &h) { return (1.0 + h.log_fitness) * h.nbr_copies; });
        double fit_tot = accumulate(fitnesses.begin(), fitnesses.end(), 0.0);

        // Random draws from the multinomial distribution
        for (std::size_t i_hap{0}; i_hap < haplotype_vector.size() - 1; i_hap++) {
            haplotype_vector.at(i_hap).nbr_copies = std::binomial_distribution<u_long>(
                children_tot, fitnesses.at(i_hap) / fit_tot)(generator);
            children_tot -= haplotype_vector.at(i_hap).nbr_copies;
            fit_tot -= fitnesses.at(i_hap);
        }
        assert(children_tot <= 2 * population_size);
        haplotype_vector.back().nbr_copies = children_tot;
    }

    void extinction() {
        // Early break if only 1 haplotype
        if (haplotype_vector.size() == 1) { return; }

        // Sort the vector of haplotypes by number of copies
        if (!is_sorted(
                haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan{})) {
            sort(haplotype_vector.begin(), haplotype_vector.end(), Haplotype::GreaterThan{});
        }

        // Remove haplotypes with 0 copies
        if (haplotype_vector.back().nbr_copies == 0) {
            auto low_bound = lower_bound(
                haplotype_vector.begin(), haplotype_vector.end(), 0, Haplotype::GreaterThan{});
            if (low_bound != haplotype_vector.end()) {
                haplotype_vector.erase(low_bound, haplotype_vector.end());
            }
        }
    }

    void fixation() {
        if (haplotype_vector.size() == 1) { return; }

        // ! Copy before the loop is necessary, since it can be updated during the loop
        auto diffs = haplotype_vector.begin()->diff_sites;

        // For each site (different to the reference) of the most common haplotype
        for (auto const &diff : diffs) {
            u_long site = diff.first;
            char codon_to = diff.second;
            bool polymorphic_site = false;
            std::set<char> poly_codons{codon_to};
            for (std::size_t hap_id{1}; hap_id < haplotype_vector.size(); hap_id++) {
                if (haplotype_vector.at(hap_id).diff_sites.count(site) == 0) {
                    // In the case the site of this haplotype is the same than the reference
                    polymorphic_site = true;
                    break;
                } else {
                    // In the case the site of this haplotype is different from the reference
                    poly_codons.insert(haplotype_vector.at(hap_id).diff_sites.at(site));
                }
            }

            // Early continue (next iteration) if the site is polymorphic
            if (polymorphic_site) { continue; }

            char codon_from = codon_seq.at(site);
            if (poly_codons.size() > 1) {
                // In the case the reference codon is not present in the alternative haplotypes,
                // but they are several alternative codons, we have to find the most common.
                std::unordered_map<char, u_long> codon_copies{};
                for (char codon : poly_codons) { codon_copies[codon] = 0; }
                for (auto const &haplotype : haplotype_vector) {
                    codon_copies.at(haplotype.diff_sites.at(site)) += haplotype.nbr_copies;
                }
                codon_to = std::max_element(
                    codon_copies.begin(), codon_copies.end(), [](auto const &p1, auto const &p2) {
                        return p1.second < p2.second;
                    })->first;
            }
            if (!is_synonymous(codon_from, codon_to)) {
                double s = fitness_map.selection_coefficient(
                    site, codonLexico.codon_to_aa[codon_to], population_size);
                for (auto &haplotype : haplotype_vector) { haplotype.log_fitness -= s; }
            }
            codon_seq[site] = codon_to;
            // Update the vector of haplotypes
            for (auto &haplotype : haplotype_vector) {
                if (haplotype.diff_sites.at(site) == codon_to) {
                    haplotype.diff_sites.erase(site);
                } else {
                    assert(poly_codons.size() > 1);
                    assert(poly_codons.find(haplotype.diff_sites.at(site)) != poly_codons.end());
                }
            }

            double t_between = time_current;
            if (!substitutions.empty()) { t_between -= substitutions.back().time_event; }
            auto sub = Substitution(codon_from, codon_to, time_current, t_between);
            substitutions.push_back(sub);
        }
    }
};