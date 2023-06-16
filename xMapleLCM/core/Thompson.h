#ifndef _thompson_hpp_INCLUDED
#define _thompson_hpp_INCLUDED

#include <algorithm>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <vector>
#include <boost/random.hpp>
#include <boost/random/discrete_distribution.hpp>
#include <boost/random/mersenne_twister.hpp>

using beta_distribution = boost::random::beta_distribution<>;

// Thompson Variables
//
typedef struct Thompson_variable {
    // gen is a Mersenne Twister random generator. We initialize it here to keep
    // the binary deterministic.
    boost::mt19937 gen;

    const size_t num_arms;
    std::vector<double> alphas;                 // Alpha values per action
    std::vector<double> betas;                  // Beta values per action
    std::vector<beta_distribution> prior_dists; // Beta distributions for each action

    Thompson_variable (size_t num_arms) : num_arms (num_arms) {
        alphas = std::vector<double> (num_arms, 0);
        betas  = std::vector<double> (num_arms, 0);

        // Initialize the prior distributions with alpha=1 beta=1
        for (size_t i = 0; i < num_arms; i++) prior_dists.push_back (beta_distribution (1, 1));
    }

    inline size_t select_lever () {
        // Sample a random value from each prior distribution and keep track
        // of the arm with the highest sampled value.
        size_t max_i = 0;
        double max = 0;
        for (size_t i = 0; i < prior_dists.size (); i++) {
            const double sample = prior_dists[i] (gen);
            if (sample > max) {
                max_i = i;
                max = sample;
            }
        }

        // Return the arm with the highest sampled value
        return max_i;
    }

    inline void update_dist(unsigned int action, bool isWin, double decayFactor) {
        // Decay wins/losses
        for (size_t i = 0; i < num_arms; i++) {
            alphas[i] *= decayFactor;
            betas [i] *= decayFactor;
        }
        
        // Bump wins/losses
        (isWin ? alphas : betas)[action]++;

        // Update the prior distribution of the chosen bandit
        // add 1 because the beta distribution behaves weirdly for parameter values less than 1
        prior_dists[action] = beta_distribution(1 + alphas[action], 1 + betas [action]);
    }

} Thompson_var;

#endif