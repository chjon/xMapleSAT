/**********************************************************************************[BCPRLManager.h]
MiniSat -- Copyright (c) 2003-2006, Niklas Een, Niklas Sorensson
           Copyright (c) 2007-2010, Niklas Sorensson

MapleSAT_Refactor, based on MapleSAT -- Copyright (c) 2022, Jonathan Chung, Vijay Ganesh, Sam Buss

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction,
including without limitation the rights to use, copy, modify, merge, publish, distribute,
sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or
substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
**************************************************************************************************/

#ifndef Minisat_BCPRLManager_h
#define Minisat_BCPRLManager_h

#include "core/RandomNumberGenerator.h"
#include "core/SolverTypes.h"
#include "core/Thompson.h"

// Define BCP switching mode
#ifndef ENABLE_PRIORITY_BCP_RL
    #define ENABLE_PRIORITY_BCP_RL false
#endif
#ifndef ENABLE_PRIORITY_BCP_RANDOM
    #define ENABLE_PRIORITY_BCP_RANDOM false
#endif

namespace Minisat {
    // Forward declarations
    class Solver;

    class BCPRLManager {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief The decay to apply to the beta distributions' shape parameters
        const double BETA_DECAY;

        /// @brief The decay to apply to the historical score
        const double SCORE_DECAY;

        /// @brief The historical (EMA) performance of the RL agent
        double historicalScore;

        Thompson_var thompson;

        // RL Score Counters

        uint64_t lbdsum;
        uint64_t prevConflicts;
        uint64_t prevDecisions;
        uint64_t prevPropagations;

    public:
        /// @brief The current variant of BCP to use
        BCPMode current_bcpmode;

        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATISTICS

        uint64_t num_delayed;

    protected:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // SOLVER REFERENCES

        RandomNumberGenerator& randomNumberGenerator;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new BCPRLManager object
         * 
         * @param s Reference to main solver object
         */
        BCPRLManager(Solver& s);

        /**
         * @brief Destroy the BCPRLManager object
         * 
         */
        ~BCPRLManager() = default;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATE MODIFICATION

        template <RLScoreType scoretype>
        inline double get_prev_round_score(const BCPRLStats& stats);

        inline void clearScores(const BCPRLStats& stats);

        BCPMode selectNextMode(BCPMode prev_bcpmode, const BCPRLStats& stats);

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // EVENT HANDLERS

        inline void handleEventRestarted(const BCPRLStats& stats);

        inline void handleEventLearntClause(uint64_t lbd);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    template <> inline double BCPRLManager::get_prev_round_score<RLScoreType::GLR>(const BCPRLStats& stats) {
        return (stats.conflicts - prevConflicts) / static_cast<double>(stats.decisions - prevDecisions);
    }

    template <> inline double BCPRLManager::get_prev_round_score<RLScoreType::LBD>(const BCPRLStats& stats) {
        return static_cast<double>(stats.conflicts - prevConflicts) / lbdsum; // 1 over avgLBD
    }

    template <> inline double BCPRLManager::get_prev_round_score<RLScoreType::PPD>(const BCPRLStats& stats) {
        return (stats.propagations - prevPropagations) / static_cast<double>(stats.decisions - prevDecisions);
    }

    inline void BCPRLManager::clearScores(const BCPRLStats& stats) {
        prevConflicts    = stats.conflicts;
        prevPropagations = stats.propagations;
        prevDecisions    = stats.decisions;
        lbdsum           = 0;
    }

    inline BCPMode BCPRLManager::selectNextMode(BCPMode prev_bcpmode, const BCPRLStats& stats) {
        if (stats.restarts > 1) {
            // Update (bump and decay) reward values
            const double prevRoundScore = get_prev_round_score<RLScoreType::LBD>(stats);
            thompson.update_dist(static_cast<size_t>(prev_bcpmode), prevRoundScore >= historicalScore, BETA_DECAY);

            // Update historical score as a weighted average
            historicalScore = historicalScore * SCORE_DECAY + prevRoundScore * (1 - SCORE_DECAY);
        }

        // Pick new BCP mode
        return static_cast<BCPMode>(thompson.select_lever ());
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // EVENT HANDLERS

    inline void BCPRLManager::handleEventRestarted(const BCPRLStats& stats) {
        if (ENABLE_PRIORITY_BCP_RL) {
            current_bcpmode = selectNextMode(current_bcpmode, stats);

            if (current_bcpmode == BCPMode::DELAYED) num_delayed++;
        } else if (ENABLE_PRIORITY_BCP_RANDOM) {
            current_bcpmode = (randomNumberGenerator.drand() > 0.5) ? BCPMode::DELAYED : BCPMode::IMMEDIATE;
        }

        #ifdef ENABLE_PRIORITY_BCP
        clearScores(stats);
        #endif
    }

    inline void BCPRLManager::handleEventLearntClause(uint64_t lbd) {
        lbdsum += lbd;
    }
}

#endif