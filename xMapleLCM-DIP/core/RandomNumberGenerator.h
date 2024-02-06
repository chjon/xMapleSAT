/*************************************************************************[RandomNumberGenerator.h]
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

#ifndef Minisat_RandomNumberGenerator_h
#define Minisat_RandomNumberGenerator_h

namespace Minisat {
    /**
     * @brief This class generates random numbers.
     * 
     */
    class RandomNumberGenerator {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // MEMBER VARIABLES

        /// @brief The seed for generating pseudorandom numbers
        double random_seed;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // CONSTRUCTORS

        /**
         * @brief Construct a new RandomNumberGenerator object
         * 
         */
        RandomNumberGenerator();

        /**
         * @brief Destroy the RandomNumberGenerator object
         * 
         */
        ~RandomNumberGenerator() = default;

        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATIC HELPER FUNCTIONS

        /**
         * @brief Get a random float 0 <= x < 1
         * 
         * @param seed The seed for generating the next pseudorandom number
         * @return A pseudorandom float 0 <= x < 1
         * 
         * @pre seed is not 0
         */
        static double drand(double& seed);

        /**
         * @brief Get a random integer 0 <= x < size
         * 
         * @param seed The seed for generating the next pseudorandom number
         * @param size The limit for the generated value
         * @return a pseudorandom integer 0 <= x < size 
         */
        static int irand(double& seed, int size);

        ///////////////////////////////////////////////////////////////////////////////////////////
        // PUBLIC API

        /**
         * @brief Get a random float 0 <= x < 1
         * 
         * @return A pseudorandom float 0 <= x < 1 
         */
        double drand();

        /**
         * @brief Get a random integer 0 <= x < size
         * 
         * @param size The limit for the generated value
         * @return a pseudorandom integer 0 <= x < size 
         */
        int irand(int size);
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    ////////////////////////////
    // STATIC HELPER FUNCTIONS

    inline double RandomNumberGenerator::drand(double& seed) {
        seed *= 1389796;
        int q = (int)(seed / 2147483647);
        seed -= (double)q * 2147483647;
        return seed / 2147483647;
    }

    inline int RandomNumberGenerator::irand(double& seed, int size) {
        return (int)(drand(seed) * size);
    }

    ///////////////
    // PUBLIC API

    inline double RandomNumberGenerator::drand() {
        return drand(random_seed);
    }

    inline int RandomNumberGenerator::irand(int size) {
        return irand(random_seed, size);
    }
}

#endif