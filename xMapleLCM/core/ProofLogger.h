/********************************************************************************[ProofLogger.h]
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

#ifndef Minisat_ProofGenerator_h
#define Minisat_ProofGenerator_h

#define BIN_DRUP

#include "core/SolverTypes.h"

namespace Minisat {
    class ProofLogger {
    private:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // STATIC CONSTANTS
    #ifdef BIN_DRUP
        static int buf_len;
        static unsigned char drup_buf[];
        static unsigned char* buf_ptr;
    #endif

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PUBLIC API

        /// @brief Output file for DRUP proof
        FILE* drup_file = nullptr;

    public:
        ///////////////////////////////////////////////////////////////////////////////////////////
        // PROOF LOGGING

        bool enabled(void);

        template<class V>
        void addClause(const V& c);

        template<class V>
        void removeClause(const V& c);

        template <class V>
        void simplifyClause(const V& original, const V& simplified);

        void flush(void);

    #ifdef BIN_DRUP
        static void byteDRUP(Lit l);

        template<class V>
        static void binDRUP(unsigned char op, const V& c, FILE* drup_file);

        static void binDRUP_strengthen(const Clause& c, Lit l, FILE* drup_file);

        static inline void binDRUP_flush(FILE* drup_file);
    #else
        template<class V>
        static void asciiDRUP(unsigned char op, const V& c, FILE* drup_file);
    #endif
    };

    ///////////////////////////////////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF INLINE FUNCTIONS

    //////////////////
    // PROOF LOGGING

    inline bool ProofLogger::enabled(void) {
        return drup_file;
    }

    template <class V>
    inline void ProofLogger::addClause(const V& c) {
        if (drup_file) {
#ifdef BIN_DRUP
        binDRUP('a', c, drup_file);
#else
        for (int i = 0; i < c.size(); i++) {
            fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
        }
        fprintf(drup_file, "0\n");
#endif
        }
    }

    template <class V>
    inline void ProofLogger::simplifyClause(const V& original, const V& simplified) {
        addClause(simplified);
        removeClause(original);
    }

    template <class V>
    inline void ProofLogger::removeClause(const V& c) {
        // Do nothing if there is no output file
        if (!drup_file) return;

    #ifdef BIN_DRUP
        binDRUP('d', c, drup_file);
    #else
        asciiDRUP('d', c, drup_file);
    #endif
    }

    template <>
    inline void ProofLogger::removeClause(const Clause& c) {
        // Do nothing if there is no output file
        if (!drup_file) return;

        if (c.mark() != 1){
        #ifdef BIN_DRUP
            binDRUP('d', c, drup_file);
        #else
            asciiDRUP('d', c, drup_file);
        #endif
        } else {
            printf("c Bug. I don't expect this to happen.\n");
        }
    }

    inline void ProofLogger::flush(void) {
    #ifdef BIN_DRUP
        if (drup_file)
            binDRUP_flush(drup_file);
    #endif
    }

#ifdef BIN_DRUP
    inline void ProofLogger::byteDRUP(Lit l) {
        unsigned int u = 2 * (var(l) + 1) + sign(l);
        do{
            *buf_ptr++ = u & 0x7f | 0x80; buf_len++;
            u = u >> 7;
        }while (u);
        *(buf_ptr - 1) &= 0x7f; // End marker of this unsigned number.
    }

    template<class V>
    inline void ProofLogger::binDRUP(unsigned char op, const V& c, FILE* drup_file) {
        assert(op == 'a' || op == 'd');
        *buf_ptr++ = op; buf_len++;
        for (int i = 0; i < c.size(); i++) byteDRUP(c[i]);
        *buf_ptr++ = 0; buf_len++;
        if (buf_len > 1048576) binDRUP_flush(drup_file);
    }

    inline void ProofLogger::binDRUP_strengthen(const Clause& c, Lit l, FILE* drup_file) {
        *buf_ptr++ = 'a'; buf_len++;
        for (int i = 0; i < c.size(); i++)
            if (c[i] != l) byteDRUP(c[i]);
        *buf_ptr++ = 0; buf_len++;
        if (buf_len > 1048576) binDRUP_flush(drup_file);
    }

    inline void ProofLogger::binDRUP_flush(FILE* drup_file) {
        // fwrite(drup_buf, sizeof(unsigned char), buf_len, drup_file);
        fwrite_unlocked(drup_buf, sizeof(unsigned char), buf_len, drup_file);
        buf_ptr = drup_buf; buf_len = 0;
    }
#else
    template<class V>
    inline void ProofLogger::asciiDRUP(unsigned char op, const V& c, FILE* drup_file) {
        fprintf(drup_file, "%c ", op);
        for (int i = 0; i < c.size(); i++)
            fprintf(drup_file, "%i ", (var(c[i]) + 1) * (-2 * sign(c[i]) + 1));
        fprintf(drup_file, "0\n");
    }
#endif
}

#endif