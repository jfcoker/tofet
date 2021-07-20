#include "global.h"

#ifndef RandomB
gsl_rng * gslRand;
#endif

bool VERBOSITY_HIGH = false;
int WARNINGS = 0;
bool RECEIVED_TERM_SIGNAL = false;

void ERROR(int code, std::string msg) {
    std::cout << "!!! ERROR !!!: " << msg << std::endl;
    exit(code);
}

void signal_handler(int s) {
    RECEIVED_TERM_SIGNAL = true;
}

template<typename T> void safe_delete(T& obj) {
    typename T::iterator it = obj.begin();
    for (; it != obj.end(); ++it) {
        delete* it;
    }
    obj.clear();
}