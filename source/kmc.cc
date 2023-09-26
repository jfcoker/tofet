///////////////////////////////////////////////////////////////////////
//  This file is part of ToFeT.
//  
//  ToFeT is free software: you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published by
//  the Free Software Foundation, either version 3 of the License, or
//  (at your option) any later version.
//  
//  ToFeT is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//  
//  You should have received a copy of the GNU Lesser General Public License
//  along with ToFeT.  If not, see <http://www.gnu.org/licenses/>.
///////////////////////////////////////////////////////////////////////
#include "kmc.h"
// TODO: merge these functions?
// TODO: Make sure timeout functionality works correctly on different system types.
// TODO: Display warning if simulation time exceeded when this is not the expected conditon for simulation finishing.
// Simple First Reaction Method
void kmc::FRM() {
    _mu = 1e50;
    double prevMu = 1e50;
    double changeInMu = 1e50;
    double dz;
    int hopReorgEnum;

    bool interrupted = false;
    if (_timeoutMinutes) {
        thread timeoutThread(&kmc::SleepUntilTimeout, this);
        timeoutThread.detach();
    }
    
    _run = 0;
    while (!interrupted) {  // entire simulation...

        if (_run >= _maxRuns) {
            cout << "!!! WARNING !!! : Mobility not converged, maxRuns reached.\n";
            WARNINGS++;
            break;
        }
        else _run++;

        _Hoppers->GenerateAll(_nHoppers, 0.0);
        if (_hopperInteractions) _Hoppers->SetHops_C(0.0);
        _Hoppers->FindFastest();
        _time=0.0;
        while (_Hoppers->GetActive()>0) {  // single run...
            _time = _Hoppers->GetFastestTime();
            hopReorgEnum = _Hoppers->GetFastestReorgEnum();
            if (hopReorgEnum >= 0) _hops[hopReorgEnum]++;
            dz = (_Hoppers->*moveFastest)();
            _sum_dz+=dz;

            auto pop = _Hoppers->GetPop();
            int popgen, poptran;
            tie(popgen, poptran) = pop;
            UpdatePhotocurrent(dz,popgen,poptran);
            if (_time > _maxTime) break;
            if (_timeoutMinutes) {
                if (_mutex.try_lock()) {
                    // Try to gain ownership of the mutex object
                    // This should only succeed if the timeout thread has released it upon reaching the end of the timout interval.
                    cout << "!!! WARNING !!! : Timeout triggered, ending KMC...\n";
                    WARNINGS++;
                    interrupted = true;
                    _mutex.unlock();
                    break;
                }
            }
            if (RECEIVED_TERM_SIGNAL) {
                cout << "!!! WARNING !!! : Received interrupt or terminate signal, ending KMC...\n";
                WARNINGS++;
                interrupted = true;
                break;
            }
        }
        _totalTimeOverAllRuns += _time;
        
        AveragePopOverRuns();

        prevMu = _mu;
        _mu = _sum_dz * 1e-16 / (_totalTimeOverAllRuns * _nHoppers * -_graph->GetFieldZ());
        changeInMu = _mu / prevMu;

        if (VERBOSITY_HIGH) {
            cout << "Run number " << _run
                 << ": Hoppers left = " << _Hoppers->GetActive()
                 << "; Mobility (cm^2/V.s) = " << _mu;
            if (_run > 1) cout << "; Fractional change of mob. = " << changeInMu;
            cout << endl << flush;
        }

        if (_run > 1 && changeInMu > _lowerTol && changeInMu < _upperTol) {
            cout << "Mobility converged" << endl;
            break;
        }

        _Hoppers->SetWaitTimes(_time);
        _Hoppers->softClear();
        _graph->ClearDCs();
    }
}
// First reaction method with all the necessary ancillary functions to handle FETs
void kmc::FRM_FET() {
    _Hoppers->SetHops_C(0.0);
    _Hoppers->FindFastest();
    _time=0.0;
    while ( _Hoppers->_run ) {
        _time  = _Hoppers->GetFastestTime();
        (_Hoppers->*moveFastest)();
    }
    _Hoppers->SetWaitTimes(_time);
}
// 
void kmc::UpdatePhotocurrent(const double & dz, const int& gen, const int& trans){
    _geometricBin = int ((log(_time) - _logDt) / _logAlpha);
    if ( _geometricBin<0 ) _geometricBin=0;
    else if (_geometricBin >= _nLogTimeBins) {
        _nLogTimeBins = _geometricBin + 1;
        _current.resize(_nLogTimeBins);

        _popgen_run.resize(_nLogTimeBins);
        _poptrans_run.resize(_nLogTimeBins);

        _popgen.resize(_nLogTimeBins);
        _poptrans.resize(_nLogTimeBins);
    }
    _current[_geometricBin] += dz;
    _popgen_run[_geometricBin] = gen;
    _poptrans_run[_geometricBin] = trans;
}

void kmc::AveragePopOverRuns() {
    for (int i = 0; i < _geometricBin - 1; i++) {
        _popgen.at(i) += _popgen_run.at(i);
        _poptrans.at(i) += _poptrans_run.at(i);
    }
}


void kmc::SleepUntilTimeout() {
    _mutex.lock();
    this_thread::sleep_for(chrono::minutes(_timeoutMinutes));
    _mutex.unlock();
}

