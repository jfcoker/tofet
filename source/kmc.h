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

/**********************************************************************
 * 'kmc' is the central kinetic Monte Carlo class.
 * At the moment it contains two 'First Reaction Method' algorithms 
 * to cater for time-of-flight and variants, and FETs. 
 *********************************************************************/

#ifndef _KMC_H
#define	_KMC_H
#include "hoppers.h"
#include "graph.h"

using namespace std;

class kmc{
    private:
        graph * _graph;
        hoppers * _Hoppers;
        double _time;   // the simulation time 
        double _totalTimeOverAllRuns;   // the total time, summed over all runs 
        int _geometricBin;  // time bins for photocurrent transients (tof)
        int _nLogTimeBins;  // number of geometric time bins
        double _maxTime; 		
        int _run;  // the number of runs in each simulation
        double _maxRuns;  // set to double so that inf can be represented
        double _tol;  // tolerance of results of simulation
        double _lowerTol, _upperTol;
        double _dt;  // width of first time bin
        double _logDt;				
        double _alpha;  // subsequent log time bins are dt * [(alpha ^ n) - (alpha ^ (n-1))] wide
        double _logAlpha;
        double _sum_dz;  // the total distance moved along z, summed over all hoppers
        double _mu; // the mobility, from total displacement over total time
        vector <unsigned int> _hops; // The cumulative number of hops, seperated by reorganisation energy used
        vector <double> _current;  // photocurrent
        vector <int> _popgen_run;  // Hoppers in generation zone
        vector <int> _poptrans_run;  // Hoppers in transport zone
        vector <int> _popgen;  // Hoppers in generation zone
        vector <int> _poptrans;  // Hoppers in transport zone
        int _nHoppers;  // initial number of hoppers.  NOTE: this is not updated as hoppers are collected
        string _mode;  // mode of simulation (FET, tof, regenerate...)
        bool _hopperInteractions;  // Coulombic interactions?
    
        void UpdatePhotocurrent(const double &, const int& , const int&);
        void AveragePopOverRuns();
        double (hoppers::*moveFastest)();  // pointer to appropriate MoveFastest_* function

       /***************************************************
        * TIMEOUT
        **************************************************/
        int _timeoutMinutes;
        std::mutex _mutex;
        void SleepUntilTimeout();

    //end of private:

    public:
        kmc(){}
        kmc(char * sim, hoppers * Hoppers, int totalHoppers, graph * Graph, int timeoutMinutes){
            _totalTimeOverAllRuns = 0.0;
            _sum_dz = 0.0;
            _graph = Graph;
            _Hoppers = Hoppers;
            _hops = vector <unsigned int> (_graph->_reorgs.size(), 0);
            _maxTime=atof(Read(sim,"maxTime").c_str());
            _timeoutMinutes = timeoutMinutes;
            _mode = Read(sim, "mode", "tof");
            if (Read(sim, "hopperInteractions", "0") == "1") {
                _hopperInteractions = true;
            }
            else {
                _hopperInteractions = false;
            }
            if (_mode == "tof" || _mode == "regenerate" || _mode == "pb") {
                _dt=atof(Read(sim,"deltaTime").c_str());
                if (_dt > _maxTime) {
                    cout << "*** ERROR *** : deltaTime > maxTime!\n";
                    exit(-1);
                }
                _alpha=atof(Read(sim,"alpha").c_str());
                _nHoppers=totalHoppers;
                _tol=atof(Read(sim,"tol","0.0").c_str());
                _upperTol=1.0+_tol;
                _lowerTol=1.0-_tol;
                _maxRuns=atof(Read(sim,"maxRuns","inf").c_str()); 
                _logAlpha = log(_alpha);
                _logDt = log(_dt);
                _nLogTimeBins = int ( (log(_maxTime)  - _logDt ) / _logAlpha );
                _current.resize(_nLogTimeBins);
                _popgen_run.resize(_nLogTimeBins);
                _poptrans_run.resize(_nLogTimeBins);
                _popgen.resize(_nLogTimeBins);
                _poptrans.resize(_nLogTimeBins);
            }
            if (_mode=="tof") { 
                if (VERBOSITY_HIGH) {
                    cout << "Setting to MoveFastest_C" << endl;
                }
                moveFastest=&hoppers::MoveFastest_C;
            }
            if (_mode=="regenerate") {
                if (_hopperInteractions) { 
                    if (VERBOSITY_HIGH) {
                        cout << "Setting moveFastest to 'MoveFastest_RCI'" << endl;
                    }
                    moveFastest=&hoppers::MoveFastest_RCI;
                }
                else {
                    if (VERBOSITY_HIGH) {
                        cout << "Setting moveFastest to 'MoveFastest_R'" << endl;
                    }
                    moveFastest=&hoppers::MoveFastest_R;
                }
            }
            if (_mode == "pb") {
                if (VERBOSITY_HIGH) {
                    cout << "Setting to MoveFastest_PB" << endl;
                }
                moveFastest=&hoppers::MoveFastest_PB;
            }
            if (_mode=="fet") {
                if ( VERBOSITY_HIGH ) {
                    cout << "Setting moveFasest to 'MoveFastest_F'" << endl;
                }
                moveFastest=&hoppers::MoveFastest_F; 
                if (Read(sim,"converged","0")=="1") {
                    _Hoppers->SetActiveHoppersConverged();
                    cout << "Assuming the charge density is already converged\n";
                }
            }
        }
        ~kmc(){
            _current.clear();
        }

        /***************************************************
         * DO'S
         * TODO: These functions could probably be merged...
         **************************************************/
        void FRM();  // simple First Reaction Method	
        void FRM_FET();  // FRM with all necessary add-ons for FET simulations
        
        /***************************************************
         * GET'S
         **************************************************/
        const double & GetSumDz() const {return _sum_dz;}
        const double & GetTotalTimeOverAllRuns() const {return _totalTimeOverAllRuns;}
        const double & GetDt() const	{return _dt;}
        const double & GetTime() const 	{return _time;}
        const double & GetAlpha() const	{return _alpha;}
        vector <double> & GetTimeBins()	{return _current;}
        const int & GetnRuns() const	{return _run;}
        const double & GetMu() const {return _mu;}
        vector <unsigned int>& GetHops() { return _hops; }
        void PrintCurrent(string dest="") {
            double t1, t2;
            // If necessary, redirect 'cout' to 'fout'
            streambuf* cout_sbuf = std::cout.rdbuf();
            ofstream   fout;
            if(dest=="file") {
                fout.open("occVert.out");
                cout.rdbuf(fout.rdbuf());
            }
            for (int i = 0; i < _nLogTimeBins; i++) {
                    if (i == 0 ) {
                        t1 = 0;
                    }
                    else {
                        t1 = _dt * pow(_alpha, i-1);
                    }
                    t2 = _dt * pow(_alpha, i);
                    if (_current.at(i) != 0) {
                        cout << '\t' << t2 << "\t" << e * _current.at(i) / ((t2 - t1) * _graph->GetDepth()) <<endl;
                    }
            }
            if (dest=="file") {
                fout.close();
                // Restore the original stream buffer  
                cout.rdbuf(cout_sbuf);
            }
        }
        void PrintPops() {
            double t1, t2;

            for (int i = 0; i < _nLogTimeBins; i++) {
                if (i == 0) {
                    t1 = 0;
                }
                else {
                    t1 = _dt * pow(_alpha, i - 1);
                }
                t2 = _dt * pow(_alpha, i);
                if (_current.at(i) != 0) { // Only print pop for timebins where at least one hop occured. If no hops occured, UpdatePhotocurrent() would not have been called, and no populations would have been stored.
                    cout << '\t' << t2 << "\t" << _popgen.at(i) << "\t" << _poptrans.at(i) << "\t" << _run * _nHoppers - (_popgen.at(i) + _poptrans.at(i)) << endl;
                }
            }
        }
    // end of public:
};
#endif	/* _KMC_H */

