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
#include "vertex.h"

/***************************
 * MISCELLANEOUS
 **************************/
// Print edge during simulation.  Useful for debugging
void vertex::PrintEdges() {
    vector <vertex *>::iterator it=_neighbours.begin();
    for (unsigned int i=0; i<_neighbours.size(); i++,it++) {
	cout << "... to (" << (*it)->GetX() << ", " << (*it)->GetY() << ", " << (*it)->GetZ() << "), type " << (*it)->GetType() 
	     << ", DE_static = "  << GetDE(i) << ", DC = " << GetDC(i) 
	     << ", J  = "  << _Js.at(i) << ", rate = " << _rates.at(i);
	     if ((*it)->IsOccupied()) cout << ", occupied";
	     cout << endl;
    } 
}
// 
void vertex::ClearDCs() {
    vector <double>::iterator it=_DCs.begin();
    for (; it!=_DCs.end(); it++) {
        (*it) = 0.0;
    }
}
//
void vertex::SetOccupied(const double & time) {
    _occupied = true;
    #ifdef printTotalOccupation
    _timeOfOccupation = time;
    _timesOccupied++;
    #endif
}
//
void vertex::SetUnoccupied(double time){
    _occupied = false;
    #ifdef printTotalOccupation
    _totalOccupationTime += time - _timeOfOccupation;
    #endif
}

/***************************
 * SETUP
 **************************/
//
void vertex::AddNeighbour(vertex *v, const double &J, const double &DE, const double &DZ, const double &RG) {
    if(find(_neighbours.begin(),_neighbours.end(),v)!=_neighbours.end()) {
        cout << "***ERROR***: Duplicated edges\n";
        exit(-1);
    }
    _neighbours.push_back(v);
    _Js.push_back(J);
    _RGs.push_back(RG);
    _DZs.push_back(DZ);
    _DEs.push_back(DE);
    _DCs.push_back(0.0);
    _rates.push_back(0.0);
}
//
void vertex::SetPos(const vec & pos){
    _pos = pos;
    _posZ = pos.getZ();
}
//
void vertex::SetType (string type){
    if (type!="c" && type!="g" && type!="-")
        ERROR(-1, "Don't understand vertex type " + type);
    _type = type;
    _electrode = (type != "-");
}
//
void vertex::SetE(double E) {
    _E = E;
}
// Modify deltaE's to reflect an applied field.
void vertex::ModifyDEsUsingField(const double & field) {
    vector <vertex *>::iterator it=_neighbours.begin();
    for (unsigned int i=0; i<_neighbours.size(); i++,it++) {
        _DEs.at(i) += field * _DZs.at(i);
    }
}
/*************************************
 * CALCULATE RATES 
 ************************************/
// Marcus hopping model
// When there are no 'hopperInteractions', can get away with simply calculating rates once:
void vertex::SetRates_DE(const double & kT) {
    double G;
    _totalRate = 0.0;
    for (unsigned int i=0; i<_neighbours.size(); i++) {
        G=_DEs.at(i) + _RGs.at(i);
        _rates[i] = ((_Js.at(i) * _Js.at(i) / hbar_eVs) * sqrt(pi / (_RGs.at(i) * kT))
                      * exp(-G * G / (4 * _RGs.at(i) * kT)));
        _totalRate += _rates[i];
    }
}
// Miller-Abrahams hopping model
// When there are no 'hopperInteractions', can get away with simply calculating rates once:
void vertex::SetRates_MA(const double& kT) {
    double DE;
    _totalRate = 0.0;
    for (unsigned int i = 0; i < _neighbours.size(); i++) {
        DE = _DEs.at(i);
        _rates[i] = _Js.at(i) * ((DE < 0.0) ? 1.0 : exp(-DE / kT));
        _totalRate += _rates[i];
    }
}
// Marcus hopping model
// When there *are* 'hopperInteractions', need to constantly update rates.
// This calculates the pre-factor in the Marcus expression 
//   (everything except the energetics)
void vertex::SetRatesPrefactor_C(const double & kT) {
    for (unsigned int i=0; i<_neighbours.size(); i++) { 
        _ratesPrefactor.push_back((_Js.at(i) * _Js.at(i) / hbar_eVs) * sqrt(pi / (_RGs.at(i) * kT)));
    }
}
// Miller-Abrahams hopping model
// When there *are* 'hopperInteractions', need to constantly update rates.
// This calculates the pre-factor in the Marcus expression 
//   (everything except the energetics)
void vertex::SetRatesPrefactor_CMA() {
    for (unsigned int i = 0; i < _neighbours.size(); i++) {
        _ratesPrefactor.push_back(_Js.at(i));
    }
}
// Whenever a charge is moved, all Coulombic energies need to be updated.
void vertex::IncrementDCs(int i, double C) {
    _DCs.at(i)+=C;
}
// Marcus hopping model
// Update the rates, given the updated _DCs.
void vertex::UpdateRates_C(const double & kT) {
    double G;
    _totalRate=0.;
    for (unsigned int i=0; i<_neighbours.size(); i++) { 
        G = _DEs.at(i) + _DCs.at(i) + _RGs.at(i);
        _rates[i]   = _ratesPrefactor[i] * exp(-G * G / (4.0 * _RGs.at(i) * kT));
        _totalRate += _rates[i];
    }
}
// Miller-Abrahams hopping model
// Update the rates, given the updated _DCs.
void vertex::UpdateRates_CMA(const double& kT) {
    double DE;
    _totalRate = 0.;
    for (unsigned int i = 0; i < _neighbours.size(); i++) {
        DE = _DEs.at(i) + _DCs.at(i);
        _rates[i] = _ratesPrefactor[i] * ((DE < 0.0) ? 1.0 : exp(-DE / kT));
        _totalRate += _rates[i];
    }
}

/************************************
 * CHOOSE DESTINATION OF HOPPER
 ***********************************/
// Choose the destination, assuming that all neighbours are unoccupied
// Return neighbour index
int vertex::ChooseNeighbour() const {
#ifdef RandomB
    double X = Uniform() * _totalRate;
#else
    double X = gsl_rng_uniform(gslRand) * _totalRate;
#endif
    for (unsigned int i = 0; i < _rates.size(); i++) {
        X -= _rates[i];
        if (X <= 0.) return i;
    }
    cout << "***ERROR***: ChooseNeighbour() in Vertex.h has not found anywhere to hop to (can't handle this yet!)\n";
    exit(-1);
    return -1;
}
// Choose the destination, but check the occupation of the neighbours first
//   (called only if an attempt is made to hop to an occupied vertex)
// Return neighbour index
int vertex::ChooseNeighbourUnoccupied(double totalRate) const {
#ifdef RandomB
    double X = Uniform() * totalRate;
#else
    double X = gsl_rng_uniform(gslRand) * totalRate;
#endif
    for (unsigned int i = 0; i < _rates.size(); i++) {
        if (!_neighbours[i]->IsOccupied()) {
            X -= _rates[i];
            if (X <= 0.) return i;
        }
    }
    cout << "***ERROR***: ChooseNeighbourUnoccupied() in Vertex.cc has not found anywhere to hop to (can't handle this yet!)\n";
    cout << scientific << "             X = " << X << endl;
    exit(-1);
    return -1;
}

// Recalculate the totalRate to unoccupied neighbours.
double vertex::CalcTotalRateToUnoccupied() {
    double totalRateToUnoccupied = 0.;
    for (unsigned int i=0; i<_neighbours.size(); ++i){
        if (!_neighbours[i]->IsOccupied()) totalRateToUnoccupied += _rates[i]; 
    }
    return totalRateToUnoccupied;
}

/*************************************************************
 * ANALYSIS
 * These functions are only called if 'printTotalOccupation'
 ************************************************************/
// Increment the Coulombic energy experienced by a hopper on this
//   vertex.
void vertex::IncrementEC(const double newEC, const double time) {
    _EC_time+= (time - _oldTime) * _EC;  // experienced old _EC for this time...
    _oldTime=time;
    _EC+=newEC;
}
// Set the Coulombic energy.  Used whenever a hopper is added / removed.
void vertex::SetEC(const double newEC, const double time) {
    _EC_time+= (time - _oldTime)*_EC;
    _oldTime=time;
    _EC=newEC;
}
//
void vertex::NormaliseTotalOccupationTime(const double maxTime, int totalHoppers) {
    _EC_time = _EC_time / _totalOccupationTime;
    _totalOccupationTime = _totalOccupationTime / (maxTime);
}

