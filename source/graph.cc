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
#include "graph.h"

// 1D minimum image distance
// Calculates shortest distance between two points, taking into account periodic boundaries at [0,size]
// For seperations of up to 1.5*size, should remain correct even if one or both points are outside the boundaries.
double min_img_dist(double s1, double s2, double size) {
    double delta = s2 - s1;
    int k = (int)delta * (2.0 / size);
    delta -= k * size;

    k = (int)delta * (2.0 / size);
    delta -= k * size;

    return delta;
}

/************************************************************
 * SET-UP GRAPH
 * Read in from ***.xyz and ***.edge files and generate graph
 ************************************************************/
// Read from ***.xyz
void graph::ReadVertices(char * filename, vector <vertex *> &vertices, bool readEnergies=false) {
    ifstream in;
    open(filename, in);
    string word;
    double x, y, z, E;
    string type;
    int counter=0;

    while (in) {
        in >> word; x=atof(word.c_str());
        in >> word; y=atof(word.c_str());
        in >> word; z=atof(word.c_str());
        in >> word; type=word;
        if (readEnergies) in >> word; E=atof(word.c_str()); 

        if (!in) break;

        vertex * newVertex = new vertex;
        vec pos(x,y,z);
        newVertex->SetPos(pos);
        newVertex->SetType(type);
	    newVertex->SetID(counter);
        if (readEnergies) newVertex->SetE(E);

        vertices.push_back(newVertex);
        counter++;
    }
    cout << "Read in " << counter << " vertices from " << filename << "\n";
    in.close();
}
// Read from ***.edge
void graph::ReadEdges(char *filename, vector <vertex *> &vertices, bool readDeltaEnergies=true) {
    if (vertices.size()<1) {
        cout << "***ERROR***: You are trying to initialise edges before vertices\n";
        exit(-1);
    }
    ifstream in;
    open(filename, in);
    string word;
    unsigned int v1, v2;
    double J, DE, DZ;
    int counter=0;
    while (in) {
        in >> word; v1=atoi(word.c_str());
        in >> word; v2=atoi(word.c_str());
        in >> word; J =atof(word.c_str());

        if (readDeltaEnergies) { in >> word; DE=atof(word.c_str()); }
        else DE = vertices[v2]->GetE() - vertices[v1]->GetE();
        
        if (_applyPBs) DZ = min_img_dist(vertices[v1]->GetZ(), vertices[v2]->GetZ(), _sizeZ);
        else DZ = vertices[v2]->GetZ() - vertices[v1]->GetZ();

        if (!in) break;

        if (v1 >= vertices.size() || v2 >= vertices.size()  || v1 == v2 ){
            cout << "***ERROR***: Trying to create an edge on non-existent vertex "
                 << v1 << "->" << v2 <<"\n";
            exit(-1);
        }

        vertices[v1] -> AddNeighbour(vertices[v2],J, DE, DZ);
        vertices[v2] -> AddNeighbour(vertices[v1],J,-DE, -DZ);
        
        counter++;
    }
    in.close();
    cout << "Read in " << counter << " edges from " << filename << "\n";
}

/***************************************************************
 * SET THE ENERGETICS AND RATES OF THE EDGES OF THE GRAPH 
 * Most of these functions simply wrap counterparts in vertex.cc
 **************************************************************/
// Modify DE's to reflect an applied field.
void graph::ModifyDEsUsingField() {
    vector <vertex *>::iterator it=_vertices.begin();
    for (; it!=_vertices.end(); it++)
	    (*it)->ModifyDEsUsingField(_fieldZ);
}
// Marcus hopping model
// When Coulombic interactions are enabled, only pre-factor 
//   is constant.
void graph::SetRatesPrefactor_C() {
    if (VERBOSITY_HIGH) cout << "Setting rates pre-factors\n";

    vector <vertex *>::iterator it=_vertices.begin();
    for (; it!=_vertices.end(); it++)
        (*it)->SetRatesPrefactor_C(_reorg,_kT);
}
// Marcus hopping model
// Without Coulombic interactions rates are constant
//   and can just be set once
void graph::SetRates_DE() {
    if (VERBOSITY_HIGH) cout << "Setting rates using DEs\n"; 

    vector <vertex *>::iterator it=_vertices.begin();
    for (; it!=_vertices.end(); it++)
        (*it)->SetRates_DE(_reorg,_kT);
}
// Miller-Abrahams hopping model
// When Coulombic interactions are enabled, only pre-factor 
//   is constant.
void graph::SetRatesPrefactor_CMA() {
    if (VERBOSITY_HIGH) cout << "Setting rates pre-factors\n";

    vector <vertex*>::iterator it = _vertices.begin();
    for (; it != _vertices.end(); it++)
        (*it)->SetRatesPrefactor_CMA();
}
// Miller-Abrahams hopping model
// Without Coulombic interactions rates are constant
//   and can just be set once
void graph::SetRates_MA() {
    if (VERBOSITY_HIGH) cout << "Setting rates using DEs\n";

    vector <vertex*>::iterator it = _vertices.begin();
    for (; it != _vertices.end(); it++)
        (*it)->SetRates_MA(_kT);
}
// Set all difference in Coulomb energies to 0.0
void graph::ClearDCs() {
    vector <vertex *>::iterator it=_vertices.begin();
    for (; it!=_vertices.end(); it++)
        (*it)->ClearDCs();
}
// Construct a lookup table for Coulombic interactions
// See hoppers.cc for discussion of the relative merits of this approach
void graph::MakeCoulombEnergyGrid() { 
    vector <double > tmp;
    vector <vertex *>::iterator it_i = _vertices.begin();
    for (; it_i != _vertices.end(); ++it_i) {
        vector <vertex *>::iterator it_j = _vertices.begin();
        tmp.clear();
        for (; it_j != _vertices.end(); ++it_j) {
            if ((*it_i) == (*it_j))
                tmp.push_back(0.0);
            else
                tmp.push_back(_coulombPrefactor/GetDistance((*it_i), (*it_j))); 
        }
	_CoulombGrid.push_back(tmp);
    }
} 
// Look-up from the energy grid (above)
double const &graph::GetCoulomb(vertex * v1, vertex * v2) {
    //return _CoulombGrid.at(v1->GetID()).at(v2->GetID());	// safe but slow
    return _CoulombGrid[v1->GetID()][v2->GetID()];  // unsafe but fast (10% speed-up)
}

/****************************************
 * OUTPUT 
 ***************************************/
// Print edges. If 'hopperInteractions' there are no DE's or rates yet.
void graph::PrintEdges() {
    cout << "mol_ID\tneigh_ID\tdelta_z\tJ";
    if (!_hopperInteractions) cout << "\tDE\trate";
    cout << endl;

    for (unsigned int i=0; i<_vertices.size(); i++) {
        cout << i << endl;
        vector <vertex *> Neighbours = (*_vertices.at(i)).GetNeighbours();
        vector <vertex *>::iterator it=Neighbours.begin();
        int neigh=0;
        for (; it<Neighbours.end(); it++,neigh++) {
            for (unsigned int j=0; j<_vertices.size(); j++) {
                if (_vertices.at(j) == *it) {      // This is baroque!
                    cout << '\t' << j << '\t';
                    cout << _vertices.at(i)->GetDZ(neigh);
                    cout << '\t' << _vertices.at(i)->GetJ(neigh);
                    if (!_hopperInteractions) {  
                        cout << '\t' << _vertices.at(i)->GetDE(neigh);
                        cout << '\t' << (*_vertices.at(i)).GetRate(neigh);
                    }
                    cout << endl;
                    break;
                }
            }
        }
    }
}
// Print vertices, including site energies if necessary.
void graph::PrintVertices(bool printSiteEnergies){
    vector <vertex *>::iterator it=_vertices.begin();
    for (; it!=_vertices.end(); it++) {
        cout << (*it)->GetX() << '\t' << (*it)->GetY() << '\t' << (*it)->GetZ() << '\t';
        
        if ((*it)->IsCollector()) cout << "c";
        else if ((*it)->IsGenerator()) cout << "g";
        else cout << "-";

        if (printSiteEnergies) cout << '\t' << (*it)->GetE();
        cout << endl;
    }
}
// Print energies.  The Coulomb contribution is averaged over the total
//   time that the vertex is occupied.
void graph::PrintEnergies() {
    cout << "> ENERGY (static + Coulomb)\n"
         << "\tx (Ang)\ty (Ang)\tz (Ang)\tE (eV)\n";
    vector <vertex *>::iterator it_outer = _vertices.begin();
    for (; it_outer!=_vertices.end(); it_outer++) {
        cout << '\t' << (*it_outer)->GetX() << '\t'
             << (*it_outer)->GetY() << '\t' 
             << (*it_outer)->GetZ() << '\t'
             << (*it_outer)->GetE() + (*it_outer)->GetEC_time() << endl;
    }
}
// Simply print all vertices that are occupied.
void graph::PrintOccupied(){
    for (unsigned int i =0; i < _vertices.size(); ++i)
        if(_vertices[i]-> IsOccupied()) cout << i << '\t' << _vertices[i] << endl;
}
//
void graph::PrintTotalOccupationTimes() {
    cout << "> TOTAL OCCUPATION TIMES AND TIMES VISITED\n"
    << "\tx (Ang)\ty (Ang)\tz (Ang)\ttime (fraction of maxTime)\ttimes visited\n";
    for (unsigned int i =0; i < _vertices.size(); ++i){
         cout << '\t' << _vertices[i]->GetX() << '\t'
              << _vertices[i]->GetY() << '\t' << _vertices[i]->GetZ() << '\t'
              << _vertices[i]->GetTotalOccupationTime() << "\t\t"   
              << _vertices[i]->GetTimesOccupied() << endl;
    }
}

/*******************************
 * MISCELLANEOUS
 ******************************/
// Get the distance between two vertices
double graph::GetDistance(vertex * v1, vertex * v2) {
    if (_applyPBs) {
        _tmpX = min_img_dist((*v1)._pos._x, (*v2)._pos._x, _sizeX);
        _tmpY = min_img_dist((*v1)._pos._y, (*v2)._pos._y, _sizeY);
        _tmpZ = min_img_dist((*v1)._pos._z, (*v2)._pos._z, _sizeZ);
    }
    else {
        _tmpX = (*v2)._pos._x - (*v1)._pos._x;
        _tmpY = (*v2)._pos._y - (*v1)._pos._y;
        _tmpZ = (*v2)._pos._z - (*v1)._pos._z;
    }
    return sqrt(_tmpX * _tmpX + _tmpY * _tmpY + _tmpZ * _tmpZ);
}
// 
int graph::CountTotalElectrodes() {
    int total=0;
    vector <vertex *>::iterator it=_vertices.begin();
    for(;it!=_vertices.end(); ++it) {
        if ( (*it)->_electrode ) total++;
    }
    return total;
}
// Read vertex ID's from 'filename'
//   and return a vector of vertices
vector <vertex *> graph::GetPreviouslyOccupied(char * filename) {
    vector <vertex *> generateOnMe;
    ifstream inFile;
    inFile.open(filename);
    if (!inFile) return generateOnMe;

    string word;
    unsigned int v;
    while (true) {
        inFile >> word;
        if (!inFile) break;

        v = atoi(word.c_str());
        if (v > _vertices.size()-1 || v < 0) {
            cout << "***ERROR*** Don't understand vertex " 
                 << v << " in inFile.out\n";
            exit(-1);
        }
        generateOnMe.push_back(_vertices.at(v));
    }
    return generateOnMe; 
}
// Return a vector of generators
vertex * graph::GetEmptyGenerator(){
    vector <vertex *>::iterator it_vert;
    vector <vertex *> plausibleCandidate;
    for (it_vert = _vertices.begin(); it_vert != _vertices.end() ; ++it_vert ){
        if ( !((*it_vert)->IsOccupied()) && (*it_vert)->IsGenerator() )
            plausibleCandidate.push_back(*it_vert);
    }

    // Make sure we pick a generator with at least 1 neighbour.
    // Otherwise we will get an infinite waitTime until next hop.
    while(plausibleCandidate.size() > 0) {

        #ifdef RandomB
        it_vert = plausibleCandidate.begin() + RandPos(plausibleCandidate.size()); //samples in range [0,n-1]
        #else
        it_vert = plausibleCandidate.begin() + gsl_rng_uniform_int(gslRand, plausibleCandidate.size()); //samples in range [0,n-1]
        #endif
        
        // Found a suitable candidate.
        if ((*it_vert)->GetNumberNeighbours() > 0)
            return *it_vert;

        // Not suitable, remove from contention.
        plausibleCandidate.erase(it_vert);

    }

    cout << "***ERROR*** Failed to find an empty generator with at least one edge!\n"
         << "            Size of _vertices=" << _vertices.size() << endl;
    exit(-1);

}
// Return a vector of collectors
vector <vertex *> graph::GetCollectors() {
    vector <vertex *> collectors;
    vector <vertex *>::iterator it = _vertices.begin();
    for (; it!=_vertices.end(); ++it) {
        if ( (*it)->IsCollector() ) collectors.push_back(*it);
    }
    return collectors;
}
// Return a vector of generators
vector <vertex *> graph::GetGenerators() {
    vector <vertex *> generators;
    vector <vertex *>::iterator it = _vertices.begin();
    for (; it!=_vertices.end(); ++it) {
        if ( (*it)->IsGenerator()) generators.push_back(*it);
    }
    return generators;
}
// Get the depth of the graph along the 'z' axis
double graph::GetDepth() {
    static double depth = -1.0;

    // We only need to calculate this the first time it is called.
    if (depth < 0.0) {
        double zMin = 1e50;
        double zMax = -1e50;
        vector <vertex*>::iterator it = _vertices.begin();
        for (; it != _vertices.end(); ++it) {
            if ((*it)->GetZ() < zMin) zMin = (*it)->GetZ();
            if ((*it)->GetZ() > zMax) zMax = (*it)->GetZ();
        }
        depth = zMax - zMin;
    }

    return depth;
}
// Wrap the same function in vertex.cc
void graph::NormaliseOccupationTimes(double maxTime, int totalHoppers ){
    vector <vertex *>::iterator it_all = _vertices.begin();
    for (; it_all!=_vertices.end(); it_all++)
        (*it_all) -> NormaliseTotalOccupationTime( maxTime, totalHoppers );
}
