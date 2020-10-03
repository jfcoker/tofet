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

/*********************************************************************
 * 'graph' is the class that defines the graph through which hoppers 
 * (charges) move.  The graph is constituted of vertices (molecules) 
 * and edges (connections between molecules).
 *********************************************************************/
#ifndef _GRAPH_H
#define	_GRAPH_H
#include "vertex.h"
#include "global.h"
#include "IO.h"

class graph{
    private:
        vector <vertex *> _vertices;
        double _Vg;  // V.Ang^-1 (sorry!)
        double _fieldZ;  // V.Ang^-1
        double _temp;  // K
        double _sizeX, _sizeY, _sizeZ; // Angs
        bool _applyPBs;
        vector <vector <double> > _CoulombGrid;
        bool _hopperInteractions; 
        double _tmpX, _tmpY, _tmpZ;
    // end of private:
    
    public:
        double _reorg; // eV
        double _kT;    // eV
        double _sourceFermiEnergy;
        double _drainFermiEnergy;
        double _coulombPrefactor;

        graph(){}

        graph(char * sim, char *xyz, char *edge){

            if (Read(sim, "mode", "tof") != "fet") 
                _fieldZ = atof(Read(sim, "fieldZ").c_str()); 
            else {
                _fieldZ = 1e50;
                _Vg = atof(Read(sim, "Vg").c_str()); 
                double Vds = atof(Read(sim, "Vds").c_str());
                _sourceFermiEnergy = _Vg;
                _drainFermiEnergy  = _Vg + Vds;
                cout << "Source Fermi energy = " << _sourceFermiEnergy
                     << ", drain Fermi energy = " << _drainFermiEnergy << endl;
            }

            _reorg = atof(Read(sim, "reorg").c_str()); 
            _temp = atof(Read(sim, "temp").c_str());
            _kT = _temp*k_eVK;

            _applyPBs = (Read(sim, "mode", "tof") == "pb");
            _hopperInteractions = (Read(sim, "hopperInteractions", "0") == "1");

            if (_applyPBs) {
                if (_hopperInteractions) {
                    cout << "Read simulation volume sizeX, sizeY, sizeZ ...\n";
                    _sizeX = atof(Read(sim, "sizeX").c_str());
                    _sizeY = atof(Read(sim, "sizeY").c_str());
                    _sizeZ = atof(Read(sim, "sizeZ").c_str());
                } 
                else {
                    cout << "Read simulation volume sizeZ ...\n";
                    _sizeZ = atof(Read(sim, "sizeZ").c_str());
                }
            }

            // If in FET mode or _hopperInteractions enabled, attempt to read site energies from .xyz, even if siteEnergies option is missing from .sim 
            bool readSiteEnergies = (Read(sim, "siteEnergies", "0") == "1" || Read(sim, "mode", "tof") == "fet" || _hopperInteractions);
            if (VERBOSITY_HIGH) {
                if (readSiteEnergies) cout << "Reading E's from ***.xyz\n";
                else cout << "Reading delta E's from ***.edge\n";
            }

            // Read input files, grabbing site energies from .xyz, or delta Es from .edge, as requested.
            // If reading site energies, calculate delta Es here as well.
            ReadVertices(xyz, _vertices, readSiteEnergies);
            ReadEdges(edge, _vertices, !readSiteEnergies);

            ModifyDEsUsingField();

            if (_hopperInteractions || Read(sim, "mode", "tof") == "fet") {

                SetRatesPrefactor_C();  // doesn't set the field!
                // The Coulomb prefactor in eV.Ang/e^2:
                _coulombPrefactor = 14.3996442 / atof(Read(sim, "dielectric").c_str());
                // The following should be uncommented if you want to use a look-up 
                // table for the Coulombic interactions.  See also 'GetSingleCoulomb'
                // in hoppers.cc
                // MakeCoulombEnergyGrid(); 
            }
            else SetRates_DE();

            if (Read(sim, "printVertices", "0") == "1") PrintVertices(readSiteEnergies);
            if (Read(sim, "printEdges", "0") == "1") PrintEdges();
        }

        ~graph(){
            vector < vertex* >::iterator it = _vertices.begin();
            for (; it!=_vertices.end(); it++)
                delete (*it);

            _vertices.clear();
        }

    /*****************************
     * SETS and ADDS
     ****************************/
    void ClearDCs(); // reset _DCs to 0
    //void AddEdge(const unsigned int &, const unsigned int &, const double &);
    void ModifyDEsUsingField();  // ... likewise for deltaE's
    void SetRatesPrefactor_C();  // set prefactors (no energies) 
    void SetRates_DE();  // set rates from deltaE's 
    void NormaliseOccupationTimes(const double, int);  
    void MakeCoulombEnergyGrid();  
    double const &GetCoulomb(vertex *, vertex *);  // ... from a grid

    /*****************************
     * PRINTS AND READS 
     ****************************/
    void ReadEdges(char *, vector <vertex *> &, bool);
    void ReadVertices(char *, vector <vertex *> &, bool);
    void PrintEdges();
    void PrintVertices(bool);
    void PrintEnergies();  // print average sum of static + coulomb energies
    void PrintOccupied();  // print all occupied molecules	
    void PrintTotalOccupationTimes();

    /*****************************
     * GETS 
     ****************************/
    vector <vertex *> GetPreviouslyOccupied(char *);  // read occupied vertices from file
    vertex * GetEmptyGenerator();  // returns random empty generator
    double GetDepth();  // get the depth of the graph in the z direction
    double GetDistance(vertex *, vertex *);  // get the distance between two vertices
    int CountTotalElectrodes();
    const double & GetFieldZ() 	const {return _fieldZ;}
    vector <vertex *> GetCollectors();
    vector <vertex *> GetGenerators(); 
};
#endif	/* _GRAPH_H */
