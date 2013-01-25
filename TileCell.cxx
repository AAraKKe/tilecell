#include <sstream>
#include <iostream>
#include <cmath>
#include <cstdlib>

#include "TROOT.h"
#include "TH1D.h"

#include "TileCell.h"

// Declaraion of static members
bool TileCell::s_isSet = false;
double TileCell::s_maxEnergy = 0;
double TileCell::s_minEnergy = 99999999;
int TileCell::s_nBins;
double TileCell::s_emin;
double TileCell::s_emax;

using namespace std;
/**
  @brief Default constructor of the class. Defines information and initialize variables

  @param      partition   The partition of the TileCell.
  @param      layer       Layer of the TileCell.
  @param      tower       Tower (division of 0.1 in eta) of the TileCell.
  @param      module      Module (one of the 64 division in phi) of the TileCell.
  @param      pmt1        Status of the PMT-1 of the TileCEll.
  @param      pmt2        Status of the PMT-2 of the TileCEll.

  The partition is numbered from 1 to 4 as

    -1 = LBA
    -2 = LBC
    -3 = EBA
    -4 = EBC

  Partions LBA and EBA correspond to positive values of eta while LBC and EBC correspond to
  negative values of eta.

  Tile TileCal is divided horizontally in layers (samples). These layers are numbered from 0 to 15 corresponding
  to 

  0 = A
  1 = B = BC = C
  2 = D
  3 = Special gap scin cells E1-E4
  >3 = individual tile used to Cx calibration.

  The tower variable is a index which tell the 0.1 division in eta where the cell is. It goes from 
  0 to 15 (eta (0,1.6) and the sign of eta is know with the partition.

  The module number (from 1 to 64) correspond to the number of division in phi.

  The status of the PMTs for each TileCell are coded as:

  0 = Good
  1 = Noisy
  2 = Some other problems
  3 = Bad
  */    
TileCell::TileCell(int partition, int layer, int tower, int module, int pmt1, int pmt2){

    /// Set the code (bitwise) to the cell. It saves sapace and memory access.
    setCode(partition, layer, tower, module);

    /// Set the status of the TileCell
    if (pmt1 == 0 && pmt2 == 0) m_isGoodCell = true;
    else m_isGoodCell = false;

    //. Initialize energy and GeV values
    m_isCalculated = false;
    /// Resize vectors
    if (!s_isSet) {
        cerr<<"ERROR <TileCell::TileCell()>: Before creating any cell you must specify the histogram properties for the energy. Use the static function setHistoInfo."<<endl;
        cerr<<"TileCell::setHistoInfo(nbins,xmin,xmax)"<<endl;
        exit(1);
    }

    /// Initialize
    /*for (int i = 0; i < s_nBins+2; i++){
      m_energy.push_back(0);
      m_error.push_back(0);
      }*/

    m_energy = NULL;
    m_error = NULL;

}

/**
  @brief The destructor of the class.
  */    
TileCell::~TileCell(){

}

/**
  @brief Sets the information for the energy histogram in the cell.

  @param  nbins   Number of bins for the histogram
  @param  emin    Minimum energy for the histogram
  @param  emax    Maximum energy for the histogram

*/    
void TileCell::setHistoInfo(int nbins, double emin, double emax){

    /// Set information   
    s_nBins = nbins;
    s_emin = emin;
    s_emax = emax;

    /// Tell the code the histogram has been set
    s_isSet = true;

    // Prompt information about cells
    cout<<"-----------------------------------------------------"<<endl;
    cout<<"The energy histogram for each cell has been defined:"<<endl;
    cout<<"Number of bins: \t"<<s_nBins<<endl;
    cout<<"Limits: \t("<<s_emin<<","<<s_emax<<")"<<endl;
    cout<<"-----------------------------------------------------"<<endl;

}

/* @brief This function prints the information about highest and lowest energies in case they are over or under
 *        the limits of the histograms of each cells
 */
void TileCell::printEnergyInfo(){
    if(s_maxEnergy > s_emax) cerr<<"WARNING: The upper limit for energy in the cells is set to "<<s_emax<<" but an energy of "<<s_maxEnergy<<" has been found while scanning. Change the upper limit or some entries can be lost in the cells"<<endl;
    if(s_minEnergy < s_emin) cerr<<"WARNING: The lower limit for energy in the cells is set to "<<s_emin<<" but an energy of "<<s_minEnergy<<" has been found while scanning. Change the lower limit or some entries can be lost in the cells"<<endl;
}

/**
  @brief This fuction store the energy of a cell into a vector. It choose the right component thinking
  in the vector as an histogram.

  @param      energy Is the value of the energy measured in the TileCell.
  @param      weight Is the weight it will have in the case some reweighting is applied.
  */    
void TileCell::setEnergy(double energy, double weight){

    /// Check if the info has been set
    if (!s_isSet) {
        cerr<<"ERROR <TileCell::TileCell()>: Before creating any cell you must specify the histogram properties for the energy. Use the static function setHistoInfo."<<endl;
        cerr<<"TileCell::setHistoInfo(nbins,xmin,xmax)"<<endl;
        exit(1);
    }
    if (m_energy == NULL){
        m_energy = new double[s_nBins+2];
        m_error =  new double[s_nBins+2];
        memset(m_energy,0,(s_nBins+2)*sizeof(double));
        memset(m_error,0,(s_nBins+2)*sizeof(double));
    }


    /// Derive energy bin taking into accouint the over and underflow
    bool inside = (energy - s_emin)*(energy - s_emax) <= 0;
    int bin = (inside) ? int((energy-s_emin)/(s_emax-s_emin)*s_nBins)+1 : (energy < s_emin) ? 0 : s_nBins + 1;

    /// Fill the vector with energy values
    /*cout<<"BIN: "<<m_code<<" "<<bin<<endl;
    cout<<"Value: "<<m_energy[bin]<<endl;
    cout<<"Peso: "<<weight<<endl;*/
    m_energy[bin] += weight;
    /*cout<<"Value after: "<<m_energy[bin]<<endl;*/

    /// Fill the vector with the errors. 
    m_error[bin] += weight*weight;

    ///To finish, update de max and min values for log.
    if (energy > s_maxEnergy) s_maxEnergy = energy;
    if (energy < s_minEnergy) s_minEnergy = energy;
}

/**
 *  @brief Return the mean energy value of the cell
 */
double TileCell::getEnergy(){
    if (m_isCalculated) return m_meanEnergy;
    setEnergyNoise();
    return m_meanEnergy;
}

/**
 *  @brief Return the error in the mean energy of the cell
 */
double TileCell::getEnergyError(){
    if (m_isCalculated) return m_energyError;
    setEnergyNoise();
    return m_energyError;
}

/**
 *  @brief Return the moise of the cell, computed as the RMS of the histogram of energy deposited in the cell
 */
double TileCell::getNoise(){
    if (m_isCalculated) return m_noise;
    setEnergyNoise();
    return m_noise;
}

/**
 *  @brief Return the error of the noise of the cell
 */
double TileCell::getNoiseError(){
    if (m_isCalculated) return m_noiseError;
    setEnergyNoise();
    return m_noiseError;
}

/**
 * @brief Returns the status of the cell as good (true) or dead (false)
 */
bool TileCell::isGood(){
    return m_isGoodCell;
}

/**
  @brief Returns a string with the values of the desired variables of the cell (partition, layer...). This can be handy when using the partition, layer... variable as names for histograms for example.

  @param  set   This is a member of the KeyString enum.

  - P, L, T or M: If this variable is one of these four parameters it will return the partition, layer, tower or module as a string.
  - PL: In this case the string "partition_layer" will be returned.
  - PLT: In this case the string "partition_layer_tower" will be returned.
  - PLTM: In this case the string "partition_layer_tower_module" will be returned.
  - Default: The default value is PL.

  @param separator It defines how each value will be separated. By default it is underscore "_".

  @param _bar If _bar is set to true at the end of each string an extra separator will be added if it want to be used as prefix. The default value is true.
  */    
string TileCell::getKey(KeyString set, string separator, bool _end){
    stringstream key;

    if (set == P) key << getPartition();
    if (set == L) key << getLayer();
    if (set == T) key << getTower();
    if (set == M) key << getModule();
    if (set == PL) key << getPartition() << separator << getLayer();
    if (set == PLT) key << getPartition() << separator << getLayer() << separator << getTower();
    if (set == PLTM) key << getPartition() << separator << getLayer() << separator << getTower() << separator << getModule();
    if (_end) key << separator;
    return key.str();
}

/**
 * @brief This function return the unique code of the cell
 */
int TileCell::getCode(){
    return m_code;
}

/**
 * @brief Return the value of the partition in case it is needed.
 */
int TileCell::getPartition(){

    return m_code & 0xFF;
}

/**
 * @brief Return the layer number in case it is needed.
 */
int TileCell::getLayer(){
    return (m_code >> 8) & 0xFF;
}

/**
 * @brief Returns the value of the tower in case it is needed
 */
int TileCell::getTower(){

    return (m_code >> 16) & 0xFF;
}

/**
 *  @brief Returns de value of the module on phi in case it is needed
 */
int TileCell::getModule(){

    return (m_code >> 24) & 0xFF;
}

/**
  @brief This funcitions takes as inputs the partition, layer, tower and module and defines an integer to store
  this information in a memory friendly way. The integar has 4 bytes, one for each property.

  Byte 0 = Partition
  Byte 1 = Layer
  Byte 2 = Tower
  Byte 3 = Module

*/    
void TileCell::setCode(int partition, int layer, int tower, int module){

    m_code = 0;
    /// Part of the code for the partition (first 4 bits)
    m_code = partition;

    /// Part of the code for the layer (next 4 bits)
    layer = layer << 8;
    m_code = m_code | layer;

    /// Part of the code for the tower (next 4 bits)
    tower = tower << 16;
    m_code = m_code | tower;

    /// Part of the code for the tower (next 4 bits)
    module = module << 24;
    m_code = m_code | module;
}

/**
 *  @brief Return the energy of the given cell as the mean value of the distribution of the distribution generated with the energies which have been added.
 *
 *  It creates an histogram and gives its mean falue. It is onlye calculated once. If this functions is called before ending the energy scan the reult is useless.
 *  
 */
void TileCell::setEnergyNoise(){

    TH1D *tmp_histo = new TH1D("h","h",s_nBins,s_emin,s_emax);

    for (int bin = 0; bin <= s_nBins+1; bin++){
        tmp_histo->SetBinContent(bin,m_energy[bin]);
        tmp_histo->SetBinError(bin,sqrt(m_error[bin]));
    }

    m_meanEnergy = tmp_histo->GetMean();
    m_energyError = tmp_histo->GetMeanError();
    m_noise = tmp_histo->GetRMS();
    m_noiseError = tmp_histo->GetRMSError();

    m_isCalculated = true;
    delete tmp_histo;
}

/**
 * @brief Returns the entire vector of energy set at the momment
 */
vector<double> TileCell::getEnergyVec(){
    vector<double> v;
    for (int i = 0; i < s_nBins+2; i++) {
        v.push_back(m_energy[i]);
    }
    return v;
}

/**
 * @brief Returns the entire vector of errors set at the momment
 */
vector<double> TileCell::getNoiseVec(){
    vector<double> v;
    for (int i = 0; i < s_nBins+2; i++) {
        v.push_back(m_error[i]);
    }
    return v;
}
