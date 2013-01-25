#ifndef TILECELL_H
#define TILECELL_H

#include <vector>
#include <string>
#include <map>



// ---------------------------------------------------------
/**
 *  @class: TileCell
 @brief: This class store basic information of Tile cells which
 is stored in TileD3Pds.


 @author: J.P. Araque <jp.araque@cern.ch>
 @version 1.0
 */    
// ---------------------------------------------------------
class TileCell {
    public:
        TileCell (int, int, int, int, int, int);
        virtual ~TileCell ();

        enum KeyString{ P,L,T,M,
            PL,
            PLT,
            PLTM};//!< This enum is used to choos which key we want to know for each cell

        static void setHistoInfo(int bins, double xmin, double xmax);
        static void printEnergyInfo();
        void setEnergy(double energy, double weight);

        std::string getKey(KeyString set = PL, std::string separator = "_", bool _end = true);
        double getEnergy();
        double getEnergyError();
        double getNoise();
        double getNoiseError();
        bool isGood();
        int getCode();
        int getPartition();
        int getLayer();
        int getTower();
        int getModule();
        std::vector<double> getEnergyVec();
        std::vector<double> getNoiseVec();


    protected:
        void setCode(int, int, int, int);
        void setEnergyNoise();


    private:
        bool m_isMC;//!< Tells you if this Cell is used to store MC or data information.
        int m_code;
        bool m_isGoodCell;//!< Tells if the Cell is active or not based on D3PD information.
        static bool s_isSet; //!< Tells if the infromation of the histogram has been set.
        bool m_isCalculated; //!< Tells if the mean energy and error of the energy distribution have been alredy calculated.

        double  *m_energy;//!< A vector to store the energy deposited in the TileCell 
        //!< in MeV separed in future bins for future histogram. 
        double  *m_error;//!< A vector to store the error for each bean of energy
        static double s_maxEnergy; //!< It store the max energy that have been added to the cell just 
        //!< in cas it exceed the expected limits.
        static double s_minEnergy; //!< It store the max energy that have been added to the cell just
        //!< in cas it exceed the expected limits.
        double GeV;//!< Just to make the changed between MeV and GeV.
        static int s_nBins; //!< Number of bins for the Energy histogram
        static double s_emin;//!< Low limit for the histogram.
        static double s_emax;//!< Max limit for the histogram.
        double m_meanEnergy; //!< Mean energy of the energy distribution
        double m_energyError;//!< Error in the energy
        double m_noise;//!< Noise asociated with the energy distribution
        double m_noiseError; //!< Error in the noise

};

typedef std::map<int,TileCell*> TileCellContainer;//!< Map container to use it when definining lots of Cells in noise studies
#endif
