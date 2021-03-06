#ifndef _CLEGAL_
#define _CLEGAL_

#include "Util.h"
#include "Row.h"
#include "Placement.h"
#include "GnuplotLivePlotter.h"
#include "GnuplotMatrixPlotter.h"
#include <math.h>
#include <stdio.h>

using namespace std;

class CLegal
{
    public:
        CLegal( Placement& placement );
       ~CLegal(void){}

        Placement& _placement;
        GnuplotLivePlotter gnuplotLivePlotter;

    //Solve the legalization, "true" indicates success
        bool solve();
    // helper functions add by Hao
        void plot();
    //Check the legalization result is overlap-free, within boundary, and on-site, "true" indicates legal
        bool check();
        void saveGlobalResult();
        void setLegalResult();
        double totalDisplacement();

        //**********************************
        //* Information of Macros or Cells *
        //**********************************
        vector<Module*> m_all_modules;
        vector<double> m_row_remain_width;
        vector<unsigned> m_cell_order; //Used as legalization order
        vector<unsigned> m_cell_order2; //Used as legalization order
        vector<Row> m_free_sites; //Record the unused sites
        vector<Row> m_left_free_sites;
        vector<CPoint> m_currentLocations;	//Record the current locations
        vector<CPoint> m_bestLocations;	//Record the best locations
        vector<Row> m_best_sites;       //Record the resulting sites of the best result

        vector<CPoint> m_bestLocations_left;
        vector<Row> m_best_sites_left;
        vector<CPoint> m_globalLocations;

        vector<int> m_macro_ids;
        vector<double> m_orig_widths;


        //Variables and functions used for 4-directional macro shifter
        vector<CPoint> m_macro_shifter_orig_positions;	//Record the original macro positions (left-bottom)
                                                        //before macro shifter
        vector<CPoint> m_macro_shifter_best_positions;	//Record the best macro positions (left-bottom)
                                                        //in 4-directional macro shifter

        vector<int> m_process_list;    //CalculateCellOrder() and CalculateNewLocation()
                                       //only act on modules in the m_process_list
        //********************************
        //*           parameter          *
        //********************************
        double m_max_module_height;
        double m_max_module_width;
        double m_average_cell_width;    //Average standard cell width
        double m_site_bottom, m_site_height; //Keep the bottom y coordinate of the sites and the
                                                                   //height of each site
        int m_unlegal_count;            //Number of unlegalized cells
        double m_chip_left_bound;	//Record the chip left boundary
        double m_chip_right_bound;	//Record the chip right boundary
};

#endif 

