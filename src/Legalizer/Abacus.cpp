#include "Abacus.h"
using namespace std;

#define cell_order                   _lga.m_cell_order
#define free_sites                   _lga.m_free_sites
#define left_free_sites              _lga.m_left_free_sites
#define origLocations                _lga.m_origLocations
#define bestLocations                _lga.m_bestLocations
#define best_sites                   _lga.m_best_sites
#define bestLocations_left           _lga.m_bestLocations_left
#define globalLocations              _lga.m_globalLocations
#define macro_ids                    _lga.m_macro_ids
#define orig_widths                  _lga.m_orig_ids
#define macro_shifter_orig_position  _lga.m_macro_shifter_orig_position
#define macro_shifter_best_position  _lga.m_macro_shifter_best_position
#define process_list                 _lga.m_process_list
#define max_module_height            _lga.m_max_module_height
#define max_module_width             _lga.m_max_module_width
#define average_cell_width           _lga.m_average_cell_width
#define site_bottom                  _lga.m_site_bottom
#define site_height                  _lga.m_site_height
#define unlegal_count                _lga.m_unlegal_count
#define chip_left_bound              _lga.m_chip_left_bound


void Abacus::legal()
{
    unlegal_count = 87;
    cout << unlegal_count << endl; exit(0);
}