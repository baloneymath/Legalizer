#include "Abacus.h"
using namespace std;

#define cell_order                   _legalizer.m_cell_order
#define free_sites                   _legalizer.m_free_sites
#define left_free_sites              _legalizer.m_left_free_sites
#define origLocations                _legalizer.m_origLocations
#define bestLocations                _legalizer.m_bestLocations
#define best_sites                   _legalizer.m_best_sites
#define bestLocations_left           _legalizer.m_bestLocations_left
#define globalLocations              _legalizer.m_globalLocations
#define macro_ids                    _legalizer.m_macro_ids
#define orig_widths                  _legalizer.m_orig_ids
#define macro_shifter_orig_position  _legalizer.m_macro_shifter_orig_position
#define macro_shifter_best_position  _legalizer.m_macro_shifter_best_position
#define process_list                 _legalizer.m_process_list
#define max_module_height            _legalizer.m_max_module_height
#define max_module_width             _legalizer.m_max_module_width
#define average_cell_width           _legalizer.m_average_cell_width
#define site_bottom                  _legalizer.m_site_bottom
#define site_height                  _legalizer.m_site_height
#define unlegal_count                _legalizer.m_unlegal_count
#define chip_left_bound              _legalizer.m_chip_left_bound
#define chip_right_bound             _legalizer.m_chip_right_bound


void Abacus::legal()
{
}

void Abacus::placeRow(unsigned rowId, Module* mod)
{
    for (unsigned i = 0; i < _rowModules[rowId].size(); ++i) {
        Cluster* c = lastCluster(rowId);
        if (i == 0) { // first module, c == NULL
            Cluster* newC = new Cluster();
            newC->_x = mod->x();
            newC->addModule(mod);
        }
        else if (c != NULL && c->_x + c->_w <= mod->x()) {
            c->addModule(mod);
            double x_min = free_sites[rowId].x();
            double x_max = x_min + free_sites[rowId].width();
            collapse(c, x_min, x_max);
        }
    }
    for (unsigned j = 0; j < _clusters[rowId].size(); ++j) {
        Cluster* cc = _clusters[rowId][j];
        double x = cc->_x;
        for (unsigned k = 0; k < cc->numModules(); ++k) {
            cc->_modules[k]->setPosition(x, free_sites[rowId].y());
            x += cc->_modules[k]->width();
        }
    }
}

void Cluster::addModule(Module* mod)
{
    _modules.push_back(mod);
    _e += mod->numPins();
    _q += mod->numPins() * (mod->x() - _w);
    _w += mod->width();
}

void Cluster::addCluster(Cluster& c)
{
    for (unsigned i = 0; i < c.numModules(); ++i) {
        _modules.push_back(c.firstModule());
    }
    _e += c._e;
    _q += c._q - (c._e * c._w);
    _w += c._w;
}

void Abacus::collapse(Cluster* c, double x_min, double x_max)
{
    vector<Cluster*>& clusters = _clusters[c->_rowId];
    c->_x = (c->_q / c->_e);
    if (c->_x < x_min) c->_x = x_min;
    if (c->_x > x_max - c->_w) c->_x = x_max - c->_w;
    Cluster* prev = c->_prev;
    if (prev != NULL) {
        if (prev->_x + prev->_w > c->_x) {
            prev->addCluster(*c);
            prev->_next = c->_next;
            c->_next->_prev = prev;
            clusters.erase(clusters.begin() + c->_siteId);
            delete c;
            collapse(prev, x_min, x_max);
        }
    }
}