#include "Abacus.h"
#include <algorithm>
#include <cmath>
using namespace std;

#define all_modules                  _legalizer.m_all_modules
#define row_remain_width             _legalizer.m_row_remain_width
#define cell_order                   _legalizer.m_cell_order
#define cell_order2                  _legalizer.m_cell_order2
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

constexpr double MAX_DBL = numeric_limits<double>::max();
constexpr double MIN_DBL = numeric_limits<double>::lowest();

void Abacus::legal()
{
    vector<double> ori_r = row_remain_width;
    vector<CPoint> ori_g = globalLocations;
    
    vector<CPoint> best1;
    double disp1 = __legal();
    best1 = bestLocations;
    
    row_remain_width = ori_r;
    cell_order = cell_order2;
    bestLocations = globalLocations;
    for (auto& cluster : _clusters) cluster.clear();
    for (auto& rowModules : _rowModules) rowModules.clear();
    flip();
    
    double disp2 = __legal();
    
    flip();
   // globalLocations = ori_g;
   // double center = (chip_left_bound + chip_right_bound) / 2;
   // for (unsigned rowId = 0; rowId < _clusters.size(); ++rowId)
   // for (unsigned i = 0; i < _clusters[rowId].size(); ++i) {
   //     Cluster& c = _clusters[rowId][i];
   //     double x = 2 * center - c._x;
   //     for (unsigned j = 0; j < c._modules.size(); ++j) {
   //         unsigned modId = c._modules[j];
   //         Module* mod = all_modules[modId];
   //         x -= mod->width();
   //         bestLocations[modId].x = x;
   //         bestLocations[modId].y = free_sites[rowId].y();
   //     }
   // }
   // for (unsigned i = 0; i < free_sites.size(); ++i) {
   //     double newx = 2 * center - (free_sites[i].x() + free_sites[i].width());
   //     free_sites[i].setPosition(newx, free_sites[i].y());
   // }
   // cerr << disp1 << ' ' << disp2 << endl;
    if (disp1 < disp2) bestLocations = best1;
}

void Abacus::flip()
{
    double center = (chip_left_bound + chip_right_bound) / 2;
    for (unsigned i = 0; i < all_modules.size(); ++i) {
        Module* mod = all_modules[i];
        globalLocations[i].x = 2 * center - (globalLocations[i].x + mod->width());
        bestLocations[i].x = 2 * center - (bestLocations[i].x + mod->width());  
    }
    for (unsigned i = 0; i < free_sites.size(); ++i) {
        double newx = 2 * center - (free_sites[i].x() + free_sites[i].width());
        free_sites[i].setPosition(newx, free_sites[i].y());
    }
}

double Abacus::__legal()
{
    // already sorted in x order
    for (unsigned i = 0; i < all_modules.size(); ++i) {
        unsigned modId = cell_order[i];
        // ignore Terminals in this program
        if (all_modules[modId]->isFixed()) continue;
        ////////////////////////////////////
        double c_best = MAX_DBL;
        double cost = MAX_DBL;
        unsigned r_best = MAX_UNSIGNED;
        
        double curHeight = globalLocations[modId].y - site_bottom;
        int startR = floor(curHeight / site_height);
        for (int r = startR; r < (int)free_sites.size(); ++r) {
            if (all_modules[modId]->width() > row_remain_width[r]) continue;
            if (1.85 * fabs(globalLocations[modId].y - free_sites[r].y()) > c_best) break;
            //if (fabs(globalLocations[modId].y - free_sites[r].y()) > c_best) break;
            //if (pow(globalLocations[modId].y - free_sites[r].y(), 2) > c_best) break;
            insertModule(modId, r);
            cost = placeRow(r, true);
            if (cost < c_best) { 
                c_best = cost;
                r_best = r;
            }
            removeModule(r);
        }
        for (int r = startR - 1; r >= 0; --r) {
            if (all_modules[modId]->width() > row_remain_width[r]) continue;
            if (1.85 * fabs(globalLocations[modId].y - free_sites[r].y()) > c_best) break;
            //if (fabs(globalLocations[modId].y - free_sites[r].y()) > c_best) break;
            //if (pow(globalLocations[modId].y - free_sites[r].y(), 2) > c_best) break;
            insertModule(modId, r);
            cost = placeRow(r, true);
            if (cost < c_best) { 
                c_best = cost;
                r_best = r;
            }
            removeModule(r);
        }

        insertModule(modId, r_best);
        row_remain_width[r_best] -= all_modules[modId]->width();
        placeRow(r_best, false);
    }
    return totalDisplacement();
}

double Abacus::totalDisplacement()
{
    double totaldis = 0;
    for(unsigned i = 0; i < all_modules.size(); ++i) {
        totaldis += CPoint::Distance(globalLocations[i] , bestLocations[i]);
    }
    return totaldis;
}

void Abacus::addModule(Cluster& c, unsigned modId)
{
    c._modules.push_back(modId);
    Module* mod = all_modules[modId];
    c._e += 1;//mod->numPins();
    c._q += /*mod->numPins() */ (globalLocations[modId].x - c._w);
    c._w += mod->width();
}

void Abacus::addCluster(Cluster& c1, Cluster& c2)
{
    for (unsigned i = 0; i < c2._modules.size(); ++i) {
        c1._modules.push_back(c2._modules[i]);
    }
    c1._e += c2._e;
    c1._q += c2._q - c2._e * c1._w;
    c1._w += c2._w;
}

void Abacus::collapse(Cluster& c, double x_min, double x_max)
{
    c._x = c._q / c._e;
    if (c._x < x_min) c._x = x_min;
    if (c._x > x_max - c._w) c._x = x_max - c._w;
    
    if (c._prev != MAX_UNSIGNED) {
        Cluster& prev = _clusters[c._rowId][c._prev];
        if (prev._x + prev._w > c._x) {
            addCluster(prev, c);
            _clusters[c._rowId].pop_back();
            collapse(prev, x_min, x_max);
        }
    }
}

double Abacus::placeRow(unsigned rowId, bool isTrial)
{
    vector<Cluster> oldClusters = _clusters[rowId];
    unsigned modId = _rowModules[rowId].back();
    bool gen = (_clusters[rowId].empty() || 
                lastCluster(rowId)._x + lastCluster(rowId)._w <= globalLocations[modId].x);
    if (gen) {
        Cluster newC = Cluster();
        newC._rowId = rowId;
        newC._x = globalLocations[modId].x;
        if (!_clusters[rowId].empty())
            newC._prev = _clusters[rowId].size() - 1;
        addModule(newC, modId);
        _clusters[rowId].push_back(newC);
    }
    else {
        Cluster& c = lastCluster(rowId);
        addModule(c, modId);
        double x_min = free_sites[rowId].x();
        double x_max = x_min + free_sites[rowId].width();
        collapse(c, x_min, x_max);
    }
    double cost = 0, placeX = 0, increaseDP = 0;
    for (unsigned i = 0; i < _clusters[rowId].size(); ++i) {
        Cluster& c = _clusters[rowId][i];
        double x = c._x;
        for (unsigned j = 0; j < c._modules.size(); ++j) {
            unsigned modId = c._modules[j];
            Module* mod = all_modules[modId];
            if (!isTrial) {
                bestLocations[modId].x = x;
                bestLocations[modId].y = free_sites[rowId].y();
            }
            else {
                increaseDP += fabs(x - bestLocations[modId].x);
            }
            x += mod->width();
        }
        placeX = x;
    }
    Module* mod = all_modules[_clusters[rowId].back()._modules.back()];
    placeX -= mod->width();
    // define cost
    cost =  (fabs(placeX - globalLocations[modId].x) +
             1.85 * fabs(free_sites[rowId].y() - globalLocations[modId].y)) +
             increaseDP;
    //cost = fabs(placeX - globalLocations[modId].x) +
    //       fabs(free_sites[rowId].y() - globalLocations[modId].y);
    //cost = pow(placeX - globalLocations[modId].x, 2) +
    //        pow(free_sites[rowId].y() - globalLocations[modId].y, 2);
    //cost = sqrt(cost);
    if (isTrial) _clusters[rowId] = oldClusters;
    return cost;
}

void Abacus::insertModule(unsigned modId, unsigned rowId)
{
    _rowModules[rowId].push_back(modId);
}

void Abacus::removeModule(unsigned rowId)
{
    _rowModules[rowId].pop_back();
}
