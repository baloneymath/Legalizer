#include "Legal.h"
#include "arghandler.h"
#include "GnuplotLivePlotter.h"
#include "GnuplotMatrixPlotter.h"
#include "Abacus.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <set>
#include <cassert>
using namespace std;

bool CLegal::solve()
{
    saveGlobalResult();
    // TODO: edit your code HERE
    // Note:
    //      1. You should save your legal solution into m_bestLocations, and call setLegalResult() tp upload it into Placement.
    //      2. Run check() to make sure that the solution is lega.
    //      3. Feel free to add any class, function, and/or data member.
    // Good luck!
    
    Abacus abacus(*this);
    abacus.legal();

    setLegalResult();
    if( check() )
    {
        cout << "total displacement: " << totalDisplacement() << endl;
        return true;
    }
    else
        return false;
}

bool sortModule( Module* a, Module* b)
{
    return a->x() < b->x();
}

bool CLegal::check()
{
    cout << "start check" << endl;
    int notInSite=0;
    int notInRow=0;
    int overLap=0;

    ///////////////////////////////////////////////////////
    //1.check all standard cell are on row and in the core region
    //////////////////////////////////////////////////////////
    for(unsigned int i=0; i<_placement.numModules(); ++i)
    {
        Module& module = _placement.module(i);
        if(module.isFixed()) continue;
        double curX = module.x();
        double curY = module.y();

        double res = ( curY - m_site_bottom ) / _placement.getRowHeight();
        //cout << curY << " " << res << endl;
        int ires = (int) res;
        if( (m_site_bottom + _placement.getRowHeight() * ires) != curY )
        {
                cerr<<"\nWarning: cell:"<<i<<" is not on row!!";
                ++notInRow;
        }
        if( (curY<_placement.boundaryBottom()) || (curX<_placement.boundaryLeft())||
                ( (curX+module.width())>_placement.boundaryRight()) ||
                ( (curY+module.height())>_placement.boundaryTop()) )
        {
            cerr<<"\nWarning: cell:"<<i<<" is not in the core!!";
            ++notInSite;
        }
    }

    ///////////////////////////////////////////
    //2. row-based overlapping checking
    ///////////////////////////////////////////

    Rectangle chip = _placement.rectangleChip();
    const double &rowHeight = _placement.getRowHeight();
    unsigned numRows = _placement.numRows();
    vector< vector<Module*> > rowModules( numRows, vector<Module*>( 0 ) );
    for(unsigned int i=0; i<_placement.numModules(); ++i)
    {
        Module& module = _placement.module(i);
        double curY = m_bestLocations[i].y;

        if( module.area() == 0 ) continue;
        if( module.isFixed() ) continue;

        double yLow = curY - chip.bottom(), yHigh= curY + module.height() - chip.bottom();
        size_t low = floor( yLow / rowHeight ), high = floor(yHigh / rowHeight);
        if( fabs( yHigh - rowHeight * floor(yHigh / rowHeight) ) < 0.01 ){ --high; }

        for( size_t i = low; i <= high; ++i ){ rowModules[ i ].push_back( &module ); }
    }
    for( size_t i = 0; i < numRows; ++i )
    {
        vector<Module*> &modules = rowModules[i];
        sort(modules.begin(), modules.end(), sortModule);
        if( modules.size() < 1 ) { continue; }
        for( size_t j = 0; j < modules.size() - 1; ++j ){
            Module &mod = *modules[ j ];
            if(mod.isFixed()) continue;
            size_t nextId = j+1;
            while( mod.x() + mod.width() - modules[ nextId ]->x() > 0.01 ){
                Module &modNext = *modules[ nextId ];
                if( mod.x() + mod.width() - modules[ nextId ]->x() > 0.01 ){
                    ++overLap;
                    cout << mod.name() << " overlap with " << modNext.name() << endl;
                }
                ++nextId; if( nextId == modules.size() ) { break; }
            }
        }
    }

    cout << endl <<
            "  # row error: "<<notInRow<<
            "\n  # site error: "<<notInSite<<
            "\n  # overlap error: "<<overLap<< endl;
    //cout << "end of check" << endl;

    if( notInRow!=0 || notInSite!=0 || overLap!=0 )
    {
        cout <<"Check failed!!" << endl;
        return false;
    }
    else
    {
        cout <<"Check success!!" << endl;
        return true;
    }
}

double CLegal::totalDisplacement()
{
    double totaldis = 0;
    for( unsigned  moduleId = 0 ; moduleId < _placement.numModules() ; moduleId++ )
    {
        Module& curModule = _placement.module(moduleId);
        double x = curModule.x();
        double y = curModule.y();

        totaldis += CPoint::Distance(m_globalLocations[moduleId] , CPoint( x, y ));
    }
    return totaldis;
}

CLegal::CLegal( Placement& placement  ) :
    _placement( placement )
{

    //Compute average cell width
    int cell_count = 0;
    double total_width = 0;
    //double max_height = 0.0;
    m_max_module_height = 0.0;
    m_max_module_width = 0.0;
    for( unsigned  moduleId = 0 ; moduleId < placement.numModules() ; moduleId++ )
    {
        Module& curModule = placement.module(moduleId);
        Module* pCur = &curModule;
        m_all_modules.push_back(pCur);

        m_max_module_height = max( m_max_module_height, curModule.height() );
        m_max_module_width = max( m_max_module_width, curModule.width() );
    //Do not include fixed cells and macros
        if( curModule.isFixed() || curModule.height() > placement.getRowHeight() )
        continue;

        cell_count++;
        total_width += curModule.width();
    }

    m_average_cell_width = total_width / cell_count;

    m_free_sites = placement.m_sites;
    m_site_bottom = m_free_sites.front().y();
    m_site_height = m_free_sites.front().height();
    m_row_remain_width.resize(m_free_sites.size());
    for (unsigned i = 0; i < m_free_sites.size(); ++i) {
        m_row_remain_width[i] = m_free_sites[i].width();
    }
 
    //initalize m_origLocations and m_bestLocations
    m_cell_order.resize( placement.numModules() );
    m_bestLocations.resize( placement.numModules() );
    m_globalLocations.resize( placement.numModules() );
    m_chip_left_bound = placement.rectangleChip().left();
    m_chip_right_bound = placement.rectangleChip().right();
    
    vector<pair<Module*, unsigned>> tmp;
    for (unsigned i = 0; i < m_cell_order.size(); ++i) {
        tmp.push_back(make_pair(m_all_modules[i], i));
    }
    auto sortByCX = [] (pair<Module*, unsigned> p1, pair<Module*, unsigned> p2) {
        return p1.first->centerX() < p2.first->centerX();
    };
    auto sortByX = [] (pair<Module*, unsigned> p1, pair<Module*, unsigned> p2) {
        return p1.first->x() < p2.first->x();
    };
    sort(tmp.begin(), tmp.end(), sortByCX);
    for (unsigned i = 0; i < m_cell_order.size(); ++i) {
        m_cell_order[i] = tmp[i].second;
    }
    m_cell_order2.resize(m_cell_order.size());
    sort(tmp.begin(), tmp.end(), sortByX);
    for (unsigned i = 0; i < m_cell_order2.size(); ++i) {
        m_cell_order2[i] = tmp[i].second;
    }
    reverse(m_cell_order2.begin(), m_cell_order2.end());
}

void CLegal::saveGlobalResult()
{
    for (unsigned moduleId = 0; moduleId < _placement.numModules(); moduleId++)
    {
        Module &curModule = _placement.module(moduleId);
        double x = curModule.x();
        double y = curModule.y();

        m_globalLocations[moduleId] = CPoint( x, y );
        m_bestLocations[moduleId] = CPoint( x, y );
    }
}

void CLegal::setLegalResult()
{
    for (unsigned moduleId = 0; moduleId < _placement.numModules(); moduleId++)
    {
        Module &curModule = _placement.module(moduleId);
        curModule.setPosition(m_bestLocations[moduleId].x,
                              m_bestLocations[moduleId].y);
    }
}

void CLegal::plot()
{
    for (unsigned i = 0; i < _placement.numModules(); ++i) {
        Module& curModule = _placement.module(i);
        gnuplotLivePlotter.addRectangle(curModule.rectangle());
    }
    gnuplotLivePlotter.show();
    gnuplotLivePlotter.clearObjects();
}
