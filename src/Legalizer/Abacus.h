#ifndef ABACUS_H
#define ABACUS_H

#include <limits>
#include "Legal.h"
using namespace std;

constexpr unsigned MAX_UNSIGNED = numeric_limits<unsigned>::max();

class Cluster {
    friend class Abacus;
    public:
        Cluster() 
            : _x(-1), _w(0), _e(0), _q(0),
            _rowId(MAX_UNSIGNED), _siteId(MAX_UNSIGNED),
            _prev(NULL), _next(NULL) {}
        ~Cluster() {}

        void addModule(Module* mod);
        void addCluster(Cluster& c);

        // get
        double x() const { return _x; }
        double w() const { return _w; }
        double e() const { return _e; }
        double q() const { return _q; }
        unsigned rowId() const { return _rowId; }
        unsigned siteId() const { return _siteId; }
        Cluster* prev() const { return _prev; }
        Cluster* next() const { return _next; }
        unsigned numModules() const { return _modules.size(); }
        Module*  firstModule() { return _modules.front(); }
        Module*  lastModule() { return _modules.back(); }

    private:    
        double _x;
        double _w;
        double _e; // weight (num of pins)
        double _q; // optimal = _q / _e
        unsigned _rowId;
        unsigned _siteId;
        Cluster *_prev, *_next;
        vector<Module*> _modules;
};

class Abacus {
    public:
        Abacus(CLegal& legalizer) : _legalizer(legalizer) {
            _clusters.resize(_legalizer.m_free_sites.size());
            _rowModules.resize(_legalizer.m_free_sites.size());
        }
        ~Abacus() {}

        void legal();

        // get 
        Cluster* lastCluster(unsigned rowId) const { 
            if (_clusters[rowId].size() > 0)
                return _clusters[rowId].back(); 
            else return NULL;
        }

    private:
        CLegal& _legalizer;
        vector<vector<Cluster*>> _clusters;
        vector<vector<Module*>>  _rowModules;

        // helper functions
        void placeRow(unsigned rowId, Module* mod, bool isTrial);
        void insertModule(Module* mod, unsigned rowId);
        void removeModule(Module* mod, unsigned rowId);
        void collapse(Cluster* c, double x_min, double x_max);
        void updateModulesLocation(Cluster* c);
        double computeCost(unsigned modId);
};

#endif