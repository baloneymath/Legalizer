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
            _prev(MAX_UNSIGNED) {}
        ~Cluster() {}

        // get
        double x() const { return _x; }
        double w() const { return _w; }
        double e() const { return _e; }
        double q() const { return _q; }
        unsigned rowId() const { return _rowId; }
        unsigned siteId() const { return _siteId; }

    private:    
        double _x;
        double _w;
        double _e; // weight (num of pins)
        double _q; // optimal = _q / _e
        unsigned _rowId;
        unsigned _siteId;
        unsigned _prev; // prev row clusterId
        vector<unsigned> _modules;
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
        Cluster& lastCluster(unsigned rowId) { 
            return _clusters[rowId].back();
        }

    private:
        CLegal& _legalizer;
        vector<vector<Cluster>>  _clusters;
        vector<vector<unsigned>>  _rowModules;

        // helper functions
        double placeRow(unsigned rowId, bool isTrial);
        void addModule(Cluster& c, unsigned modId);
        void addCluster(Cluster& c1, Cluster& c2);
        void insertModule(unsigned modId, unsigned rowId);
        void removeModule(unsigned modId, unsigned rowId);
        void collapse(Cluster& c, double x_min, double x_max);
};

#endif