#ifndef ABACUS_H
#define ABACUS_H

#include <limits>
#include "Legal.h"
using namespace std;

constexpr double MAX_DOUBlE = numeric_limits<double>::max();
constexpr unsigned MAX_UNSIGNED = numeric_limits<unsigned>::max();

struct Cluster {
    Cluster() 
        : _x(-1), _width(0), _edgeWeight(0), _q(0),
          _rowId(MAX_UNSIGNED), _siteId(MAX_UNSIGNED),
          _prev(NULL), _next(NULL) {}
    ~Cluster() {}

    void addModule(Module* mod) { _modules.push_back(mod); }
    void addCluster(Cluster& c) {
        for (unsigned i = 0; i < c._modules.size(); ++i) {
            Module* mod = c._modules[i];
            _modules.push_back(mod);
            _width += mod->width();
        }
        _next = c._next;
    }
    
    double _x;
    double _width;
    double _edgeWeight;
    double _q; // optimal = _q / _edgeWeight
    unsigned _rowId;
    unsigned _siteId;
    Cluster *_prev, *_next;
    vector<Module*> _modules;
};

class Abacus {
    public:
        Abacus(CLegal& lga)
            : _lga(lga) {}
        ~Abacus() {}

        void legal();

    private:
        CLegal& _lga;
        vector<Cluster*> _lastClusters;
};

#endif